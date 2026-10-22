#!/usr/bin/env python
"""
Summaries of the Linux runs for report/data/linux_*.json.

    summarise.py validation DIR --label NAME --output FILE
        condense the validation_*.json written by
        scripts/validate_against_lalsimulation.py and
        scripts/validate_nrsur7dq4.py into DIR (per model: draws compared,
        bit-identical draws, largest and median errors)
    summarise.py bench DIR --label NAME --output FILE
        timings written by scripts/linux/bench.sh into DIR
    summarise.py environment
        print the environment record (platform, CPU, toolchain, packages)

Every summary records the environment it was produced in. If FILE exists,
the entry for NAME is replaced and other labels are kept.
"""
import argparse
import json
import os
import platform
import re
import subprocess
from pathlib import Path

import numpy as np

PACKAGES = ("python-lalsimulation", "liblalsimulation", "liblal", "gsl",
            "libopenblas", "libblas", "numpy", "python")


def conda_packages():
    prefix = Path(os.environ.get("CONDA_PREFIX", "/opt/conda/envs/lal"))
    found = dict()
    for meta in sorted((prefix / "conda-meta").glob("*.json")):
        record = json.loads(meta.read_text())
        if record["name"] in PACKAGES:
            found[record["name"]] = f"{record['version']}-{record['build']}"
    return found


def cpu():
    info = dict(machine=platform.machine())
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return info
    flags = re.search(r"^(flags|Features)\s*:\s*(.*)$", text, re.M)
    if flags:
        names = set(flags.group(2).split())
        wanted = ("sse4_1", "avx", "avx2", "fma", "avx512f", "asimd", "fphp", "sve")
        info["features"] = sorted(names & set(wanted))
    model = re.search(r"^model name\s*:\s*(.*)$", text, re.M)
    if model:
        info["model"] = model.group(1)
    info["n_cpu"] = os.cpu_count()
    return info


def environment():
    def output(*command):
        try:
            return subprocess.run(command, capture_output=True, text=True,
                                  check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    return dict(
        platform=os.environ.get("WAVEFORMRS_PLATFORM", f"linux-{platform.machine()}"),
        kernel=platform.release(),
        glibc=" ".join(platform.libc_ver()),
        cpu=cpu(),
        rustc=output("rustc", "--version"),
        rustflags=os.environ.get("RUSTFLAGS", ""),
        packages=conda_packages(),
        omp_num_threads=os.environ.get("OMP_NUM_THREADS"),
    )


def summarise_model(path):
    data = json.loads(path.read_text())
    records = data["records"]
    errors = np.array([max(r["max_rel_err_plus"], r["max_rel_err_cross"]) for r in records])
    mismatches = np.array([max(r["mismatch_plus"], r["mismatch_cross"]) for r in records])
    # NRSur7dq4 records also time a fresh evaluation (time_rust reuses a
    # workspace that caches the modes of the same binary)
    speedup = np.array([r["time_lal"] / r.get("time_rust_fresh", r["time_rust"]) for r in records])
    summary = dict(
        n_draws=data["n_draws"],
        n_compared=data["n_compared"],
        n_bit_identical=int(np.sum(errors == 0)),
        max_relative_error=float(errors.max()),
        median_relative_error=float(np.median(errors)),
        p90_relative_error=float(np.quantile(errors, 0.9)),
        max_mismatch=float(mismatches.max()),
        median_mismatch=float(np.median(mismatches)),
        median_speedup=float(np.median(speedup)),
        failures=data["failures"],
    )
    for key in ("n_lal_errno_fallback", "n_epoch_mismatch", "n_length_mismatch"):
        if key in data:
            summary[key] = data[key]
    # the draw with the largest error, to find it again
    worst = records[int(np.argmax(errors))]
    summary["worst_draw"] = {
        key: worst[key] for key in ("mass_1", "mass_2", "spin_1", "spin_2", "f_min", "f_ref")
    }
    return data["approximant"], summary


def validation(args):
    models = dict()
    for path in sorted(Path(args.directory).glob("validation_*.json")):
        name, summary = summarise_model(path)
        models[name] = summary
    output = Path(args.output)
    combined = json.loads(output.read_text()) if output.exists() else dict()
    combined.setdefault("description", (
        "waveformrs against the lalsimulation of the same Linux container "
        "(conda-forge build), random draws of scripts/validate_against_lalsimulation.py "
        "and scripts/validate_nrsur7dq4.py with OMP_NUM_THREADS=1. Relative errors are "
        "max|h_rs - h_lal| / max|h_lal| over both polarisations; bit-identical draws "
        "have zero error. Speed-ups are only meaningful on native (arm64) runs."))
    combined[args.label] = dict(environment=environment(), models=models)
    output.write_text(json.dumps(combined, indent=1) + "\n")
    print(f"{'model':<22}{'identical':>10}{'max rel':>11}{'median':>11}{'max mm':>11}")
    for name, s in models.items():
        print(f"{name:<22}{s['n_bit_identical']:>5}/{s['n_compared']:<4}"
              f"{s['max_relative_error']:11.2e}{s['median_relative_error']:11.2e}"
              f"{s['max_mismatch']:11.2e}")


GRID_LINE = re.compile(
    r"^(?P<name>.+?)\s+per-frequency\s+(?P<old>[\d.]+) us\s+grid\s+(?P<new>[\d.]+) us"
    r"\s+x[\d.]+\s+max rel diff (?P<diff>\S+)$")


def bench(args):
    directory = Path(args.directory)
    entry = dict(environment=environment())
    timings = directory / "timings.json"
    if timings.exists():
        data = json.loads(timings.read_text())
        entry["single"] = {
            name: dict(rust_us=t * 1e6, lal_us=data["lal_sequence"][name] * 1e6,
                       speedup=data["lal_sequence"][name] / t)
            if name in data["lal_sequence"] else dict(rust_us=t * 1e6)
            for name, t in data["rust"].items()
        }
        entry["lal_uniform_us"] = {name: t * 1e6 for name, t in data["lal_uniform"].items()}
    elif (directory / "benchmark.json").exists():
        data = json.loads((directory / "benchmark.json").read_text())
        entry["single"] = {name: dict(rust_us=v["best"] * 1e6) for name, v in data.items()}
    serial = directory / "serial_timings.json"
    if serial.exists():
        data = json.loads(serial.read_text())
        entry["serial"] = dict(
            n_draws=data["n_draws"], n_frequencies=data["n_frequencies"],
            lal_python_call_overhead_us=data["lal_python_call_overhead"] * 1e6,
            models={
                name: dict(
                    rust_us=r["seconds"] * 1e6,
                    rust_grid_us=r["grid_seconds"] * 1e6 if r.get("grid_seconds") else None,
                    lal_us=data["lal"][name]["seconds"] * 1e6,
                    speedup=data["lal"][name]["seconds"] / r["seconds"],
                )
                for name, r in data["rust"].items() if name in data["lal"]
            })
    grid = directory / "grid.txt"
    if grid.exists():
        entry["grid"] = {
            m["name"]: dict(per_frequency_us=float(m["old"]), grid_us=float(m["new"]),
                            max_rel_diff=float(m["diff"]))
            for m in map(GRID_LINE.match, grid.read_text().splitlines()) if m
        }
    batch = directory / "batch.txt"
    if batch.exists():
        lines = batch.read_text().splitlines()
        threads = [int(t) for t in re.findall(r"(\d+) threads?\b", next(
            (line for line in lines if line.startswith("model")), ""))]
        rows = dict()
        for line in lines:
            m = re.match(r"^(?P<name>.+?)\s+(?P<n>\d+)\s+(?P<cells>\d+\s+\(.*)$", line)
            if m:
                rates = [float(r) for r in re.findall(r"(\d+)\s+\(\s*[\d.]+x\)", m["cells"])]
                rows[m["name"]] = dict(n=int(m["n"]), waveforms_per_second=dict(
                    zip(map(str, threads), rates)))
        entry["batch"] = dict(threads=threads, models=rows)
    output = Path(args.output)
    combined = json.loads(output.read_text()) if output.exists() else dict()
    combined.setdefault("description", (
        "Single-thread timings in microseconds. single: examples/benchmark.rs (construction "
        "+ polarizations, 5921 frequencies, best of 10 batches) against "
        "SimInspiralChooseFDWaveformSequence (report/lal_benchmark.py); serial: random "
        "binaries one after another on a fixed grid (report/serial_benchmark.py); grid: "
        "examples/grid_benchmark.rs; batch: examples/batch_benchmark.rs, waveforms per second "
        "on 1 to all threads of the machine. Emulated (linux/amd64 on Apple silicon) timings are "
        "not meaningful as absolute numbers."))
    combined[args.label] = entry
    output.write_text(json.dumps(combined, indent=1) + "\n")
    for name, s in entry.get("single", {}).items():
        lal = f"{s['lal_us']:10.1f}{s['speedup']:8.2f}" if "lal_us" in s else ""
        print(f"{name:<34}{s['rust_us']:10.1f}{lal}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validation", "bench"):
        p = sub.add_parser(name)
        p.add_argument("directory")
        p.add_argument("--label", required=True)
        p.add_argument("--output", required=True)
    sub.add_parser("environment")
    args = parser.parse_args()
    if args.command == "validation":
        validation(args)
    elif args.command == "bench":
        bench(args)
    else:
        print(json.dumps(environment(), indent=1))


if __name__ == "__main__":
    main()
