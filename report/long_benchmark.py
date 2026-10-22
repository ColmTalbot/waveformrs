#!/usr/bin/env python
"""
Timings on long frequency grids: low-mass binaries whose signal lasts up to
~150 s from 20 Hz, on grids with Delta f = 1/128 Hz (up to ~520k
frequencies), single thread.

For every frequency-domain model except NRSur7dq4 this times

- lalsimulation: ``SimInspiralChooseFDWaveformSequence`` on the grid
  (multibanding is never used for sequences) and, for IMRPhenomXP, XHM and
  XPHM, ``SimInspiralChooseFDWaveform`` on the uniform grid of the same band
  with multibanding off and (XHM, XPHM) with its default multibanding;
- waveformrs (``examples/long_benchmark.rs``): the per-frequency methods
  (construction plus polarizations), the grid paths (``*_on_grid`` on a
  prepared ``FrequencyGrid`` with reused buffers) and ``polarizations_uniform``
  for IMRPhenomXP/XPHM.

Other processes load the machine, so the timings are interleaved: each round
runs every LAL call once and then the Rust benchmark (which alternates its
own paths), and the best of the rounds is kept. With ``--memory`` the peak
resident memory of one evaluation of each case is also measured, for both
codes, in fresh processes.

    OMP_NUM_THREADS=1 python report/long_benchmark.py [--rounds 5] [--memory]

Writes ``report/data/long_timings.json``.
"""
import argparse
import datetime
import json
import os
import platform
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import lal
import lalsimulation as ls
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "report" / "data" / "long_timings.json"
NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0.0, PhenomXPHMThresholdMband=0.0)

# 3 Msun: ~150 s from 20 Hz (XLALSimInspiralChirpTimeBound)
ALIGNED = dict(mass_1=1.6, mass_2=1.4, spin_1=[0.0, 0.0, 0.05], spin_2=[0.0, 0.0, 0.02],
               luminosity_distance=150.0, inclination=0.6, phi_ref=1.3, f_ref=20.0)
# 3.1 Msun, ~130 s from 20 Hz, in-plane spins
PRECESSING = dict(mass_1=1.8, mass_2=1.3, spin_1=[0.3, 0.1, 0.05], spin_2=[-0.1, 0.2, 0.02],
                  luminosity_distance=200.0, inclination=1.1, phi_ref=0.4, f_ref=20.0)
GRIDS = {
    "1 kHz": dict(f_min=20.0, f_max=1024.0, delta_f=1 / 128),
    "2 kHz": dict(f_min=20.0, f_max=2048.0, delta_f=1 / 128),
    "4 kHz": dict(f_min=20.0, f_max=4096.0, delta_f=1 / 128),
}
MODELS = [
    ("TaylorF2", None), ("IMRPhenomD", None), ("IMRPhenomXAS", None), ("IMRPhenomXHM", None),
    ("IMRPhenomXP", 102), ("IMRPhenomXP", 223), ("IMRPhenomXP", 320),
    ("IMRPhenomXPHM", 102), ("IMRPhenomXPHM", 223), ("IMRPhenomXPHM", 320),
]
# The models that LAL can multiband on uniform grids, and those with a
# uniform-grid path timed for both codes.
MULTIBANDING = {"IMRPhenomXHM", "IMRPhenomXPHM"}
UNIFORM = {"IMRPhenomXP", "IMRPhenomXHM", "IMRPhenomXPHM"}


def cases(grids):
    out = list()
    for grid_name in grids:
        for model, version in MODELS:
            binary = ALIGNED if version is None else PRECESSING
            label = model if version is None else f"{model} ({version})"
            out.append(dict(binary, label=f"{label}, {grid_name}", model=model,
                            version=version, grid=grid_name, **GRIDS[grid_name]))
    return out


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in flags.items():
        getattr(ls, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def flags(case, multibanding=False):
    out = dict() if multibanding else dict(NO_MULTIBANDING)
    if case["version"] is not None:
        out["PhenomXPrecVersion"] = case["version"]
    return out


def frequencies(case):
    n = int(round((case["f_max"] - case["f_min"]) / case["delta_f"])) + 1
    return case["f_min"] + np.arange(n) * case["delta_f"]


class LalCalls:
    """The LAL calls of one case, with their arguments prepared once."""

    def __init__(self, case):
        self.case = case
        self.approximant = ls.GetApproximantFromString(case["model"])
        freqs = frequencies(case)
        self.sequence = lal.CreateREAL8Vector(len(freqs))
        self.sequence.data = freqs
        self.masses = (case["mass_1"] * lal.MSUN_SI, case["mass_2"] * lal.MSUN_SI)
        self.distance = case["luminosity_distance"] * 1e6 * lal.PC_SI
        self.params = lal_dict(flags(case))
        self.params_mb = lal_dict(flags(case, multibanding=True))

    def run_sequence(self):
        c = self.case
        return ls.SimInspiralChooseFDWaveformSequence(
            c["phi_ref"], *self.masses, *c["spin_1"], *c["spin_2"], c["f_ref"],
            self.distance, c["inclination"], self.params, self.approximant, self.sequence)

    def run_uniform(self, multibanding=False):
        c = self.case
        return ls.SimInspiralChooseFDWaveform(
            *self.masses, *c["spin_1"], *c["spin_2"], self.distance, c["inclination"],
            c["phi_ref"], 0.0, 0.0, 0.0, c["delta_f"], c["f_min"], c["f_max"], c["f_ref"],
            self.params_mb if multibanding else self.params, self.approximant)

    def paths(self):
        out = dict(sequence=self.run_sequence)
        if self.case["model"] in UNIFORM:
            out["uniform"] = self.run_uniform
        if self.case["model"] in MULTIBANDING:
            out["uniform_multibanding"] = lambda: self.run_uniform(multibanding=True)
        return out


def seconds(call):
    start = time.perf_counter()
    call()
    return time.perf_counter() - start


def run_rust(case_list, repeats, path):
    path.write_text(json.dumps(dict(repeats=repeats, cases=[
        {k: v for k, v in case.items() if k != "grid"} for case in case_list])))
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--release", "--example", "long_benchmark", "--", str(path)],
        cwd=ROOT, capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def max_rss_bytes(command):
    """Peak resident memory of `command` run in a fresh wrapper process."""
    wrapper = ("import resource, subprocess, sys; subprocess.run(sys.argv[1:], check=True, "
               "stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); "
               "print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)")
    out = subprocess.run([sys.executable, "-c", wrapper, *command], capture_output=True,
                         text=True, check=True, cwd=ROOT)
    value = int(out.stdout.strip())
    # bytes on macOS, kilobytes on Linux
    return value if platform.system() == "Darwin" else 1024 * value


LAL_MEMORY = """
import json, sys
sys.path.insert(0, {report!r})
from long_benchmark import LalCalls
calls = LalCalls(json.loads({case!r}))
calls.paths()[{path!r}]()
"""


def measure_memory(case_list, tmp):
    """Peak RSS of one evaluation of each path (both codes), and of a
    process that only imports lalsimulation or loads the Rust binary."""
    binary = ROOT / "target" / "release" / "examples" / "long_benchmark"
    out = dict()
    baseline = tmp / "baseline.json"
    baseline.write_text(json.dumps(dict(repeats=0, cases=[])))
    out["baseline"] = dict(
        rust=max_rss_bytes([str(binary), str(baseline)]),
        lal=max_rss_bytes([sys.executable, "-c", "import lal, lalsimulation"]))
    for i, case in enumerate(case_list):
        path = tmp / f"memory_{i}.json"
        path.write_text(json.dumps(dict(repeats=0, cases=[
            {k: v for k, v in case.items() if k != "grid"}])))
        entry = dict(rust_all_paths=max_rss_bytes([str(binary), str(path)]))
        for name in LalCalls(case).paths():
            code = LAL_MEMORY.format(report=str(ROOT / "report"), case=json.dumps(case),
                                     path=name)
            entry[f"lal_{name}"] = max_rss_bytes([sys.executable, "-c", code])
        out[case["label"]] = entry
        print(f"memory {case['label']:<34} rust {entry['rust_all_paths'] / 1e6:7.0f} MB, lal "
              + ", ".join(f"{k[4:]} {v / 1e6:.0f} MB" for k, v in entry.items()
                          if k.startswith("lal_")), flush=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--rust-repeats", type=int, default=2,
                        help="repeats of each Rust path per round")
    parser.add_argument("--grids", nargs="*", default=list(GRIDS))
    parser.add_argument("--memory", action="store_true")
    parser.add_argument("--output", default=str(OUTPUT))
    args = parser.parse_args()
    if os.environ.get("OMP_NUM_THREADS") != "1":
        raise SystemExit("run with OMP_NUM_THREADS=1")
    subprocess.run(["cargo", "build", "--quiet", "--release", "--example", "long_benchmark"],
                   cwd=ROOT, check=True)
    case_list = cases(args.grids)
    lal_calls = [LalCalls(case) for case in case_list]
    for calls in lal_calls:  # warm up (LAL loads tables on first use)
        for path in calls.paths().values():
            path()
    lal_best = [dict() for _ in case_list]
    rust_best = {case["label"]: dict() for case in case_list}
    load = [os.getloadavg()[0]]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for round_index in range(args.rounds):
            for best, calls in zip(lal_best, lal_calls):
                for name, path in calls.paths().items():
                    best[name] = min(best.get(name, np.inf), seconds(path))
            rust = run_rust(case_list, args.rust_repeats, tmp / "cases.json")
            for label, entry in rust.items():
                for name, value in entry.items():
                    if name in ("per_call", "grid", "uniform") and value is not None:
                        rust_best[label][name] = min(rust_best[label].get(name, np.inf), value)
                    elif name not in ("per_call", "grid", "uniform"):
                        rust_best[label][name] = value
            load.append(os.getloadavg()[0])
            print(f"round {round_index + 1}/{args.rounds} done (load {load[-1]:.1f})", flush=True)
        memory = measure_memory(case_list, tmp) if args.memory else None

    rustc = subprocess.run(["rustc", "--version"], capture_output=True, text=True).stdout.strip()
    out = dict(
        lalsimulation_version=ls.__version__,
        rustc=rustc,
        machine=platform.platform(),
        processor=platform.processor() or platform.machine(),
        date=datetime.date.today().isoformat(),
        rounds=args.rounds,
        rust_repeats_per_round=args.rust_repeats,
        load_average=load,
        grids={name: dict(GRIDS[name], n_frequencies=len(frequencies(GRIDS[name])))
               for name in args.grids},
        binaries=dict(aligned=ALIGNED, precessing=PRECESSING),
        cases=[dict(label=case["label"], model=case["model"], version=case["version"],
                    grid=case["grid"], lal=lal, rust=rust_best[case["label"]])
               for case, lal in zip(case_list, lal_best)],
        memory_bytes=memory,
    )
    Path(args.output).write_text(json.dumps(out, indent=1))

    def ms(value):
        return f"{value * 1e3:9.1f}" if value is not None else f"{'-':>9}"

    print(f"\n{'case':<34}{'LAL seq':>9}{'Rust':>9}{'grid':>9}{'x':>6}{'x grid':>7}"
          f"{'LAL uni':>9}{'LAL MB':>9}{'Rust uni':>9}   [ms]")
    for entry in out["cases"]:
        lal, rust = entry["lal"], entry["rust"]
        print(f"{entry['label']:<34}{ms(lal['sequence'])}{ms(rust['per_call'])}"
              f"{ms(rust['grid'])}{lal['sequence'] / rust['per_call']:6.1f}"
              f"{lal['sequence'] / rust['grid']:7.1f}{ms(lal.get('uniform'))}"
              f"{ms(lal.get('uniform_multibanding'))}{ms(rust.get('uniform'))}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
