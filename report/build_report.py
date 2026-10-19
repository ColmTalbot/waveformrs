#!/usr/bin/env python
"""
Build the technical report from the data in ``report/data``.

Writes

* ``report/waveformrs-report.html``: standalone HTML document,
* ``report/artifact/index.html``: the same page as a fragment (no html/head/body),
* ``report/waveformrs-report.pdf``: printed with headless Chrome (optional).

Usage
-----
    python report/build_report.py [--chrome PATH] [--no-pdf]

Run ``scripts/validate_against_lalsimulation.py`` and ``report/lal_benchmark.py``
first (see the "Reproducing the results" section of the report).
"""
import argparse
import datetime
import html
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = HERE / "data"
sys.path.insert(0, str(HERE))

import make_figures  # noqa: E402

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

# Rows of the summary and timing tables: (validation name, benchmark name, label, status)
MODELS = [
    ("TaylorF2", "TaylorF2", "TaylorF2", "fixed"),
    ("IMRPhenomD", "IMRPhenomD", "IMRPhenomD", "fixed"),
    ("IMRPhenomXAS", "IMRPhenomXAS", "IMRPhenomXAS", "new"),
    ("IMRPhenomXHM", "IMRPhenomXHM", "IMRPhenomXHM", "new"),
    ("IMRPhenomXP_NNLO", "IMRPhenomXP (102)", "IMRPhenomXP, NNLO (102)", "new"),
    ("IMRPhenomXP_MSA", "IMRPhenomXP (223)", "IMRPhenomXP, MSA (223)", "new"),
    ("IMRPhenomXP", "IMRPhenomXP (320)", "IMRPhenomXP, SpinTaylor (320)", "new"),
    ("IMRPhenomXPHM_NNLO", "IMRPhenomXPHM (102)", "IMRPhenomXPHM, NNLO (102)", "new"),
    ("IMRPhenomXPHM_MSA", "IMRPhenomXPHM (223)", "IMRPhenomXPHM, MSA (223)", "new"),
    ("IMRPhenomXPHM", "IMRPhenomXPHM (320)", "IMRPhenomXPHM, SpinTaylor (320)", "new"),
    ("IMRPhenomXPNR", "IMRPhenomXPNR", "IMRPhenomXPNR (6.2.0, no parity fix)", "new"),
    ("IMRPhenomXPNR_parityfix", "IMRPhenomXPNR", "IMRPhenomXPNR with the parity fix (patched LAL)", "new"),
    ("NRSur7dq4", "NRSur7dq4", "NRSur7dq4 (time domain)", "new"),
]
# Rows of the serial-evaluation table: (serial benchmark name, label)
SERIAL = [
    ("TaylorF2", "TaylorF2"),
    ("IMRPhenomD", "IMRPhenomD"),
    ("IMRPhenomXAS", "IMRPhenomXAS"),
    ("IMRPhenomXHM", "IMRPhenomXHM"),
    ("IMRPhenomXP (102)", "IMRPhenomXP, NNLO (102)"),
    ("IMRPhenomXP (223)", "IMRPhenomXP, MSA (223)"),
    ("IMRPhenomXP (320)", "IMRPhenomXP, SpinTaylor (320)"),
    ("IMRPhenomXPHM (102)", "IMRPhenomXPHM, NNLO (102)"),
    ("IMRPhenomXPHM (223)", "IMRPhenomXPHM, MSA (223)"),
    ("IMRPhenomXPHM (320)", "IMRPhenomXPHM, SpinTaylor (320)"),
    ("IMRPhenomXPNR", "IMRPhenomXPNR"),
    ("NRSur7dq4", "NRSur7dq4 (time domain)"),
]
STATUS = {
    "fixed": '<span class="chip warn">fixed</span>',
    "new": '<span class="chip info">new</span>',
}

TEST_LABELS = {
    "taylorf2_matches_lalsimulation": "TaylorF2 polarizations",
    "imrphenomd_matches_lalsimulation": "IMRPhenomD polarizations",
    "imrphenomxas_matches_lalsimulation": "IMRPhenomXAS polarizations",
    "imrphenomxhm_matches_lalsimulation": "IMRPhenomXHM polarizations",
    "imrphenomxhm_modes_match_lalsimulation": "IMRPhenomXHM individual modes (both signs of m)",
    "imrphenomxp_nnlo_matches_lalsimulation": "IMRPhenomXP, NNLO angles",
    "imrphenomxp_msa_matches_lalsimulation": "IMRPhenomXP, MSA angles (incl. MSA → NNLO fallback)",
    "imrphenomxp_spintaylor_matches_lalsimulation": "IMRPhenomXP, SpinTaylor angles",
    "imrphenomxp_spintaylor_fallback_matches_lalsimulation": "IMRPhenomXP, SpinTaylor → MSA fallback",
    "imrphenomxp_uniform_matches_lalsimulation": "IMRPhenomXP, uniform grids",
    "imrphenomxp_spintaylor_angles_match_lalsimulation": "SpinTaylor Euler angles (SimIMRPhenomXPSpinTaylorAngles)",
    "imrphenomxphm_spintaylor_matches_lalsimulation": "IMRPhenomXPHM, SpinTaylor angles",
    "imrphenomxphm_nnlo_matches_lalsimulation": "IMRPhenomXPHM, NNLO angles",
    "imrphenomxphm_msa_matches_lalsimulation": "IMRPhenomXPHM, MSA angles",
    "imrphenomxphm_spintaylor_fallback_matches_lalsimulation": "IMRPhenomXPHM, SpinTaylor → MSA fallback",
    "imrphenomxphm_uniform_matches_lalsimulation": "IMRPhenomXPHM, uniform grids",
    "imrphenomxphm_mode_arrays_match_lalsimulation": "IMRPhenomXPHM, restricted mode arrays",
    "spintaylor_t4_matches_lalsimulation": "SpinTaylorT4 orbit evolutions",
    "nrsur7dq4_polarizations_match_lalsimulation": "NRSur7dq4 polarizations (ChooseTDWaveform)",
    "nrsur7dq4_modes_match_lalsimulation": "NRSur7dq4 modes (ChooseTDModes)",
    "nrsur7dq4_dynamics_match_lalsimulation": "NRSur7dq4 dynamics",
    "hdf5_reader_matches_h5py": "HDF5 reader against h5py (datasets)",
    "imrphenomxhm_122019_matches_lalsimulation": "IMRPhenomXHM release 122019 polarizations",
    "imrphenomxhm_122019_modes_match_lalsimulation": "IMRPhenomXHM release 122019 modes",
    "imrphenomxpnr_matches_lalsimulation": "IMRPhenomXPNR polarizations (incl. every regime switch)",
    "imrphenomxpnr_parity_fix_matches_patched_lalsimulation": "IMRPhenomXPNR with the parity fix (6.2.0 + lalsuite !2593)",
    "imrphenomxpnr_uniform_matches_lalsimulation": "IMRPhenomXPNR, uniform grids (with and without the fix)",
    "imrphenomxpnr_choose_fd_modes_matches_lalsimulation": "IMRPhenomXPNR modes (ChooseFDModes)",
    "imrphenomxpnr_sequence_modes_match_lalsimulation": "IMRPhenomXPNR modes on frequency sequences",
    "choose_fd_modes_matches_lalsimulation": "Modes of IMRPhenomXHM/XPHM (ChooseFDModes)",
    "modes_on_frequency_sequences_match_lalsimulation": "Modes on frequency sequences",
    "choose_fd_model_matches_lalsimulation": "Common interface against the reference suites",
    "conditioned_fd_matches_lalsimulation": "SimInspiralFD conditioning of the FD models",
    "nrsur7dq4_fd_matches_lalsimulation": "SimInspiralFD, NRSur7dq4",
    "nrsur7dq4_td_matches_lalsimulation": "SimInspiralTD, NRSur7dq4",
    "taylorf2_td_matches_lalsimulation": "SimInspiralTD, TaylorF2",
    "imrphenomd_td_matches_lalsimulation": "SimInspiralTD, IMRPhenomD",
    "imrphenomxas_td_matches_lalsimulation": "SimInspiralTD, IMRPhenomXAS",
    "imrphenomxhm_td_matches_lalsimulation": "SimInspiralTD, IMRPhenomXHM",
    "imrphenomxp_td_matches_lalsimulation": "SimInspiralTD, IMRPhenomXP",
    "imrphenomxphm_td_matches_lalsimulation": "SimInspiralTD, IMRPhenomXPHM",
    "imrphenomxpnr_td_matches_lalsimulation": "SimInspiralTD, IMRPhenomXPNR",
}


def sci(value, digits=1):
    """Scientific notation with a real superscript, e.g. 2.2×10⁻¹¹."""
    if value is None:
        return "–"
    if value <= 0:
        return "0"
    mantissa, exponent = f"{value:.{digits}e}".split("e")
    return f"{mantissa}×10<sup>{int(exponent)}</sup>"


def duration(seconds):
    if seconds >= 1e-3:
        return f"{seconds * 1e3:.2f} ms"
    return f"{seconds * 1e6:.0f} µs"


def factor(value):
    return f"{value:.1f}×"


def test_output():
    """Output of all integration tests, one test at a time."""
    env = dict(os.environ)
    env.setdefault("NRSUR7DQ4_DATA", str(ROOT / "lal_data" / "NRSur7dq4_v1.0.h5"))
    return subprocess.run(
        ["cargo", "test", "--release", "--tests", "--", "--nocapture", "--test-threads=1"],
        cwd=ROOT, capture_output=True, text=True, check=True, env=env,
    ).stdout


def run_long_tests(output):
    """Largest errors of the long-waveform tests, by model and path."""
    pattern = re.compile(
        r"(?:test \w+ \.\.\. )?(\w+(?: \(\d+\))?) (\S+) \S+ \[(\d+) frequencies\]\s+"
        r"(per-freq|grid)\s+plus ([0-9.eE+-]+) .*cross ([0-9.eE+-]+)"
    )
    models = dict()
    for line in output.splitlines():
        match = pattern.search(line)
        if not match:
            continue
        model, case, n, path, plus, cross = match.groups()
        entry = models.setdefault(model, dict(cases=set(), n=0, per_freq=0.0, grid=0.0))
        entry["cases"].add(case)
        entry["n"] = max(entry["n"], int(n))
        key = "per_freq" if path == "per-freq" else "grid"
        entry[key] = max(entry[key], float(plus), float(cross))
    return models


def run_tests(output):
    """Parse the per-case output of the LALSimulation comparison tests."""
    suites = dict()
    for output in [output]:
        current = None
        for line in output.splitlines():
            match = re.match(r"test (\w+) \.\.\. ?(.*)", line)
            if match:
                current = match.group(1)
                suites.setdefault(current, dict(n=0, err=0.0, mismatch=None))
                line = match.group(2)
            if current is None:
                continue
            err = re.search(r"max (?:rel(?:ative)?|abs) err(?:or)? ([0-9.eE+-]+)", line)
            angle = re.search(r"max \|angle difference\| ([0-9.eE+-]+)", line)
            datasets = re.search(r"([0-9]+) datasets match h5py", line)
            if datasets:
                suites[current]["n"] += int(datasets.group(1))
                continue
            mismatch = re.search(r"mismatch ([0-9.eE+-]+)", line)
            if err or angle:
                value = float((err or angle).group(1))
                suite = suites[current]
                suite["n"] += 1
                suite["err"] = max(suite["err"], value)
                if mismatch:
                    m = float(mismatch.group(1))
                    suite["mismatch"] = m if suite["mismatch"] is None else max(suite["mismatch"], m)
    return {name: value for name, value in suites.items() if value["n"] > 0}


def python_test_count():
    """Number of tests of the Python bindings (collected by pytest)."""
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q"],
        cwd=ROOT / "python", capture_output=True, text=True, check=True,
    )
    match = re.search(r"(\d+) tests? collected", result.stdout)
    return int(match.group(1))


def git_log():
    out = subprocess.run(
        ["git", "log", "--reverse", "--format=%h%x09%s", "main..HEAD"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout
    return [line.split("\t", 1) for line in out.strip().splitlines()]


def build(args):
    make_figures.style()
    results = make_figures.load()
    make_figures.accuracy(results)
    make_figures.speedup(results)
    make_figures.summary(results)

    summary = {row["approximant"]: row for row in json.loads((DATA / "summary.json").read_text())}
    timings = json.loads((DATA / "timings.json").read_text())
    baseline_time = json.loads((DATA / "timings_baseline_main.json").read_text())
    baseline_acc = json.loads((DATA / "baseline_main_accuracy.json").read_text())
    rust, lal_seq, lal_uni = timings["rust"], timings["lal_sequence"], timings["lal_uniform"]

    values = dict()
    values["css"] = (HERE / "report.css").read_text()
    values["lalsim_version"] = timings["lalsimulation_version"]
    values["lal_version"] = "7.7.0"
    values["rustc"] = subprocess.run(
        ["rustc", "--version"], capture_output=True, text=True
    ).stdout.split(" (")[0]
    values["date"] = datetime.date.today().strftime("%-d %B %Y")
    commits = [
        (sha, subject)
        for sha, subject in git_log()
        if not subject.startswith(("DOC: add the technical report", "DOC: update the technical report"))
    ]
    values["n_commits"] = str(len(commits))

    n_total = sum(summary[m[0]]["n_draws"] for m in MODELS if m[0] in summary)
    values["n_draws_total"] = f"{n_total:,}"
    values["max_rel_all"] = sci(max(summary[m[0]]["max_relative_error"] for m in MODELS))
    values["max_mm_all"] = sci(max(summary[m[0]]["max_mismatch"] for m in MODELS))

    # reference files generated with lalsimulation (the HDF5 checksum dump
    # comes from h5py)
    references = [
        json.loads(path.read_text())
        for path in sorted((ROOT / "tests" / "data").glob("*.json"))
    ]
    references = [ref for ref in references if "cases" in ref]
    values["n_ref_files"] = str(len(references))
    values["n_ref_cases"] = str(sum(len(ref["cases"]) for ref in references))

    values["xphm320_rust"] = duration(rust["IMRPhenomXPHM (320)"])
    values["xphm320_lal"] = duration(lal_seq["IMRPhenomXPHM (320)"])
    values["xphm320_speedup"] = factor(lal_seq["IMRPhenomXPHM (320)"] / rust["IMRPhenomXPHM (320)"]) + " faster"
    values["xphm320_lal_mb"] = duration(lal_uni["IMRPhenomXPHM (320) (uniform, default multibanding)"])
    values["xphm320_rust_uniform"] = duration(rust["IMRPhenomXPHM (320) (uniform)"])
    values["tf2_vs_main"] = factor(baseline_time["TaylorF2 new + waveform_modes_us"] * 1e-6 / rust["TaylorF2"])
    values["d_vs_main"] = factor(baseline_time["IMRPhenomD new + waveform_modes_us"] * 1e-6 / rust["IMRPhenomD"])

    rows = []
    for name, bench, label, status in MODELS:
        row = summary[name]
        speed = lal_seq[bench] / rust[bench]
        rows.append(
            f"<tr><td>{label}</td><td>{STATUS[status]}</td><td class=\"num\">{row['n_compared']}/{row['n_draws']}</td>"
            f"<td class=\"num\">{sci(row['max_relative_error'])}</td><td class=\"num\">{sci(row['max_mismatch'])}</td>"
            f"<td class=\"num\">{duration(rust[bench])}</td><td class=\"num\">{duration(lal_seq[bench])}</td>"
            f"<td class=\"num\">{factor(speed)}</td></tr>"
        )
    values["summary_rows"] = "\n".join(rows)

    rows = []
    for name, key in [("TaylorF2", "TaylorF2"), ("IMRPhenomD", "IMRPhenomD")]:
        records = baseline_acc[name]
        mm = np.array([r["mismatch"] for r in records])
        amp = np.array([r["amplitude_error"] for r in records])
        rows.append(
            f"<tr><td>{name}</td><td class=\"num\">{sci(np.median(mm))}</td><td class=\"num\">{sci(np.percentile(mm, 90))}</td>"
            f"<td class=\"num\">{sci(mm.max())}</td><td class=\"num\">{int(np.sum(mm > 1e-3))} / {len(mm)}</td>"
            f"<td class=\"num\">{sci(amp.max())}</td><td class=\"num\">{duration(baseline_time[f'{key} new + waveform_modes_us'] * 1e-6)}</td></tr>"
        )
    values["baseline_rows"] = "\n".join(rows)

    output = test_output()
    suites = run_tests(output)
    long_models = run_long_tests(output)
    rows = []
    for test, label in TEST_LABELS.items():
        if test not in suites:
            continue
        suite = suites[test]
        bit = ' <span class="chip good">bit-identical</span>' if suite["err"] == 0 else ""
        mismatch = sci(suite["mismatch"]) if suite["mismatch"] and suite["mismatch"] > 0 else ("≤ 0" if suite["mismatch"] is not None else "–")
        rows.append(
            f"<tr><td>{label}<br><code>{test}</code></td><td class=\"num\">{suite['n']}</td>"
            f"<td class=\"num\">{sci(suite['err'])}{bit}</td><td class=\"num\">{mismatch}</td></tr>"
        )
    values["suite_rows"] = "\n".join(rows)

    rows = []
    for name, bench, label, status in MODELS:
        row = summary[name]
        rows.append(
            f"<tr><td>{label}</td><td class=\"num\">{row['n_compared']}/{row['n_draws']}</td>"
            f"<td class=\"num\">{sci(row['max_relative_error'])}</td><td class=\"num\">{sci(row['median_relative_error'])}</td>"
            f"<td class=\"num\">{sci(row['max_mismatch'])}</td><td class=\"num\">{sci(row['median_mismatch'])}</td></tr>"
        )
    values["random_rows"] = "\n".join(rows)

    rows = []
    for name, bench, label, status in MODELS:
        key = f"{bench} new + waveform_modes_us"
        original = duration(baseline_time[key] * 1e-6) if key in baseline_time else "–"
        rows.append(
            f"<tr><td>{label}</td><td class=\"num\">{original}</td><td class=\"num\">{duration(rust[bench])}</td>"
            f"<td class=\"num\">{duration(lal_seq[bench])}</td><td class=\"num\">{factor(lal_seq[bench] / rust[bench])}</td></tr>"
        )
    values["timing_rows"] = "\n".join(rows)

    rows = []
    for name, label in [("IMRPhenomXHM", "IMRPhenomXHM"), ("IMRPhenomXPHM (320)", "IMRPhenomXPHM, SpinTaylor (320)")]:
        t_rust = rust[f"{name} (uniform)"]
        t_off = lal_uni[f"{name} (uniform)"]
        t_on = lal_uni[f"{name} (uniform, default multibanding)"]
        rows.append(
            f"<tr><td>{label}</td><td class=\"num\">{duration(t_rust)}</td><td class=\"num\">{duration(t_off)}</td>"
            f"<td class=\"num\">{duration(t_on)}</td><td class=\"num\">{factor(t_on / t_rust)}</td></tr>"
        )
    values["uniform_rows"] = "\n".join(rows)

    for figure in ["accuracy", "speedup"]:
        svg = (HERE / "figures" / f"{figure}.svg").read_text()
        svg = re.sub(r"<metadata>.*?</metadata>\s*", "", svg, flags=re.DOTALL)
        label = {"accuracy": "Relative error against total mass for each model", "speedup": "Speed-up over LALSimulation for each model"}[figure]
        svg = svg.replace("<svg ", f'<svg role="img" aria-label="{label}" ', 1)
        values[f"{figure}_svg"] = svg

    kinds = {"FEAT": "info", "BUG": "bad", "PERF": "good", "TEST": "", "DOC": "", "MAINT": "", "EXAMPLE": "", "STY": ""}
    rows = []
    for sha, subject in commits:
        kind, _, rest = subject.partition(": ")
        chip = f'<span class="chip {kinds.get(kind, "")}">{html.escape(kind)}</span>' if kind in kinds else ""
        text = html.escape(rest if chip else subject)
        text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
        rows.append(f"<tr><td class=\"mono\">{sha}</td><td>{chip}</td><td>{text}</td></tr>")
    values["commit_rows"] = "\n".join(rows)

    # NRSur7dq4
    values["nrsur_speedup"] = factor(lal_seq["NRSur7dq4"] / rust["NRSur7dq4"])
    values["hdf5_lines"] = f"{len((ROOT / 'src' / 'hdf5.rs').read_text().splitlines()):,}"
    nrsur = json.loads((DATA / "nrsur_timings.json").read_text())
    extrinsic = [c["extrinsic"] for c in nrsur["rust"]["configurations"]]
    values["nrsur_extrinsic"] = f"{min(extrinsic) * 1e3:.2f}–{max(extrinsic) * 1e3:.2f} ms"
    values["nrsur_extrinsic_range"] = values["nrsur_extrinsic"]
    values["nrsur_load"] = duration(nrsur["rust"]["load_cached_file"])
    rows = []
    for r, l in zip(nrsur["rust"]["configurations"], nrsur["lal"]["configurations"]):
        rows.append(
            f"<tr><td class=\"num\">{r['total_mass']:g}</td><td class=\"num\">{r['sampling_rate']:g}</td>"
            f"<td class=\"num\">{r['f_min']:g}</td><td class=\"num\">{r['samples']:,}</td>"
            f"<td class=\"num\">{duration(l['polarizations'])}</td><td class=\"num\">{duration(r['fresh'])}</td>"
            f"<td class=\"num\">{duration(r['reused'])}</td><td class=\"num\">{duration(r['extrinsic'])}</td>"
            f"<td class=\"num\">{duration(l['modes'])}</td><td class=\"num\">{duration(r['modes'])}</td></tr>"
        )
    values["nrsur_rows"] = "\n".join(rows)

    # many waveforms in sequence
    serial = json.loads((DATA / "serial_timings.json").read_text())
    values["serial_n"] = str(serial["n_draws"])
    values["serial_nf"] = f"{serial['n_frequencies']:,}"
    values["serial_overhead"] = duration(serial["lal_python_call_overhead"])
    rows, grid_speedups = [], []
    for name, label in SERIAL:
        if name not in serial["rust"]:
            continue
        r, l = serial["rust"][name], serial["lal"][name]["seconds"]
        grid = r.get("grid_seconds")
        if grid:
            grid_speedups.append(l / grid)
        rows.append(
            f"<tr><td>{label}</td><td class=\"num\">{duration(r['seconds'])}</td>"
            f"<td class=\"num\">{duration(grid) if grid else '–'}</td><td class=\"num\">{duration(l)}</td>"
            f"<td class=\"num\">{factor(l / r['seconds'])}</td>"
            f"<td class=\"num\">{factor(l / grid) if grid else '–'}</td></tr>"
        )
    values["serial_rows"] = "\n".join(rows)
    values["serial_best_range"] = f"{factor(min(grid_speedups))}–{factor(max(grid_speedups))}"
    xp = [serial["rust"].get(n, {}).get("grid_seconds") for n, _ in SERIAL if n.startswith("IMRPhenomXP")]
    if all(xp):
        speeds = [serial["rust"][n]["seconds"] / serial["rust"][n]["grid_seconds"]
                  for n, _ in SERIAL if n.startswith("IMRPhenomXP")]
        values["xp_grid_text"] = (
            f"Their grid paths are {factor(min(speeds))}–{factor(max(speeds))} faster than the "
            "per-call paths. They work in passes over the grid. The Euler angles are computed for "
            "each m′, the NNLO angles from the stored powers without transcendental functions, and "
            "the MSA angles in blocks of 32 frequencies. An MSA angle already computed for another "
            "m′ at a bit-identical argument is reused, which skips 29% of the evaluations for "
            "IMRPhenomXPHM. The ε phase is folded into each mode's phase before one vectorised "
            "conversion, and the twist-up keeps its sums in registers. What remains is the MSA "
            "angles, whose spin-cubic roots and corrections must use the C library to stay stable "
            "for nearly aligned spins; the SpinTaylor PN integration, about 0.5 ms; and the "
            "IMRPhenomXHM mode evaluations."
        )
    else:
        values["xp_grid_text"] = (
            "They have no grid path yet: nearly all their time is spent per frequency, in the "
            "Euler angles (most of all the MSA angles) and in the co-precessing modes."
        )

    # long waveforms
    long = json.loads((DATA / "long_timings.json").read_text())
    sizes = [m["n"] for m in long_models.values()]
    values["long_min_n"] = f"{min(sizes) // 1000}k"
    values["long_max_n"] = f"{max(sizes) // 1000}k"
    rows = []
    for model, entry in long_models.items():
        rows.append(
            f"<tr><td>{model}</td><td class=\"num\">{len(entry['cases'])}</td><td class=\"num\">{entry['n']:,}</td>"
            f"<td class=\"num\">{sci(entry['per_freq'])}</td><td class=\"num\">{sci(entry['grid'])}</td></tr>"
        )
    values["long_accuracy_rows"] = "\n".join(rows)
    grid_label = max(long["grids"], key=lambda g: long["grids"][g]["n_frequencies"])
    values["long_grid_label"] = grid_label
    values["long_grid_n"] = f"{long['grids'][grid_label]['n_frequencies']:,}"
    values["long_rounds"] = str(long["rounds"])
    rows = []
    for case in long["cases"]:
        if case["grid"] != grid_label:
            continue
        lal_t, rust_t = case["lal"], case["rust"]
        mb = lal_t.get("uniform_multibanding")
        rows.append(
            f"<tr><td>{html.escape(case['model'])}{' (' + str(case['version']) + ')' if case['version'] else ''}</td>"
            f"<td class=\"num\">{duration(lal_t['sequence'])}</td><td class=\"num\">{duration(rust_t['per_call'])}</td>"
            f"<td class=\"num\">{duration(rust_t['grid'])}</td><td class=\"num\">{duration(mb) if mb else '–'}</td>"
            f"<td class=\"num\">{factor(lal_t['sequence'] / rust_t['per_call'])}</td></tr>"
        )
    values["long_timing_rows"] = "\n".join(rows)
    memory = (long.get("memory_bytes") or {}).get(f"IMRPhenomXPHM (320), {grid_label}")
    if not memory:
        raise SystemExit("long_timings.json has no memory measurements: run report/long_benchmark.py --memory")
    values["long_mem_rust"] = f"{memory['rust_all_paths'] / 2**20:.0f} MB"
    values["long_mem_lal"] = f"{memory['lal_sequence'] / 2**30:.2f} GB"

    # IMRPhenomXPNR
    xpnr_cases = 0
    for name in ["imrphenomxpnr", "imrphenomxpnr_uniform", "imrphenomxpnr_modes", "imrphenomxpnr_sequence_modes"]:
        xpnr_cases += len(json.loads((ROOT / "tests" / "data" / f"{name}.json").read_text())["cases"])
    values["xpnr_cases"] = str(xpnr_cases)
    values["xpnr_rust"] = duration(rust["IMRPhenomXPNR"])
    values["xpnr_lal"] = duration(lal_seq["IMRPhenomXPNR"])
    values["xpnr_lal_mb"] = duration(lal_uni["IMRPhenomXPNR (uniform, default multibanding)"])
    fixed = json.loads((DATA / "validation_IMRPhenomXPNR_parityfix.json").read_text())
    values["xpnr_fix_n"] = str(fixed["n_compared"])
    values["xpnr_fix_rel"] = sci(fixed["max_relative_error"])
    values["xpnr_fix_mm"] = sci(fixed["max_mismatch"])
    build = fixed["lalsimulation_build"]
    values["xpnr_fix_build"] = html.escape(build.get("vcsStatus", "") or "")
    xpnr_serial = serial["rust"]["IMRPhenomXPNR"]
    values["xpnr_seq_time"] = duration(xpnr_serial["seconds"])
    values["xpnr_grid_time"] = duration(xpnr_serial["grid_seconds"])
    values["xpnr_lal_time"] = duration(serial["lal"]["IMRPhenomXPNR"]["seconds"])
    values["xpnr_grid_speedup"] = factor(xpnr_serial["seconds"] / xpnr_serial["grid_seconds"])

    # Python bindings
    python = json.loads((DATA / "python_timings.json").read_text())
    values["py_machine"] = html.escape(f"{python['machine']}, {python['cpus']} CPUs")
    rows = []
    for label, entry in python["overhead"].items():
        rows.append(
            f"<tr><td>{html.escape(label)}</td><td class=\"num\">{duration(entry['waveformrs'])}</td>"
            f"<td class=\"num\">{duration(entry['lalsimulation'])}</td></tr>"
        )
    values["py_overhead_rows"] = "\n".join(rows)
    rows = []
    for label, entry in python["throughput"].items():
        grid = entry.get("on_grid")
        best = grid or entry["waveformrs"]
        rows.append(
            f"<tr><td>{html.escape(label)}</td><td class=\"num\">{duration(entry['waveformrs'])}</td>"
            f"<td class=\"num\">{duration(grid) if grid else '–'}</td>"
            f"<td class=\"num\">{duration(entry['lalsimulation'])}</td>"
            f"<td class=\"num\">{factor(entry['lalsimulation'] / best)}</td></tr>"
        )
    values["py_throughput_rows"] = "\n".join(rows)
    thread_keys = list(next(iter(python["batch"].values()))["waveformrs"].keys())
    values["py_batch_header"] = "".join(
        f"<th class=\"num\">{'all threads' if k == '0' else k + (' thread' if k == '1' else ' threads')}</th>"
        for k in thread_keys
    )
    rows = []
    for label, entry in python["batch"].items():
        cells = "".join(f"<td class=\"num\">{entry['waveformrs'][k]:,.0f}</td>" for k in thread_keys)
        rows.append(
            f"<tr><td>{html.escape(label)}</td><td class=\"num\">{entry['lalsimulation']:,.0f}</td>{cells}</tr>"
        )
    values["py_batch_rows"] = "\n".join(rows)
    values["py_n_tests"] = str(python_test_count())

    # batch scaling
    batch = json.loads((DATA / "batch_scaling.json").read_text())
    values["batch_nf"] = f"{batch['frequencies']:,}"
    values["batch_header"] = "".join(
        f"<th class=\"num\">{t} thread{'s' if t > 1 else ''}</th>" for t in batch["threads"]
    )
    rows = []
    for row in batch["rows"]:
        rates = row["waveforms_per_second"]
        cells = "".join(
            f"<td class=\"num\">{rate:,.0f}{'' if i == 0 else f' ({rate / rates[0]:.1f}×)'}</td>"
            for i, rate in enumerate(rates)
        )
        rows.append(f"<tr><td>{html.escape(row['model'])}</td>{cells}</tr>")
    values["batch_rows"] = "\n".join(rows)

    # conversions
    conversions = json.loads((DATA / "conversion_timings.json").read_text())
    rows = []
    for row in conversions["rows"]:
        lal_t = row["lal"]
        what, label = row["what"], row["label"]
        if what.startswith("("):  # FFT rows: "(plan + FFT)", "forward real FFT, N = ..."
            what, label = label, what.strip("()")
        rows.append(
            f"<tr><td>{html.escape(what)}</td><td>{html.escape(label)}</td>"
            f"<td class=\"num\">{duration(row['rust'])}</td><td class=\"num\">{duration(lal_t) if lal_t else '–'}</td>"
            f"<td class=\"num\">{factor(lal_t / row['rust']) if lal_t else '–'}</td></tr>"
        )
    values["conversion_rows"] = "\n".join(rows)

    # other platforms
    comparison = json.loads((DATA / "linux_reference_comparison.json").read_text())
    arm = comparison["lalsimulation linux-arm64 vs osx-arm64"]
    amd = comparison["lalsimulation linux-amd64 vs osx-arm64"]
    rows = []
    for suite in sorted(arm):
        rows.append(
            f"<tr><td><code>{suite}</code></td><td class=\"num\">{sci(arm[suite]['max_relative_error'])}</td>"
            f"<td class=\"num\">{sci(amd.get(suite, {}).get('max_relative_error'))}</td></tr>"
        )
    values["linux_lal_rows"] = "\n".join(rows)
    arm_models = json.loads((DATA / "linux_validation_arm64.json").read_text())["linux-arm64"]["models"]
    amd_models = json.loads((DATA / "linux_validation_amd64.json").read_text())["linux-amd64"]["models"]
    rows = []
    for model in sorted(set(arm_models) | set(amd_models)):
        cells = []
        for table in (arm_models, amd_models):
            entry = table.get(model)
            cells.append(
                f"{sci(entry['median_relative_error'])} / {sci(entry['max_relative_error'])}" if entry else "–"
            )
        rows.append(f"<tr><td>{model}</td><td class=\"num\">{cells[0]}</td><td class=\"num\">{cells[1]}</td></tr>")
    values["linux_validation_rows"] = "\n".join(rows)
    linux = json.loads((DATA / "linux_timings_arm64.json").read_text())["linux-arm64"]
    rows = []
    for model, entry in linux["single"].items():
        # the workspace row repeats one binary, so it times the cached modes
        if "lal_us" not in entry or "workspace" in model:
            continue
        grid_entry = linux.get("grid", {}).get(model, {})
        grid_t = grid_entry.get("grid_us")
        rows.append(
            f"<tr><td>{model}</td><td class=\"num\">{duration(entry['rust_us'] * 1e-6)}</td>"
            f"<td class=\"num\">{duration(entry['lal_us'] * 1e-6)}</td><td class=\"num\">{factor(entry['speedup'])}</td>"
            f"<td class=\"num\">{duration(grid_t * 1e-6) if grid_t else '–'}</td></tr>"
        )
    values["linux_timing_rows"] = "\n".join(rows)

    page = (HERE / "template.html").read_text()
    sections = re.findall(r'<section id="([a-z]+)">', page)
    for number, section in enumerate(sections, start=1):
        page = page.replace("{{sec:" + section + "}}", str(number))
    for key, value in values.items():
        page = page.replace("{{" + key + "}}", value)
    missing = re.findall(r"\{\{(\w+)\}\}", page)
    if missing:
        raise SystemExit(f"unfilled placeholders: {missing}")

    (HERE / "artifact").mkdir(exist_ok=True)
    (HERE / "artifact" / "index.html").write_text(page)
    document = (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        "</head>\n<body>\n" + page + "\n</body>\n</html>\n"
    )
    out = HERE / "waveformrs-report.html"
    out.write_text(document)
    print(f"wrote {out} and {HERE / 'artifact' / 'index.html'}")

    if not args.no_pdf:
        pdf = HERE / "waveformrs-report.pdf"
        subprocess.run(
            [
                args.chrome, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                "--virtual-time-budget=20000", f"--user-data-dir={tempfile.mkdtemp()}",
                f"--print-to-pdf={pdf}", out.as_uri(),
            ],
            check=True, capture_output=True,
        )
        print(f"wrote {pdf}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrome", default=CHROME)
    parser.add_argument("--no-pdf", action="store_true")
    build(parser.parse_args())


if __name__ == "__main__":
    main()
