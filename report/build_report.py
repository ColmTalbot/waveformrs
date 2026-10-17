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


def run_tests():
    """Parse the per-case output of the LALSimulation comparison tests."""
    suites = dict()
    for test in ["lalsimulation", "spintaylor"]:
        output = subprocess.run(
            ["cargo", "test", "--release", "--test", test, "--", "--nocapture", "--test-threads=1"],
            cwd=ROOT, capture_output=True, text=True, check=True,
        ).stdout
        current = None
        for line in output.splitlines():
            match = re.match(r"test (\w+) \.\.\. ?(.*)", line)
            if match:
                current = match.group(1)
                suites.setdefault(current, dict(n=0, err=0.0, mismatch=None))
                line = match.group(2)
            if current is None:
                continue
            err = re.search(r"max rel(?:ative)? err(?:or)? ([0-9.eE+-]+)", line)
            angle = re.search(r"max \|angle difference\| ([0-9.eE+-]+)", line)
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
        if not subject.startswith("DOC: add the technical report")
    ]
    values["n_commits"] = str(len(commits))

    n_total = sum(summary[m[0]]["n_draws"] for m in MODELS if m[0] in summary)
    values["n_draws_total"] = f"{n_total:,}"
    values["max_rel_all"] = sci(max(summary[m[0]]["max_relative_error"] for m in MODELS))
    values["max_mm_all"] = sci(max(summary[m[0]]["max_mismatch"] for m in MODELS))

    ref_files = sorted((ROOT / "tests" / "data").glob("*.json"))
    values["n_ref_files"] = str(len(ref_files))
    values["n_ref_cases"] = str(
        sum(len(json.loads(path.read_text())["cases"]) for path in ref_files)
    )

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

    suites = run_tests()
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

    page = (HERE / "template.html").read_text()
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
