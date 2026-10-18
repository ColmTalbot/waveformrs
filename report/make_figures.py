#!/usr/bin/env python
"""
Make the figures for the technical report from the validation results in
``report/data`` (written by ``scripts/validate_against_lalsimulation.py``).

The figures are written as SVG with the text and axis colours replaced by
CSS variables so that they follow the light/dark theme of the HTML report.

Usage
-----
    python report/make_figures.py
"""
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGURES = HERE / "figures"

# Sentinel colours replaced by CSS variables after rendering.
FOREGROUND = "#010203"
MUTED = "#040506"
GRID = "#070809"
REPLACEMENTS = {
    FOREGROUND: "var(--fg)",
    MUTED: "var(--muted)",
    GRID: "var(--rule)",
}

ORDER = [
    "TaylorF2",
    "IMRPhenomD",
    "IMRPhenomXAS",
    "IMRPhenomXHM",
    "IMRPhenomXP_NNLO",
    "IMRPhenomXP_MSA",
    "IMRPhenomXP",
    "IMRPhenomXPHM_NNLO",
    "IMRPhenomXPHM_MSA",
    "IMRPhenomXPHM",
    "IMRPhenomXPNR",
    "NRSur7dq4",
]
LABELS = {
    "IMRPhenomXP_NNLO": "IMRPhenomXP (102)",
    "IMRPhenomXP_MSA": "IMRPhenomXP (223)",
    "IMRPhenomXP": "IMRPhenomXP (320)",
    "IMRPhenomXPHM_NNLO": "IMRPhenomXPHM (102)",
    "IMRPhenomXPHM_MSA": "IMRPhenomXPHM (223)",
    "IMRPhenomXPHM": "IMRPhenomXPHM (320)",
    "NRSur7dq4": "NRSur7dq4 (bit-identical)",
}
# Mid-tone colours that read on both light and dark backgrounds.
COLOURS = [
    "#3b7dd8",
    "#d0612f",
    "#2f9e6e",
    "#9b59c9",
    "#c79a1c",
    "#d1477a",
    "#2a9bb3",
    "#7a8b2e",
    "#b8662d",
    "#5a6fd6",
    "#8c6d3f",
    "#3f8c7a",
]


def style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 9,
            "text.color": FOREGROUND,
            "axes.labelcolor": FOREGROUND,
            "axes.edgecolor": MUTED,
            "axes.linewidth": 0.6,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelcolor": FOREGROUND,
            "ytick.labelcolor": FOREGROUND,
            "grid.color": GRID,
            "grid.linewidth": 0.5,
            "axes.grid": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "figure.facecolor": "none",
            "axes.facecolor": "none",
            "savefig.facecolor": "none",
            "legend.frameon": False,
            "svg.hashsalt": "waveformrs",
        }
    )


def load():
    results = dict()
    for name in ORDER:
        path = DATA / f"validation_{name}.json"
        if path.exists():
            with open(path) as ff:
                results[name] = json.load(ff)
    return results


def save(fig, name):
    FIGURES.mkdir(exist_ok=True)
    path = FIGURES / f"{name}.svg"
    fig.savefig(path, bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    text = path.read_text()
    for sentinel, variable in REPLACEMENTS.items():
        text = re.sub(sentinel, variable, text, flags=re.IGNORECASE)
    # inline-friendly: drop the XML prolog and fixed size so CSS can scale it
    text = re.sub(r"<\?xml[^>]*\?>\s*", "", text)
    text = re.sub(r"<!DOCTYPE[^>]*>\s*", "", text, flags=re.DOTALL)
    text = re.sub(r'(<svg[^>]*?) width="[^"]*" height="[^"]*"', r"\1", text, count=1)
    path.write_text(text)
    print(f"wrote {path}")


def total_mass(record):
    return record["mass_1"] + record["mass_2"]


def rust_time(record):
    """Time of a fresh Rust call (the NRSur7dq4 harness also times a reused
    workspace, which can skip the surrogate when the binary repeats)."""
    return record.get("time_rust_fresh", record["time_rust"])


def accuracy(results):
    names = list(results)
    ncols = 4
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(8.0, 1.75 * nrows), sharex=True, sharey=True
    )
    axes = np.atleast_2d(axes)
    for ax, name in zip(axes.flat, names):
        records = results[name]["records"]
        mass = np.array([total_mass(r) for r in records])
        error = np.array(
            [max(r["max_rel_err_plus"], r["max_rel_err_cross"]) for r in records]
        )
        error = np.maximum(error, 1e-17)
        colour = COLOURS[ORDER.index(name) % len(COLOURS)]
        ax.scatter(mass, error, s=5, color=colour, linewidths=0, alpha=0.85)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(LABELS.get(name, name), fontsize=8.5, loc="left")
        ax.axhline(1e-10, color=MUTED, linewidth=0.6, linestyle="--")
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    axes[0, 0].set_ylim(1e-17, 1e-9)
    for ax in axes[-1]:
        ax.set_xlabel(r"total mass [$M_\odot$]")
    for ax in axes[:, 0]:
        ax.set_ylabel("max relative error")
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    save(fig, "accuracy")


def speedup(results):
    names = list(results)
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(names) + 0.7))
    rng = np.random.default_rng(1)
    for ii, name in enumerate(names):
        records = results[name]["records"]
        ratio = np.array([r["time_lal"] / rust_time(r) for r in records])
        colour = COLOURS[ORDER.index(name) % len(COLOURS)]
        yy = ii + rng.uniform(-0.18, 0.18, len(ratio))
        ax.scatter(ratio, yy, s=5, color=colour, linewidths=0, alpha=0.7)
        median = np.median(ratio)
        ax.plot([median, median], [ii - 0.32, ii + 0.32], color=FOREGROUND, lw=1.2)
        ax.annotate(
            f"×{median:.1f}",
            (median, ii + 0.32),
            xytext=(3, -1),
            textcoords="offset points",
            fontsize=7.5,
            va="bottom",
        )
    ax.set_yticks(range(len(names)), [LABELS.get(n, n) for n in names])
    ax.set_ylim(len(names) - 0.5, -0.6)
    ax.set_xscale("log")
    ax.axvline(1, color=MUTED, lw=0.8)
    ax.set_xlabel("speed-up over lalsimulation (single thread, same frequency grid)")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    save(fig, "speedup")


def summary(results):
    rows = list()
    for name, result in results.items():
        records = result["records"]
        ratio = np.array([r["time_lal"] / rust_time(r) for r in records])
        rows.append(
            dict(
                approximant=name,
                label=LABELS.get(name, name),
                n_draws=result["n_draws"],
                n_compared=result.get("n_compared", result["n_draws"]),
                max_mismatch=result["max_mismatch"],
                median_mismatch=result["median_mismatch"],
                max_relative_error=result["max_relative_error"],
                median_relative_error=result["median_relative_error"],
                median_speedup=float(np.median(ratio)),
                min_speedup=float(np.min(ratio)),
                max_speedup=float(np.max(ratio)),
                median_time_rust=float(np.median([rust_time(r) for r in records])),
                median_time_lal=float(np.median([r["time_lal"] for r in records])),
                failures=[
                    dict(
                        total_mass=f["mass_1"] + f["mass_2"],
                        rust_error=f.get("rust_error"),
                        lal_error=f.get("lal_error"),
                    )
                    for f in result.get("failures", [])
                ],
            )
        )
    with open(DATA / "summary.json", "w") as ff:
        json.dump(rows, ff, indent=1)
    print(f"wrote {DATA / 'summary.json'}")


def main():
    style()
    results = load()
    accuracy(results)
    speedup(results)
    summary(results)


if __name__ == "__main__":
    main()
