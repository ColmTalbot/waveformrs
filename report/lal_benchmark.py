#!/usr/bin/env python
"""
Time lalsimulation on the configuration of ``examples/benchmark.rs`` and
write ``report/data/timings.json`` with both sets of timings.

Run single-threaded (lalsimulation's IMRPhenomX* loops use OpenMP):

    OMP_NUM_THREADS=1 python report/lal_benchmark.py [--output FILE]
"""
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import lal
import lalsimulation as ls
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
os.environ.setdefault("LAL_DATA_PATH", str(ROOT / "lal_data"))

FREQUENCIES = np.arange(80, 6001) * 0.25
ALIGNED = dict(masses=(40.0, 20.0), spins=(0, 0, 0.3, 0, 0, -0.2))
PRECESSING = dict(masses=(40.0, 20.0), spins=(0.4, -0.3, 0.2, 0.1, 0.5, -0.3))
DISTANCE, PHI_REF, INCLINATION, F_REF = 400.0, 0.1, 0.4, 20.0
NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0.0, PhenomXPHMThresholdMband=0.0)

# name -> (lalsimulation approximant, binary, flags)
CASES = {
    "TaylorF2": ("TaylorF2", ALIGNED, dict()),
    "IMRPhenomD": ("IMRPhenomD", ALIGNED, dict()),
    "IMRPhenomXAS": ("IMRPhenomXAS", ALIGNED, dict()),
    "IMRPhenomXHM": ("IMRPhenomXHM", ALIGNED, NO_MULTIBANDING),
    "IMRPhenomXP (102)": ("IMRPhenomXP", PRECESSING, dict(PhenomXPrecVersion=102)),
    "IMRPhenomXP (223)": ("IMRPhenomXP", PRECESSING, dict(PhenomXPrecVersion=223)),
    "IMRPhenomXP (320)": ("IMRPhenomXP", PRECESSING, dict(PhenomXPrecVersion=320)),
    "IMRPhenomXPHM (102)": (
        "IMRPhenomXPHM",
        PRECESSING,
        dict(PhenomXPrecVersion=102, **NO_MULTIBANDING),
    ),
    "IMRPhenomXPHM (223)": (
        "IMRPhenomXPHM",
        PRECESSING,
        dict(PhenomXPrecVersion=223, **NO_MULTIBANDING),
    ),
    "IMRPhenomXPHM (320)": (
        "IMRPhenomXPHM",
        PRECESSING,
        dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
    ),
    # (the frequency-sequence driver switches multibanding off)
    "IMRPhenomXPNR": ("IMRPhenomXPNR", PRECESSING, dict()),
}
# Uniform grids (SimInspiralChooseFDWaveform, which stops at Mf = 0.3), with
# multibanding off and with the default multibanding thresholds.
UNIFORM = {
    "IMRPhenomXHM (uniform)": ("IMRPhenomXHM", ALIGNED, NO_MULTIBANDING),
    "IMRPhenomXHM (uniform, default multibanding)": ("IMRPhenomXHM", ALIGNED, dict()),
    "IMRPhenomXPHM (320) (uniform)": (
        "IMRPhenomXPHM",
        PRECESSING,
        dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
    ),
    "IMRPhenomXPHM (320) (uniform, default multibanding)": (
        "IMRPhenomXPHM",
        PRECESSING,
        dict(PhenomXPrecVersion=320),
    ),
    "IMRPhenomXPNR (uniform)": ("IMRPhenomXPNR", PRECESSING, NO_MULTIBANDING),
    "IMRPhenomXPNR (uniform, default multibanding)": ("IMRPhenomXPNR", PRECESSING, dict()),
}


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in flags.items():
        getattr(ls, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def best_time(func):
    func()
    start = time.perf_counter()
    calls = 0
    while calls < 3 or time.perf_counter() - start < 0.05:
        func()
        calls += 1
    number = max(3, int(0.2 / ((time.perf_counter() - start) / calls)))
    best = np.inf
    for _ in range(7):
        start = time.perf_counter()
        for _ in range(number):
            func()
        best = min(best, (time.perf_counter() - start) / number)
    return best


def sequence(approximant, binary, flags):
    frequencies = lal.CreateREAL8Vector(len(FREQUENCIES))
    frequencies.data = FREQUENCIES
    m1, m2 = binary["masses"]
    params = lal_dict(flags)
    approx = ls.GetApproximantFromString(approximant)
    return lambda: ls.SimInspiralChooseFDWaveformSequence(
        PHI_REF, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *binary["spins"], F_REF,
        DISTANCE * 1e6 * lal.PC_SI, INCLINATION, params, approx, frequencies,
    )


def uniform(approximant, binary, flags):
    m1, m2 = binary["masses"]
    params = lal_dict(flags)
    approx = ls.GetApproximantFromString(approximant)
    return lambda: ls.SimInspiralChooseFDWaveform(
        m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *binary["spins"],
        DISTANCE * 1e6 * lal.PC_SI, INCLINATION, PHI_REF, 0.0, 0.0, 0.0,
        0.25, 20.0, 1500.0, F_REF, params, approx,
    )


def nrsur(delta_t=1 / 4096):
    """NRSur7dq4 with SimInspiralChooseTDWaveform from 20 Hz."""
    m1, m2 = PRECESSING["masses"]
    return lambda: ls.SimInspiralChooseTDWaveform(
        m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *PRECESSING["spins"],
        DISTANCE * 1e6 * lal.PC_SI, INCLINATION, PHI_REF, 0.0, 0.0, 0.0,
        delta_t, F_REF, F_REF, None, ls.NRSur7dq4,
    )


def rust_timings():
    env = dict(os.environ, NRSUR7DQ4_DATA=str(ROOT / "lal_data" / "NRSur7dq4_v1.0.h5"))
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--release", "--example", "benchmark", "--", "--json"],
        capture_output=True, text=True, cwd=ROOT, check=True, env=env,
    )
    return {name: value["best"] for name, value in json.loads(result.stdout).items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT / "report" / "data" / "timings.json"))
    args = parser.parse_args()
    rust = rust_timings()
    lal_sequence = {name: best_time(sequence(*case)) for name, case in CASES.items()}
    lal_uniform = {name: best_time(uniform(*case)) for name, case in UNIFORM.items()}
    if "NRSur7dq4" in rust:
        lal_sequence["NRSur7dq4"] = best_time(nrsur())
        lal_sequence["NRSur7dq4 (workspace)"] = lal_sequence["NRSur7dq4"]
    out = dict(
        lalsimulation_version=ls.__version__,
        n_frequencies=len(FREQUENCIES),
        description=(
            "Best time per call (s), single thread. Rust: examples/benchmark.rs "
            "(construction + polarizations). LAL: SimInspiralChooseFDWaveformSequence "
            "on the same 5921 frequencies (multibanding off), and "
            "SimInspiralChooseFDWaveform on the uniform 0.25 Hz grid (stops at Mf = 0.3) "
            "with multibanding off and with the default thresholds."
        ),
        rust=rust,
        lal_sequence=lal_sequence,
        lal_uniform=lal_uniform,
    )
    path = Path(args.output)
    path.write_text(json.dumps(out, indent=1))
    print(f"{'model':<24}{'rust [us]':>12}{'LAL [us]':>12}{'speed-up':>10}")
    for name, t_rust in rust.items():
        t_lal = lal_sequence.get(name)
        if t_lal is not None:
            print(f"{name:<24}{t_rust * 1e6:12.1f}{t_lal * 1e6:12.1f}{t_lal / t_rust:10.2f}")
    for name, t_lal in lal_uniform.items():
        print(f"{name:<52} LAL {t_lal * 1e6:10.1f} us")


if __name__ == "__main__":
    main()
