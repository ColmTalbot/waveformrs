#!/usr/bin/env python
"""
Throughput when generating many waveforms one after another on a fixed
frequency grid, as in a stochastic sampler: random binaries, the same grid
for every call, single thread.

Writes ``report/data/serial_timings.json`` with the mean time per waveform
for waveformrs (``examples/serial_benchmark.rs``) and lalsimulation
(``SimInspiralChooseFDWaveformSequence``, multibanding off) on the same
draws.

    OMP_NUM_THREADS=1 LAL_DATA_PATH=lal_data python report/serial_benchmark.py [--n-draws 300]

NRSur7dq4 is included when ``lal_data/NRSur7dq4_v1.0.h5`` exists.
"""
import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

import lal
import lalsimulation as ls
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
# NRSur7dq4 data (git-ignored), for both codes; LAL reads LAL_DATA_PATH at import
LAL_DATA = ROOT / "lal_data"
ENV = dict(os.environ, NRSUR7DQ4_DATA=str(LAL_DATA / "NRSur7dq4_v1.0.h5"))
os.environ.setdefault("LAL_DATA_PATH", str(LAL_DATA))
NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0.0, PhenomXPHMThresholdMband=0.0)
MODELS = {
    "TaylorF2": ("TaylorF2", dict()),
    "IMRPhenomD": ("IMRPhenomD", dict()),
    "IMRPhenomXAS": ("IMRPhenomXAS", dict()),
    "IMRPhenomXHM": ("IMRPhenomXHM", NO_MULTIBANDING),
    "IMRPhenomXP (102)": ("IMRPhenomXP", dict(PhenomXPrecVersion=102)),
    "IMRPhenomXP (223)": ("IMRPhenomXP", dict(PhenomXPrecVersion=223)),
    "IMRPhenomXP (320)": ("IMRPhenomXP", dict(PhenomXPrecVersion=320)),
    "IMRPhenomXPHM (102)": ("IMRPhenomXPHM", dict(PhenomXPrecVersion=102, **NO_MULTIBANDING)),
    "IMRPhenomXPHM (223)": ("IMRPhenomXPHM", dict(PhenomXPrecVersion=223, **NO_MULTIBANDING)),
    "IMRPhenomXPHM (320)": ("IMRPhenomXPHM", dict(PhenomXPrecVersion=320, **NO_MULTIBANDING)),
    # (SimInspiralChooseFDWaveformSequence switches multibanding off; the
    # Rust side without the parity fix of lalsuite !2593, as in 6.2.0)
    "IMRPhenomXPNR": ("IMRPhenomXPNR", dict()),
}
ALIGNED = {"TaylorF2", "IMRPhenomD", "IMRPhenomXAS", "IMRPhenomXHM"}


def draws(n, seed):
    """Binary black holes as in a broad BBH prior; spins isotropic, |chi| <= 0.9."""
    rng = np.random.default_rng(seed)
    out = list()
    for _ in range(n):
        total_mass = float(np.exp(rng.uniform(np.log(20.0), np.log(150.0))))
        mass_ratio = float(rng.uniform(1 / 6, 1.0))
        mass_1 = total_mass / (1 + mass_ratio)
        spins = list()
        for _ in range(2):
            magnitude = rng.uniform(0, 0.9)
            cos_theta = rng.uniform(-1, 1)
            phi = rng.uniform(0, 2 * np.pi)
            sin_theta = np.sqrt(1 - cos_theta**2)
            spins.append([float(x) for x in magnitude * np.array(
                [sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta])])
        out.append(dict(
            mass_1=mass_1, mass_2=total_mass - mass_1, spin_1=spins[0], spin_2=spins[1],
            luminosity_distance=float(np.exp(rng.uniform(np.log(100), np.log(3000)))),
            inclination=float(np.arccos(rng.uniform(-1, 1))),
            phi_ref=float(rng.uniform(0, 2 * np.pi)), f_ref=20.0,
        ))
    return out


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in flags.items():
        getattr(ls, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def lal_serial(name, cases, frequencies):
    approximant, flags = MODELS[name]
    approx = ls.GetApproximantFromString(approximant)
    params = lal_dict(flags)
    fv = lal.CreateREAL8Vector(len(frequencies))
    fv.data = frequencies
    aligned = name in ALIGNED

    def call(case):
        s1 = [0.0, 0.0, case["spin_1"][2]] if aligned else case["spin_1"]
        s2 = [0.0, 0.0, case["spin_2"][2]] if aligned else case["spin_2"]
        return ls.SimInspiralChooseFDWaveformSequence(
            case["phi_ref"], case["mass_1"] * lal.MSUN_SI, case["mass_2"] * lal.MSUN_SI,
            *s1, *s2, case["f_ref"], case["luminosity_distance"] * 1e6 * lal.PC_SI,
            case["inclination"], params, approx, fv,
        )

    failures = 0
    for case in cases[:5]:
        try:
            call(case)
        except RuntimeError:
            pass
    start = time.perf_counter()
    for case in cases:
        try:
            call(case)
        except RuntimeError:
            failures += 1
    return dict(seconds=(time.perf_counter() - start) / len(cases), failures=failures)


def lal_nrsur(cases, delta_t=1 / 4096):
    """NRSur7dq4 with SimInspiralChooseTDWaveform from the start of the
    surrogate (f_min = f_ref = 0)."""
    def call(case):
        return ls.SimInspiralChooseTDWaveform(
            case["mass_1"] * lal.MSUN_SI, case["mass_2"] * lal.MSUN_SI,
            *case["spin_1"], *case["spin_2"], case["luminosity_distance"] * 1e6 * lal.PC_SI,
            case["inclination"], case["phi_ref"], 0.0, 0.0, 0.0, delta_t, 0.0, 0.0,
            None, ls.NRSur7dq4,
        )

    failures = 0
    for case in cases[:5]:
        try:
            call(case)
        except RuntimeError:
            pass
    start = time.perf_counter()
    for case in cases:
        try:
            call(case)
        except RuntimeError:
            failures += 1
    return dict(seconds=(time.perf_counter() - start) / len(cases), failures=failures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-draws", type=int, default=300)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--delta-f", type=float, default=0.125)
    parser.add_argument("--f-min", type=float, default=20.0)
    parser.add_argument("--f-max", type=float, default=1024.0)
    parser.add_argument("--output", default=str(ROOT / "report" / "data" / "serial_timings.json"))
    args = parser.parse_args()

    frequencies = np.arange(args.f_min, args.f_max + args.delta_f / 2, args.delta_f)
    cases = draws(args.n_draws, args.seed)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "draws.json"
        path.write_text(json.dumps(dict(frequencies=frequencies.tolist(), draws=cases)))
        result = subprocess.run(
            ["cargo", "run", "--quiet", "--release", "--example", "serial_benchmark", "--", str(path)],
            cwd=ROOT, capture_output=True, text=True, check=True, env=ENV,
        )
    rust = json.loads(result.stdout)
    lal_results = {name: lal_serial(name, cases, frequencies) for name in MODELS}
    names = list(MODELS)
    if "NRSur7dq4" in rust:
        lal_results["NRSur7dq4"] = lal_nrsur(cases)
        names.append("NRSur7dq4")

    # Cost of a call through the Python bindings with (almost) no work.
    tiny = np.array([args.f_min, args.f_min + args.delta_f])
    overhead = lal_serial("TaylorF2", cases, tiny)["seconds"]

    out = dict(
        lalsimulation_version=ls.__version__,
        n_draws=len(cases),
        n_frequencies=len(frequencies),
        grid=dict(f_min=args.f_min, f_max=args.f_max, delta_f=args.delta_f),
        lal_python_call_overhead=overhead,
        rust=rust,
        lal=lal_results,
    )
    Path(args.output).write_text(json.dumps(out, indent=1))
    print(f"{len(cases)} draws, {len(frequencies)} frequencies; LAL call overhead "
          f"{overhead * 1e6:.1f} us")
    print(f"{'model':<22}{'rust [us]':>12}{'grid [us]':>12}{'LAL [us]':>12}"
          f"{'speed-up':>10}{'grid':>8}")
    for name in names:
        t_rust, t_lal = rust[name]["seconds"], lal_results[name]["seconds"]
        t_grid = rust[name].get("grid_seconds")
        grid_text = f"{t_grid * 1e6:12.1f}" if t_grid else f"{'-':>12}"
        grid_speed = f"{t_lal / t_grid:8.2f}" if t_grid else f"{'-':>8}"
        print(f"{name:<22}{t_rust * 1e6:12.1f}{grid_text}{t_lal * 1e6:12.1f}"
              f"{t_lal / t_rust:10.2f}{grid_speed}"
              f"   errors rust/lal {rust[name]['failures']}/{lal_results[name]['failures']}")


if __name__ == "__main__":
    main()
