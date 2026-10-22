#!/usr/bin/env python
"""
Compare the NRSur7dq4 polarizations of waveformrs with lalsimulation's
``SimInspiralChooseTDWaveform`` over random draws.

The draws span mass ratios 1 <= q <= 6 (with the lighter body first in
about a third of them, exercising the label switch), isotropic spins with
magnitudes up to 0.99, total masses 20-300 Msun (log-uniform) and sampling
rates of 1-8 kHz. ``f_min`` is drawn above the surrogate's starting
frequency for the binary (or zero, the whole surrogate, in 15% of the
draws) and ``f_ref`` is ``f_min`` or above it.

For each draw this records the maximum relative error ``max|h_rs - h_lal| /
max|h_lal|`` of both polarizations, whether the epochs and lengths agree,
the noise-weighted mismatch of each polarization (aLIGO design
sensitivity from 10 Hz, after an FFT of the zero-padded series, without
maximisation over time or phase) and the wall-clock times of both codes
(single thread; the Rust times reuse a workspace, as a sampler would).

Usage
-----
    OMP_NUM_THREADS=1 LAL_DATA_PATH=/path/to/lal_data \
        python scripts/validate_nrsur7dq4.py --n-draws 200 --outdir report/data

The Rust side is driven through
``cargo run --release --example nrsur7dq4_validation``.
"""
import argparse
import json
import subprocess
import time
from pathlib import Path

import lal
import lalsimulation as lalsim
import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def ordered(case):
    """Masses (Msun) and spins with the heavier body first, as LAL does."""
    m1 = case["mass_1"] * lal.MSUN_SI
    m2 = case["mass_2"] * lal.MSUN_SI
    s1, s2 = case["spin_1"], case["spin_2"]
    if m1 < m2:
        m1, m2 = m2, m1
        s1, s2 = [-s2[0], -s2[1], s2[2]], [-s1[0], -s1[1], s1[2]]
    return m1 / lal.MSUN_SI, m2 / lal.MSUN_SI, s1, s2


def start_frequency(case):
    """The surrogate's starting frequency for this binary, by bisection on
    f_ref: PrecessingNRSurDynamics fails below it."""
    mass_a, mass_b, s1, s2 = ordered(case)
    m_tot_sec = (mass_a + mass_b) * lal.MTSUN_SI

    def valid(f_ref):
        try:
            lalsim.PrecessingNRSurDynamics(
                mass_a / mass_b, *s1, *s2, (m_tot_sec * f_ref) * np.pi,
                1.0, 0.0, 0.0, 0.0, 0.0, None, lalsim.NRSur7dq4,
            )
            return True
        except RuntimeError:
            return False

    # omega_ref > 0.201 is rejected too: stay below it
    lo, hi = 0.1, 0.19 / (np.pi * m_tot_sec)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if valid(mid) else (mid, hi)
    return hi


def draw(rng):
    total_mass = float(np.exp(rng.uniform(np.log(20.0), np.log(300.0))))
    q = rng.uniform(1.0, 6.0)
    heavy, light = total_mass * q / (1 + q), total_mass / (1 + q)
    mass_1, mass_2 = (light, heavy) if rng.uniform() < 1 / 3 else (heavy, light)
    spins = list()
    for _ in range(2):
        magnitude = rng.uniform(0, 0.99)
        cos_theta = rng.uniform(-1, 1)
        phi = rng.uniform(0, 2 * np.pi)
        sin_theta = np.sqrt(1 - cos_theta**2)
        spins.append(
            (magnitude * np.array([sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta])).tolist()
        )
    rate = 4096.0 if total_mass < 60 else 2048.0 if total_mass < 150 else 1024.0
    if rng.uniform() < 0.25:
        rate *= 2
    case = dict(
        mass_1=float(mass_1),
        mass_2=float(mass_2),
        spin_1=spins[0],
        spin_2=spins[1],
        luminosity_distance=float(np.exp(rng.uniform(np.log(100), np.log(3000)))),
        inclination=float(np.arccos(rng.uniform(-1, 1))),
        phi_ref=float(rng.uniform(0, 2 * np.pi)),
        delta_t=1.0 / rate,
    )
    f_start = start_frequency(case)
    if rng.uniform() < 0.15:
        f_min = 0.0
        f_ref = 0.0 if rng.uniform() < 0.5 else f_start * rng.uniform(1.001, 1.5)
    else:
        f_min = f_start * rng.uniform(1.001, 1.6)
        f_ref = f_min if rng.uniform() < 0.5 else f_min * rng.uniform(1.0, 1.5)
    case.update(f_min=float(f_min), f_ref=float(f_ref), f_start=float(f_start))
    return case


def lal_waveform(case):
    args = (
        case["mass_1"] * lal.MSUN_SI,
        case["mass_2"] * lal.MSUN_SI,
        *case["spin_1"],
        *case["spin_2"],
        case["luminosity_distance"] * 1e6 * lal.PC_SI,
        case["inclination"],
        case["phi_ref"],
        0.0,
        0.0,
        0.0,
        case["delta_t"],
        case["f_min"],
        case["f_ref"],
        None,
        lalsim.NRSur7dq4,
    )
    try:
        hplus, hcross = lalsim.SimInspiralChooseTDWaveform(*args)
    except RuntimeError as error:
        return None, str(error)
    timings = list()
    for _ in range(3):
        start = time.perf_counter()
        lalsim.SimInspiralChooseTDWaveform(*args)
        timings.append(time.perf_counter() - start)
    return (
        dict(
            epoch=float(hplus.epoch),
            hplus=hplus.data.data.copy(),
            hcross=hcross.data.data.copy(),
            seconds=float(np.median(timings)),
        ),
        None,
    )


def rust_waveforms(cases):
    payload = json.dumps(dict(cases=cases))
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--release", "--example", "nrsur7dq4_validation"],
        input=payload,
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    )
    return json.loads(result.stdout)["results"]


def mismatch(a, b, delta_t, f_low=10.0):
    n = 1 << int(np.ceil(np.log2(2 * len(a))))
    fa = np.fft.rfft(a, n) * delta_t
    fb = np.fft.rfft(b, n) * delta_t
    frequencies = np.fft.rfftfreq(n, delta_t)
    keep = frequencies >= f_low
    psd = np.array([lalsim.SimNoisePSDaLIGOZeroDetHighPower(f) for f in frequencies[keep]])
    inner = lambda x, y: np.sum((x[keep] * np.conj(y[keep])).real / psd)
    aa, bb = inner(fa, fa), inner(fb, fb)
    if aa == 0 and bb == 0:
        return 0.0
    return float(1 - inner(fa, fb) / np.sqrt(aa * bb))


def max_relative_error(a, b):
    scale = np.max(np.abs(b))
    return float(np.max(np.abs(a - b)) / scale) if scale > 0 else 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-draws", type=int, default=200)
    parser.add_argument("--seed", type=int, default=150914)
    parser.add_argument("--outdir", default=str(ROOT / "report" / "data"))
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    cases = [draw(rng) for _ in range(args.n_draws)]
    rust = rust_waveforms(cases)
    records = list()
    failures = list()
    for case, rs in zip(cases, rust):
        lal_result, lal_error = lal_waveform(case)
        if rs["error"] is not None or lal_result is None:
            failures.append(dict(case, rust_error=rs["error"], lal_error=lal_error))
            continue
        hp, hc = np.array(rs["hplus"]), np.array(rs["hcross"])
        same_length = len(hp) == len(lal_result["hplus"])
        record = dict(
            case,
            n_samples=len(lal_result["hplus"]),
            same_length=same_length,
            same_epoch=rs["epoch"] == lal_result["epoch"],
            time_lal=lal_result["seconds"],
            time_rust=rs["seconds"],
            time_rust_fresh=rs["seconds_fresh"],
        )
        if same_length:
            record.update(
                max_rel_err_plus=max_relative_error(hp, lal_result["hplus"]),
                max_rel_err_cross=max_relative_error(hc, lal_result["hcross"]),
                bit_identical=bool(
                    np.array_equal(hp, lal_result["hplus"])
                    and np.array_equal(hc, lal_result["hcross"])
                ),
                mismatch_plus=mismatch(hp, lal_result["hplus"], case["delta_t"]),
                mismatch_cross=mismatch(hc, lal_result["hcross"], case["delta_t"]),
            )
        records.append(record)
    compared = [r for r in records if r["same_length"]]
    errors = np.array([max(r["max_rel_err_plus"], r["max_rel_err_cross"]) for r in compared])
    mismatches = np.array([max(r["mismatch_plus"], r["mismatch_cross"]) for r in compared])
    # Fresh calls: the reused-workspace timing repeats the same binary, so it
    # measures the cache for extrinsic-only changes, not a new waveform.
    speedup = np.array([r["time_lal"] / r["time_rust_fresh"] for r in records])
    summary = dict(
        approximant="NRSur7dq4",
        lalsimulation_version=lalsim.__version__,
        n_draws=len(cases),
        n_compared=len(compared),
        n_bit_identical=int(sum(r["bit_identical"] for r in compared)),
        n_epoch_mismatch=int(sum(not r["same_epoch"] for r in records)),
        n_length_mismatch=int(sum(not r["same_length"] for r in records)),
        failures=failures,
        max_relative_error=float(errors.max()),
        median_relative_error=float(np.median(errors)),
        max_mismatch=float(mismatches.max()),
        median_mismatch=float(np.median(mismatches)),
        median_time_lal=float(np.median([r["time_lal"] for r in records])),
        median_time_rust=float(np.median([r["time_rust_fresh"] for r in records])),
        median_time_rust_cached=float(np.median([r["time_rust"] for r in records])),
        median_speedup=float(np.median(speedup)),
        records=records,
    )
    with open(outdir / "validation_NRSur7dq4.json", "w") as ff:
        json.dump(summary, ff)
    print(
        f"NRSur7dq4: {len(compared)}/{len(cases)} compared, "
        f"{summary['n_bit_identical']} bit-identical, "
        f"max rel err {errors.max():.2e} (median {np.median(errors):.2e}), "
        f"max mismatch {mismatches.max():.2e}, "
        f"{summary['n_epoch_mismatch']} epoch / {summary['n_length_mismatch']} length "
        f"differences, {len(failures)} failures, "
        f"median time lal {1e3 * summary['median_time_lal']:.2f} ms / rust "
        f"{1e3 * summary['median_time_rust']:.2f} ms (x{np.median(speedup):.1f}; "
        f"{1e3 * summary['median_time_rust_cached']:.3f} ms when only extrinsic parameters change)"
    )


if __name__ == "__main__":
    main()
