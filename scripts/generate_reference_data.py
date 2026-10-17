#!/usr/bin/env python
"""
Generate reference frequency-domain waveforms from LALSimulation.

The output JSON files in ``tests/data`` are consumed by the Rust integration
tests (``tests/lalsimulation.rs``).  Each file contains a list of cases, each
case holding the input parameters, a frequency array and the plus/cross
polarizations computed with ``SimInspiralChooseFDWaveformSequence`` (or with
``SimInspiralChooseFDWaveform`` on a uniform grid starting at zero frequency
for cases with ``delta_f > 0``).  The ``imrphenomxp_spintaylor_angles`` suite
instead stores the Euler angles of ``SimIMRPhenomXPSpinTaylorAngles``.

Usage
-----
    python scripts/generate_reference_data.py [--n-frequencies 256]

Requires ``python-lalsimulation`` (tested with lalsimulation 6.2.0).
"""
import argparse
import json
from pathlib import Path

import lal
import lalsimulation as lalsim
import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent / "tests" / "data"


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in (flags or dict()).items():
        getattr(lalsim, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def generate_uniform(case):
    """Polarizations from SimInspiralChooseFDWaveform on the grid k * delta_f."""
    s1 = case["spin_1"]
    s2 = case["spin_2"]
    hplus, hcross = lalsim.SimInspiralChooseFDWaveform(
        case["mass_1"] * lal.MSUN_SI,
        case["mass_2"] * lal.MSUN_SI,
        *s1,
        *s2,
        case["luminosity_distance"] * 1e6 * lal.PC_SI,
        case["inclination"],
        case["phi_ref"],
        0.0,
        0.0,
        0.0,
        case["delta_f"],
        case["f_min"],
        case["f_max"],
        case["f_ref"],
        lal_dict(case.get("lal_flags")),
        lalsim.GetApproximantFromString(case["approximant"]),
    )
    case = dict(case)
    case["frequencies"] = (np.arange(hplus.data.length) * case["delta_f"]).tolist()
    case["hplus_real"] = hplus.data.data.real.tolist()
    case["hplus_imag"] = hplus.data.data.imag.tolist()
    case["hcross_real"] = hcross.data.data.real.tolist()
    case["hcross_imag"] = hcross.data.data.imag.tolist()
    return case


def initialise_phenomx_globals():
    """XLALSimIMRPhenomXPSpinTaylorAngles relies on the global powers_of_lalpi
    that only the IMRPhenomX waveform drivers initialise: without a previous
    waveform call the 3PN orbital angular momentum misses its pi^2 term."""
    lalsim.SimInspiralChooseFDWaveform(
        30 * lal.MSUN_SI, 30 * lal.MSUN_SI, 0, 0, 0, 0, 0, 0, 1e6 * lal.PC_SI,
        0, 0, 0, 0, 0, 1.0, 20.0, 64.0, 20.0, None, lalsim.IMRPhenomXAS,
    )


def generate_spintaylor_angles(case):
    """Euler angles from SimIMRPhenomXPSpinTaylorAngles, every stride-th point."""
    initialise_phenomx_globals()
    alpha, cosbeta, gamma = lalsim.SimIMRPhenomXPSpinTaylorAngles(
        case["mass_1"] * lal.MSUN_SI,
        case["mass_2"] * lal.MSUN_SI,
        *case["spin_1"],
        *case["spin_2"],
        case["f_min"],
        case["f_max"],
        case["delta_f"],
        case["f_ref"],
        case["phi_ref"],
        lal_dict(dict(PhenomXPrecVersion=320)),
    )
    stride = case["stride"]
    case = dict(case)
    case["length"] = alpha.length
    case["alpha"] = alpha.data[::stride].tolist()
    case["cos_beta"] = cosbeta.data[::stride].tolist()
    case["gamma"] = gamma.data[::stride].tolist()
    return case


def generate(case, n_frequencies):
    if "stride" in case:
        return generate_spintaylor_angles(case)
    if case.get("delta_f", 0.0) > 0.0:
        return generate_uniform(case)
    f_min = case["f_min"]
    f_max = case["f_max"]
    frequencies = np.geomspace(f_min, f_max, n_frequencies)
    if case["f_ref"] > 0:
        frequencies = np.unique(np.concatenate([frequencies, [case["f_ref"]]]))
    freq_vector = lal.CreateREAL8Vector(len(frequencies))
    freq_vector.data = frequencies
    s1 = case.get("spin_1", [0.0, 0.0, case.get("chi_1", 0.0)])
    s2 = case.get("spin_2", [0.0, 0.0, case.get("chi_2", 0.0)])
    hplus, hcross = lalsim.SimInspiralChooseFDWaveformSequence(
        case["phi_ref"],
        case["mass_1"] * lal.MSUN_SI,
        case["mass_2"] * lal.MSUN_SI,
        *s1,
        *s2,
        case["f_ref"],
        case["luminosity_distance"] * 1e6 * lal.PC_SI,
        case["inclination"],
        lal_dict(case.get("lal_flags")),
        lalsim.GetApproximantFromString(case["approximant"]),
        freq_vector,
    )
    case = dict(case)
    case["spin_1"] = list(s1)
    case["spin_2"] = list(s2)
    case["frequencies"] = frequencies.tolist()
    case["hplus_real"] = hplus.data.data.real.tolist()
    case["hplus_imag"] = hplus.data.data.imag.tolist()
    case["hcross_real"] = hcross.data.data.real.tolist()
    case["hcross_imag"] = hcross.data.data.imag.tolist()
    return case


def aligned_cases(approximant, f_max_fn, lal_flags=None):
    """A fixed set of aligned-spin configurations spanning the parameter space."""
    configs = [
        # mass_1, mass_2, chi_1, chi_2, distance, inclination, phi_ref, f_min, f_ref
        (30.0, 30.0, 0.0, 0.0, 400.0, 0.0, 0.0, 20.0, 20.0),
        (36.0, 29.0, 0.3, -0.4, 500.0, 0.7, 1.2, 20.0, 30.0),
        (60.0, 15.0, 0.8, 0.5, 1000.0, 2.3, 4.0, 15.0, 20.0),
        (80.0, 10.0, -0.7, 0.2, 2000.0, 1.0, 0.3, 10.0, 25.0),
        (1.6, 1.3, 0.05, 0.02, 40.0, 0.4, 2.0, 30.0, 40.0),
        (8.0, 1.4, 0.9, 0.0, 200.0, 1.5, 5.5, 20.0, 50.0),
        (150.0, 100.0, 0.95, 0.95, 3000.0, 2.9, 3.0, 10.0, 15.0),
        (25.0, 20.0, -0.9, -0.9, 300.0, 0.2, 0.7, 20.0, 20.0),
    ]
    cases = list()
    for m1, m2, c1, c2, dl, iota, phi, fmin, fref in configs:
        cases.append(
            dict(
                approximant=approximant,
                mass_1=m1,
                mass_2=m2,
                chi_1=c1,
                chi_2=c2,
                luminosity_distance=dl,
                inclination=iota,
                phi_ref=phi,
                f_min=fmin,
                f_ref=fref,
                f_max=f_max_fn(m1 + m2),
                lal_flags=lal_flags or dict(),
            )
        )
    return cases


def precessing_cases(approximant, f_max_fn, lal_flags=None):
    configs = [
        # mass_1, mass_2, spin_1, spin_2, distance, inclination, phi_ref, f_min, f_ref
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 500.0, 0.7, 1.2, 20.0, 20.0),
        (60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 1000.0, 2.3, 4.0, 15.0, 20.0),
        (80.0, 10.0, [-0.7, 0.4, -0.2], [0.0, 0.5, 0.2], 2000.0, 1.0, 0.3, 10.0, 25.0),
        (30.0, 30.0, [0.5, 0.5, 0.5], [0.5, -0.5, -0.5], 400.0, 1.6, 2.0, 20.0, 30.0),
        (12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 200.0, 0.3, 5.5, 20.0, 20.0),
        (150.0, 100.0, [0.1, 0.8, 0.4], [0.2, 0.0, -0.8], 3000.0, 2.9, 3.0, 10.0, 15.0),
        (40.0, 4.0, [0.8, 0.0, 0.5], [0.0, 0.0, 0.0], 800.0, 1.2, 0.9, 15.0, 20.0),
        # effectively aligned spins: lalsimulation falls back to NNLO angles
        (25.0, 20.0, [0.0, 0.0, 0.5], [0.0, 0.0, -0.3], 300.0, 0.9, 0.7, 20.0, 20.0),
    ]
    cases = list()
    for m1, m2, s1, s2, dl, iota, phi, fmin, fref in configs:
        cases.append(
            dict(
                approximant=approximant,
                mass_1=m1,
                mass_2=m2,
                spin_1=s1,
                spin_2=s2,
                luminosity_distance=dl,
                inclination=iota,
                phi_ref=phi,
                f_min=fmin,
                f_ref=fref,
                f_max=f_max_fn(m1 + m2),
                lal_flags=lal_flags or dict(),
            )
        )
    return cases


def uniform_precessing_cases():
    """IMRPhenomXP on uniform grids (SimInspiralChooseFDWaveform)."""
    configs = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, delta_f, f_max, PrecVersion
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 20.0, 1.0, 0.0, 320),
        (60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0, 0.5, 512.0, 320),
        (12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 20.0, 0.0, 1.0, 1024.0, 320),
        (12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 20.0, 25.0, 1.0, 1024.0, 102),
    ]
    return [
        dict(
            approximant="IMRPhenomXP",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=500.0,
            inclination=1.1,
            phi_ref=0.6,
            f_min=fmin,
            f_ref=fref,
            delta_f=df,
            f_max=fmax,
            lal_flags=dict(PhenomXPrecVersion=version),
        )
        for m1, m2, s1, s2, fmin, fref, df, fmax, version in configs
    ]


def spintaylor_angle_cases():
    """SimIMRPhenomXPSpinTaylorAngles (PrecVersion 320)."""
    configs = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_max, delta_f, f_ref, phi_ref
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 1024.0, 0.125, 20.0, 1.2),
        (80.0, 10.0, [-0.7, 0.4, -0.2], [0.0, 0.5, 0.2], 10.0, 600.0, 0.0625, 25.0, 0.3),
        (12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 20.0, 2048.0, 0.125, 20.0, 5.5),
        (25.0, 25.0, [0.0, 0.6, 0.1], [-0.3, 0.0, 0.4], 15.0, 800.0, 0.25, 30.0, 2.0),
    ]
    return [
        dict(
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            f_min=fmin,
            f_max=fmax,
            delta_f=df,
            f_ref=fref,
            phi_ref=phi,
            stride=16,
        )
        for m1, m2, s1, s2, fmin, fmax, df, fref, phi in configs
    ]


def f_max_isco(total_mass):
    """Schwarzschild ISCO GW frequency, TaylorF2's default cutoff."""
    return 1.0 / (6.0**1.5 * np.pi * total_mass * lal.MTSUN_SI)


def f_max_imr(mf_max):
    def func(total_mass):
        return min(mf_max / (total_mass * lal.MTSUN_SI), 4096.0)

    return func


NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0.0, PhenomXPHMThresholdMband=0.0)

SUITES = dict(
    taylorf2=lambda: aligned_cases("TaylorF2", f_max_isco)
    + [
        dict(case, f_ref=0.0)
        for case in aligned_cases("TaylorF2", f_max_isco)[1::2]
    ],
    imrphenomd=lambda: aligned_cases("IMRPhenomD", f_max_imr(0.2)),
    imrphenomxas=lambda: aligned_cases("IMRPhenomXAS", f_max_imr(0.3)),
    imrphenomxhm=lambda: aligned_cases(
        "IMRPhenomXHM", f_max_imr(0.3), lal_flags=NO_MULTIBANDING
    ),
    imrphenomxp_nnlo=lambda: precessing_cases(
        "IMRPhenomXP", f_max_imr(0.3), lal_flags=dict(PhenomXPrecVersion=102)
    ),
    imrphenomxp_spintaylor=lambda: precessing_cases(
        "IMRPhenomXP", f_max_imr(0.3), lal_flags=dict(PhenomXPrecVersion=320)
    ),
    imrphenomxp_uniform=uniform_precessing_cases,
    imrphenomxp_spintaylor_angles=spintaylor_angle_cases,
    imrphenomxphm_spintaylor=lambda: precessing_cases(
        "IMRPhenomXPHM",
        f_max_imr(0.3),
        lal_flags=dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
    ),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-frequencies", type=int, default=256)
    parser.add_argument("--suites", nargs="*", default=list(SUITES))
    args = parser.parse_args()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name in args.suites:
        cases = [generate(case, args.n_frequencies) for case in SUITES[name]()]
        output = dict(
            lalsimulation_version=lalsim.__version__,
            lal_version=lal.__version__,
            cases=cases,
        )
        path = DATA_DIR / f"{name}.json"
        with open(path, "w") as ff:
            json.dump(output, ff)
        print(f"Wrote {len(cases)} cases to {path}")


if __name__ == "__main__":
    main()
