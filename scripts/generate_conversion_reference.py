#!/usr/bin/env python
"""
Generate the LALSimulation reference data of the time/frequency-domain
conversions (``tests/conversions.rs``).

Suites (``tests/data/conversions_*.json``):

- ``conversions_td_<model>``: ``SimInspiralTD`` for the frequency-domain
  approximants IMRPhenomD, IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP,
  IMRPhenomXPHM (SpinTaylor angles) and TaylorF2: conditioned waveforms
  transformed to the time domain;
- ``conversions_fd_nrsur7dq4``: ``SimInspiralFD`` for NRSur7dq4 (time
  domain, transformed to the frequency domain);
- ``conversions_td_nrsur7dq4``: ``SimInspiralTD`` for NRSur7dq4 (tapered
  only, no Fourier transform);
- ``conversions_fd``: ``SimInspiralFD`` for frequency-domain approximants
  (conditioned in the frequency domain, no Fourier transform).

Each case stores the epoch (``gpsSeconds``, ``gpsNanoSeconds``), the
sampling interval, the full length and every ``stride``-th sample of the
polarizations, with ``stride`` chosen so that at most ``--max-samples``
samples are kept per polarization (the files stay below a few MB); the
Rust tests compute the full transforms and compare the stored samples.

Multibanding is switched off for IMRPhenomXHM and IMRPhenomXPHM
(``PhenomXHMThresholdMband = PhenomXPHMThresholdMband = 0``), which the
Rust port does not implement.

Usage
-----
    OMP_NUM_THREADS=1 LAL_DATA_PATH=/path/to/lal_data \\
        python scripts/generate_conversion_reference.py [--suites ...]

``OMP_NUM_THREADS=1`` is required: some IMRPhenomX(P)HM code paths are not
reproducible with several OpenMP threads.
"""
import argparse
import json
from pathlib import Path

import lal
import lalsimulation as lalsim
import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent / "tests" / "data"

NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0, PhenomXPHMThresholdMband=0)


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in (flags or dict()).items():
        getattr(lalsim, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def arguments(case):
    return (
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
    )


def downsample(case, data, max_samples):
    stride = max(1, -(-len(data) // max_samples))
    case["length"] = len(data)
    case["stride"] = stride
    return data[::stride]


def store_epoch(case, series):
    case["epoch_seconds"] = int(series.epoch.gpsSeconds)
    case["epoch_nanoseconds"] = int(series.epoch.gpsNanoSeconds)


def generate_td(case, max_samples):
    """SimInspiralTD (time-domain polarizations)."""
    hp, hc = lalsim.SimInspiralTD(
        *arguments(case),
        case["delta_t"],
        case["f_min"],
        case["f_ref"],
        lal_dict(case.get("lal_flags")),
        lalsim.GetApproximantFromString(case["approximant"]),
    )
    case = dict(case)
    store_epoch(case, hp)
    case["delta_t_out"] = hp.deltaT
    case["hplus"] = downsample(case, hp.data.data, max_samples).tolist()
    case["hcross"] = downsample(case, hc.data.data, max_samples).tolist()
    case["max_abs_hplus"] = float(np.max(np.abs(hp.data.data)))
    case["max_abs_hcross"] = float(np.max(np.abs(hc.data.data)))
    return case


def generate_fd(case, max_samples):
    """SimInspiralFD (frequency-domain polarizations)."""
    hp, hc = lalsim.SimInspiralFD(
        *arguments(case),
        case["delta_f"],
        case["f_min"],
        case["f_max"],
        case["f_ref"],
        lal_dict(case.get("lal_flags")),
        lalsim.GetApproximantFromString(case["approximant"]),
    )
    case = dict(case)
    store_epoch(case, hp)
    case["delta_f_out"] = hp.deltaF
    for name, series in [("hplus", hp), ("hcross", hc)]:
        data = downsample(case, series.data.data, max_samples)
        case[f"{name}_real"] = data.real.tolist()
        case[f"{name}_imag"] = data.imag.tolist()
        case[f"max_abs_{name}"] = float(np.max(np.abs(series.data.data)))
    return case


def binary(m1, m2, s1, s2, distance, iota, phi):
    return dict(
        mass_1=m1,
        mass_2=m2,
        spin_1=s1,
        spin_2=s2,
        luminosity_distance=distance,
        inclination=iota,
        phi_ref=phi,
    )


ALIGNED = [
    binary(36.0, 29.0, [0, 0, 0.3], [0, 0, -0.4], 500.0, 0.7, 1.2),
    binary(60.0, 15.0, [0, 0, 0.8], [0, 0, 0.5], 1000.0, 2.3, 4.0),
    binary(10.0, 5.0, [0, 0, -0.5], [0, 0, 0.2], 300.0, 1.4, 0.3),
    binary(150.0, 100.0, [0, 0, 0.9], [0, 0, 0.9], 3000.0, 2.9, 3.0),
    binary(20.0, 30.0, [0, 0, 0.1], [0, 0, -0.6], 400.0, 0.2, 5.5),
]

# f_min, f_ref, delta_t for the aligned binaries: f_ref = 0 (the models use
# their starting frequency, which the conditioning lowers), f_ref = f_min,
# f_ref > f_min, non-dyadic sampling, and f_min above the Kerr ISCO of the
# heavy binary (lowered to it)
ALIGNED_SAMPLING = [
    (20.0, 0.0, 1 / 4096),
    (15.0, 15.0, 1 / 4096),
    (20.0, 0.0, 1 / 1000),
    (10.5, 30.0, 1 / 1024),
    (25.0, 40.0, 1 / 2048),
]

PRECESSING = [
    binary(36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 500.0, 0.7, 1.2),
    binary(60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 1000.0, 2.3, 4.0),
    binary(12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 200.0, 0.3, 5.5),
    binary(30.0, 40.0, [0.5, 0.5, 0.5], [0.5, -0.5, -0.5], 400.0, 1.6, 2.0),
]

PRECESSING_SAMPLING = [
    (20.0, 0.0, 1 / 4096),
    (15.0, 20.0, 1 / 2048),
    (20.0, 20.0, 1 / 1000),
    (20.0, 0.0, 1 / 2048),
]


def td_cases(approximant, binaries, sampling, flags=None):
    return [
        dict(
            approximant=approximant,
            **b,
            f_min=f_min,
            f_ref=f_ref,
            delta_t=delta_t,
            lal_flags=flags or dict(),
        )
        for b, (f_min, f_ref, delta_t) in zip(binaries, sampling)
    ]


def taylorf2_cases():
    return [
        dict(
            approximant="TaylorF2",
            **binary(1.4, 1.3, [0, 0, 0.0], [0, 0, 0.0], 40.0, 0.4, 2.0),
            f_min=30.0,
            f_ref=0.0,
            delta_t=1 / 4096,
        ),
        dict(
            approximant="TaylorF2",
            **binary(10.0, 5.0, [0, 0, 0.3], [0, 0, -0.2], 300.0, 1.0, 0.5),
            f_min=30.0,
            f_ref=40.0,
            delta_t=1 / 3000,
        ),
        dict(
            approximant="TaylorF2",
            **binary(8.0, 1.4, [0, 0, 0.5], [0, 0, 0.0], 200.0, 2.0, 1.1),
            f_min=25.0,
            f_ref=0.0,
            delta_t=1 / 1000,
        ),
    ]


def nrsur_fd_cases():
    generic = binary(40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 100.0, 0.7, 0.1)
    configs = [
        # binary, f_min, f_ref, delta_f, f_max
        (generic, 20.0, 20.0, 0.0, 2048.0),
        (generic, 20.0, 20.0, 0.25, 2048.0),
        # f_max / delta_f not a power of two: f_max raised to 1024
        (generic, 20.0, 20.0, 0.25, 1000.0),
        # non-dyadic delta_f
        (generic, 20.0, 20.0, 0.1, 1024.0),
        # delta_f too large: the start of the waveform is cut
        (generic, 20.0, 20.0, 2.0, 1024.0),
        # the whole surrogate
        (binary(50.0, 30.0, [0.2, -0.6, 0.3], [0.5, 0.2, 0.1], 600.0, 2.0, 1.0), 0.0, 30.0, 0.0, 1024.0),
        # m1 < m2, high sampling rate
        (binary(20.0, 45.0, [0.3, 0.4, -0.2], [-0.5, 0.1, 0.6], 500.0, 1.4, 3.3), 20.0, 30.0, 0.0, 4096.0),
        (binary(96.0, 24.0, [0.7, -0.5, 0.4], [0.3, 0.3, -0.8], 1500.0, 1.2, 2.5), 12.0, 20.0, 0.125, 512.0),
    ]
    return [
        dict(approximant="NRSur7dq4", **b, f_min=f_min, f_ref=f_ref, delta_f=df, f_max=f_max)
        for b, f_min, f_ref, df, f_max in configs
    ]


def nrsur_td_cases():
    return [
        dict(
            approximant="NRSur7dq4",
            **binary(40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 100.0, 0.7, 0.1),
            f_min=20.0,
            f_ref=20.0,
            delta_t=1 / 4096,
        ),
        dict(
            approximant="NRSur7dq4",
            **binary(50.0, 30.0, [0.2, -0.6, 0.3], [0.5, 0.2, 0.1], 600.0, 2.0, 1.0),
            f_min=0.0,
            f_ref=30.0,
            delta_t=1 / 2048,
        ),
    ]


def fd_cases():
    """SimInspiralFD for frequency-domain approximants."""
    aligned = ALIGNED[0]
    precessing = PRECESSING[0]
    return [
        dict(approximant="IMRPhenomD", **aligned, f_min=20.0, f_ref=0.0, delta_f=0.0, f_max=2048.0),
        dict(approximant="IMRPhenomD", **aligned, f_min=20.0, f_ref=30.0, delta_f=0.25, f_max=1000.0),
        dict(approximant="IMRPhenomXAS", **aligned, f_min=20.0, f_ref=0.0, delta_f=0.125, f_max=1024.0),
        dict(
            approximant="IMRPhenomXHM",
            **ALIGNED[1],
            f_min=15.0,
            f_ref=0.0,
            delta_f=0.0,
            f_max=1024.0,
            lal_flags=dict(NO_MULTIBANDING),
        ),
        dict(
            approximant="IMRPhenomXPHM",
            **precessing,
            f_min=20.0,
            f_ref=0.0,
            delta_f=0.0,
            f_max=2048.0,
            lal_flags=dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
        ),
        dict(approximant="TaylorF2", **ALIGNED[2], f_min=30.0, f_ref=0.0, delta_f=0.0, f_max=1024.0),
    ]


SUITES = dict(
    conversions_td_imrphenomd=(generate_td, lambda: td_cases("IMRPhenomD", ALIGNED, ALIGNED_SAMPLING)),
    conversions_td_imrphenomxas=(generate_td, lambda: td_cases("IMRPhenomXAS", ALIGNED, ALIGNED_SAMPLING)),
    conversions_td_imrphenomxhm=(
        generate_td,
        lambda: td_cases("IMRPhenomXHM", ALIGNED, ALIGNED_SAMPLING, NO_MULTIBANDING),
    ),
    conversions_td_imrphenomxp=(
        generate_td,
        lambda: td_cases("IMRPhenomXP", PRECESSING[:2], PRECESSING_SAMPLING[:2], dict(PhenomXPrecVersion=320)),
    ),
    conversions_td_imrphenomxphm=(
        generate_td,
        lambda: td_cases(
            "IMRPhenomXPHM",
            PRECESSING,
            PRECESSING_SAMPLING,
            dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
        ),
    ),
    conversions_td_taylorf2=(generate_td, taylorf2_cases),
    conversions_fd_nrsur7dq4=(generate_fd, nrsur_fd_cases),
    conversions_td_nrsur7dq4=(generate_td, nrsur_td_cases),
    conversions_fd=(generate_fd, fd_cases),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suites", nargs="*", default=list(SUITES))
    parser.add_argument("--max-samples", type=int, default=4096)
    args = parser.parse_args()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name in args.suites:
        generate, cases = SUITES[name]
        cases = [generate(case, args.max_samples) for case in cases()]
        output = dict(
            lalsimulation_version=lalsim.__version__,
            lal_version=lal.__version__,
            cases=cases,
        )
        path = DATA_DIR / f"{name}.json"
        with open(path, "w") as ff:
            json.dump(output, ff)
        print(f"Wrote {len(cases)} cases to {path} ({path.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
