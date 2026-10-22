#!/usr/bin/env python
"""
Reference data from LALSimulation for long, low-mass waveforms.

The cases have up to ~520k frequencies (Delta f = 1/128 Hz from 20 Hz to
4 kHz and similar), and binaries whose signal lasts up to ~130 s from
f_min (a few last for several minutes). Storing whole waveforms would take
hundreds of MB, so each case stores the definition of its grid together with
LAL's polarizations at a subset of the indices only: every k-th point, the
first and last points of the band, and log-spaced points that resolve the
low frequencies, where the inspiral amplitude is largest. The Rust tests
(``tests/long_waveforms.rs``) evaluate the whole grid, since the SpinTaylor
angles depend on it, and compare at the stored indices.

Each case is generated with one of

- ``sequence``: ``SimInspiralChooseFDWaveformSequence`` on the frequencies
  ``f_min + k delta_f``, ``k = 0 .. length - 1``;
- ``uniform``: ``SimInspiralChooseFDWaveform`` (multibanding off) on its
  grid ``k delta_f``, ``k = 0 .. length - 1``.

Multibanding is switched off for IMRPhenomXHM and IMRPhenomXPHM (it is not
used for frequency sequences anyway).

The ``angles`` entries hold every ``stride``-th value of the SpinTaylor Euler
angles of ``SimIMRPhenomXPSpinTaylorAngles`` for a long integration of the
PN equations.

Usage
-----
    OMP_NUM_THREADS=1 python scripts/generate_long_reference.py [--points 1024]

``OMP_NUM_THREADS=1`` is required: the MSA angles of lalsimulation are not
reproducible with several OpenMP threads (see
``report/lalsimulation-openmp-race.md``).

Requires ``python-lalsimulation`` (tested with lalsimulation 6.2.0).
"""
import argparse
import json
import os
import time
from pathlib import Path

import lal
import lalsimulation as lalsim
import numpy as np

OUTPUT = Path(__file__).resolve().parent.parent / "tests" / "data" / "long_waveforms.json"
NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0.0, PhenomXPHMThresholdMband=0.0)

# Binaries. Durations are XLALSimInspiralChirpTimeBound from f_min, an upper
# bound (the signals are ~15% shorter).
ALIGNED = [
    # 3 Msun, ~150 s from 20 Hz: the 128-s segment of a binary neutron star
    dict(name="bns_20Hz_4kHz", mass_1=1.6, mass_2=1.4, chi_1=0.05, chi_2=0.02,
         f_min=20.0, f_max=4096.0, delta_f=1 / 128, f_ref=20.0),
    # 6 Msun, q = 1/2, ~55 s from 20 Hz
    dict(name="m6_q05_20Hz_2kHz", mass_1=4.0, mass_2=2.0, chi_1=0.3, chi_2=-0.2,
         f_min=20.0, f_max=2048.0, delta_f=1 / 128, f_ref=35.0),
    # 2.6 Msun from 10 Hz, ~1200 s: a longer inspiral (phase ~6e4 rad)
    dict(name="bns_10Hz_512Hz", mass_1=1.4, mass_2=1.2, chi_1=-0.03, chi_2=0.04,
         f_min=10.0, f_max=512.0, delta_f=1 / 1024, f_ref=10.0),
    # 10 Msun, q = 1/4, strong aligned spin, ~25 s: a shorter case
    dict(name="m10_q025_20Hz_2kHz", mass_1=8.0, mass_2=2.0, chi_1=0.7, chi_2=0.1,
         f_min=20.0, f_max=2048.0, delta_f=1 / 32, f_ref=20.0),
]

PRECESSING = [
    # 3.1 Msun, ~130 s from 20 Hz, reference at f_min
    dict(name="p31_20Hz_2kHz", mass_1=1.8, mass_2=1.3,
         spin_1=[0.3, 0.1, 0.05], spin_2=[-0.1, 0.2, 0.02],
         f_min=20.0, f_max=2048.0, delta_f=1 / 128, f_ref=20.0),
    # 6 Msun, q = 1/2, larger in-plane spins, reference above f_min
    dict(name="p6_q05_20Hz_1kHz", mass_1=4.0, mass_2=2.0,
         spin_1=[0.5, -0.3, 0.2], spin_2=[0.1, 0.4, -0.3],
         f_min=20.0, f_max=1024.0, delta_f=1 / 128, f_ref=30.0),
]

# SimIMRPhenomXPSpinTaylorAngles from 10 Hz with f_ref = 20 Hz: the PN
# equations are integrated backward over ~900 s (~2.3 million output
# samples of the coarse step) and forward through merger. Every stride-th
# value of the ~260k angles is stored.
ANGLES = [
    dict(name="p31_angles_10Hz_2kHz", mass_1=1.8, mass_2=1.3,
         spin_1=[0.3, 0.1, 0.05], spin_2=[-0.1, 0.2, 0.02],
         f_min=10.0, f_max=2048.0, delta_f=1 / 128, f_ref=20.0, phi_ref=0.4, stride=256),
]

ALIGNED_MODELS = ["TaylorF2", "IMRPhenomD", "IMRPhenomXAS", "IMRPhenomXHM"]
PRECESSING_MODELS = ["IMRPhenomXP", "IMRPhenomXPHM"]
VERSIONS = [102, 223, 320]


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in flags.items():
        getattr(lalsim, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def subset(length, first, last, points):
    """Every k-th index of [first, last] (about `points` of them), the first
    and last 16 indices of the band and `points // 4` log-spaced indices from
    `first`."""
    stride = max(1, (last - first) // points)
    indices = set(range(first, last + 1, stride))
    indices.update(range(first, min(first + 16, last + 1)))
    indices.update(range(max(first, last - 15), last + 1))
    span = last - first + 1
    indices.update(int(first + x) for x in np.geomspace(1, span, points // 4) - 1)
    return sorted(i for i in indices if first <= i <= last)


def cases():
    out = list()
    for binary in ALIGNED:
        for approximant in ALIGNED_MODELS:
            flags = NO_MULTIBANDING if approximant == "IMRPhenomXHM" else dict()
            out.append(dict(
                binary, approximant=approximant, generator="sequence", lal_flags=flags,
                spin_1=[0.0, 0.0, binary["chi_1"]], spin_2=[0.0, 0.0, binary["chi_2"]],
                luminosity_distance=150.0, inclination=0.6, phi_ref=1.3,
            ))
    for binary in PRECESSING:
        for approximant in PRECESSING_MODELS:
            for version in VERSIONS:
                flags = dict(PhenomXPrecVersion=version)
                if approximant == "IMRPhenomXPHM":
                    flags.update(NO_MULTIBANDING)
                for generator in ["sequence", "uniform"]:
                    out.append(dict(
                        binary, approximant=approximant, generator=generator,
                        lal_flags=flags, luminosity_distance=200.0,
                        inclination=1.1, phi_ref=0.4,
                    ))
    for case in out:
        case.pop("chi_1", None)
        case.pop("chi_2", None)
    return out


def generate(case, points):
    approximant = lalsim.GetApproximantFromString(case["approximant"])
    params = lal_dict(case["lal_flags"])
    m1 = case["mass_1"] * lal.MSUN_SI
    m2 = case["mass_2"] * lal.MSUN_SI
    distance = case["luminosity_distance"] * 1e6 * lal.PC_SI
    start = time.perf_counter()
    if case["generator"] == "sequence":
        length = int(round((case["f_max"] - case["f_min"]) / case["delta_f"])) + 1
        frequencies = case["f_min"] + np.arange(length) * case["delta_f"]
        freq_vector = lal.CreateREAL8Vector(length)
        freq_vector.data = frequencies
        hplus, hcross = lalsim.SimInspiralChooseFDWaveformSequence(
            case["phi_ref"], m1, m2, *case["spin_1"], *case["spin_2"], case["f_ref"],
            distance, case["inclination"], params, approximant, freq_vector,
        )
    else:
        hplus, hcross = lalsim.SimInspiralChooseFDWaveform(
            m1, m2, *case["spin_1"], *case["spin_2"], distance, case["inclination"],
            case["phi_ref"], 0.0, 0.0, 0.0, case["delta_f"], case["f_min"], case["f_max"],
            case["f_ref"], params, approximant,
        )
        length = hplus.data.length
    seconds = time.perf_counter() - start
    hp = hplus.data.data
    hc = hcross.data.data
    assert len(hp) == length
    nonzero = np.flatnonzero(np.abs(hp) + np.abs(hc))
    indices = subset(length, int(nonzero[0]), int(nonzero[-1]), points)
    case = dict(case)
    case["length"] = length
    case["indices"] = indices
    case["hplus_real"] = hp.real[indices].tolist()
    case["hplus_imag"] = hp.imag[indices].tolist()
    case["hcross_real"] = hc.real[indices].tolist()
    case["hcross_imag"] = hc.imag[indices].tolist()
    return case, seconds


def initialise_phenomx_globals():
    """XLALSimIMRPhenomXPSpinTaylorAngles relies on the global powers_of_lalpi
    that only the IMRPhenomX waveform drivers initialise (see
    scripts/generate_reference_data.py)."""
    lalsim.SimInspiralChooseFDWaveform(
        30 * lal.MSUN_SI, 30 * lal.MSUN_SI, 0, 0, 0, 0, 0, 0, 1e6 * lal.PC_SI,
        0, 0, 0, 0, 0, 1.0, 20.0, 64.0, 20.0, None, lalsim.IMRPhenomXAS,
    )


def generate_angles(case):
    """Euler angles from SimIMRPhenomXPSpinTaylorAngles, every stride-th value."""
    initialise_phenomx_globals()
    start = time.perf_counter()
    alpha, cos_beta, gamma = lalsim.SimIMRPhenomXPSpinTaylorAngles(
        case["mass_1"] * lal.MSUN_SI, case["mass_2"] * lal.MSUN_SI,
        *case["spin_1"], *case["spin_2"], case["f_min"], case["f_max"], case["delta_f"],
        case["f_ref"], case["phi_ref"], lal_dict(dict(PhenomXPrecVersion=320)),
    )
    seconds = time.perf_counter() - start
    stride = case["stride"]
    case = dict(case)
    case["length"] = alpha.length
    case["alpha"] = alpha.data[::stride].tolist()
    case["cos_beta"] = cos_beta.data[::stride].tolist()
    case["gamma"] = gamma.data[::stride].tolist()
    return case, seconds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=1024,
                        help="approximate number of stored indices per case")
    parser.add_argument("--output", default=str(OUTPUT))
    parser.add_argument("--only", nargs="*", help="names of the binaries to generate")
    args = parser.parse_args()
    if os.environ.get("OMP_NUM_THREADS") != "1":
        raise SystemExit("run with OMP_NUM_THREADS=1")
    angles = list()
    for case in ANGLES:
        if args.only and case["name"] not in args.only:
            continue
        result, seconds = generate_angles(case)
        print(f"{case['name']:<22} SpinTaylor angles {result['length']:>8} frequencies "
              f"{len(result['alpha']):>5} stored, {seconds:.2f} s", flush=True)
        angles.append(result)
    out = list()
    for case in cases():
        if args.only and case["name"] not in args.only:
            continue
        result, seconds = generate(case, args.points)
        print(f"{case['name']:<22} {case['approximant']:<14} "
              f"{str(case['lal_flags'].get('PhenomXPrecVersion', '')):<4} "
              f"{case['generator']:<9} {result['length']:>8} frequencies "
              f"{len(result['indices']):>5} stored, {seconds:.2f} s", flush=True)
        out.append(result)
    Path(args.output).write_text(json.dumps(
        dict(lalsimulation_version=lalsim.__version__, cases=out, angles=angles)))
    print(f"wrote {args.output} ({Path(args.output).stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
