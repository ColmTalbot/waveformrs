#!/usr/bin/env python
"""
Generate reference frequency-domain waveforms from LALSimulation.

The output JSON files in ``tests/data`` are consumed by the Rust integration
tests (``tests/lalsimulation.rs``).  Each file contains a list of cases, each
case holding the input parameters, a frequency array and the plus/cross
polarizations computed with ``SimInspiralChooseFDWaveformSequence`` (or with
``SimInspiralChooseFDWaveform`` on a uniform grid starting at zero frequency
for cases with ``delta_f > 0``).  The ``imrphenomxp_spintaylor_angles`` suite
instead stores the Euler angles of ``SimIMRPhenomXPSpinTaylorAngles`` and the
``imrphenomxhm_modes`` suite the individual modes of
``SimIMRPhenomXHMFrequencySequenceOneMode``.

Usage
-----
    OMP_NUM_THREADS=1 python scripts/generate_reference_data.py [--n-frequencies 256]

``OMP_NUM_THREADS=1`` is required for the MSA suites (``imrphenomxp_msa``,
``imrphenomxp_spintaylor_fallback``): the MSA angles of lalsimulation write
per-frequency roots into the shared precession struct inside an OpenMP
loop, so multi-threaded runs are not reproducible.

Requires ``python-lalsimulation`` (tested with lalsimulation 6.2.0).
"""
import argparse
import json
from pathlib import Path

import lal
import lalsimulation as lalsim
import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent / "tests" / "data"


def lal_dict(flags, mode_array=None):
    params = lal.CreateDict()
    for key, value in (flags or dict()).items():
        getattr(lalsim, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    if mode_array is not None:
        modes = lalsim.SimInspiralCreateModeArray()
        for ell, emm in mode_array:
            lalsim.SimInspiralModeArrayActivateMode(modes, ell, emm)
        lalsim.SimInspiralWaveformParamsInsertModeArray(params, modes)
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
        lal_dict(case.get("lal_flags"), case.get("mode_array")),
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


def generate_modes(case, n_frequencies):
    """Individual IMRPhenomXHM modes from SimIMRPhenomXHMFrequencySequenceOneMode."""
    frequencies = np.geomspace(case["f_min"], case["f_max"], n_frequencies)
    if case["f_ref"] > 0:
        frequencies = np.unique(np.concatenate([frequencies, [case["f_ref"]]]))
    freq_vector = lal.CreateREAL8Vector(len(frequencies))
    freq_vector.data = frequencies
    case = dict(case)
    case["frequencies"] = frequencies.tolist()
    case["modes"] = list()
    for ell, emm in case.pop("mode_list"):
        hlm = lalsim.SimIMRPhenomXHMFrequencySequenceOneMode(
            freq_vector,
            case["mass_1"] * lal.MSUN_SI,
            case["mass_2"] * lal.MSUN_SI,
            case["chi_1"],
            case["chi_2"],
            ell,
            emm,
            case["luminosity_distance"] * 1e6 * lal.PC_SI,
            case["phi_ref"],
            case["f_ref"],
            lal_dict(case.get("lal_flags")),
        )
        case["modes"].append(
            dict(
                ell=ell,
                emm=emm,
                real=hlm.data.data.real.tolist(),
                imag=hlm.data.data.imag.tolist(),
            )
        )
    return case


def generate(case, n_frequencies):
    if "stride" in case:
        return generate_spintaylor_angles(case)
    if "mode_list" in case:
        return generate_modes(case, n_frequencies // 2)
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
        lal_dict(case.get("lal_flags"), case.get("mode_array")),
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


def uniform_precessing_cases(approximant="IMRPhenomXP", extra_flags=None, extra=()):
    """Precessing models on uniform grids (SimInspiralChooseFDWaveform)."""
    configs = list(extra) + [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, delta_f, f_max, PrecVersion
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 20.0, 1.0, 0.0, 320),
        (60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0, 0.5, 512.0, 320),
        (12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 20.0, 0.0, 1.0, 1024.0, 320),
        (12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 20.0, 25.0, 1.0, 1024.0, 102),
    ]
    return [
        dict(
            approximant=approximant,
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
            lal_flags=dict(PhenomXPrecVersion=version, **(extra_flags or dict())),
        )
        for m1, m2, s1, s2, fmin, fref, df, fmax, version in configs
    ]


def xphm_mode_array_cases():
    """IMRPhenomXPHM with restricted mode arrays (which also change the
    frequency range of the SpinTaylor angles) and equal masses and aligned
    spin components (vanishing odd co-precessing modes)."""
    base = dict(
        approximant="IMRPhenomXPHM",
        mass_1=40.0,
        mass_2=12.0,
        spin_1=[0.5, -0.2, 0.3],
        spin_2=[0.1, 0.6, -0.2],
        luminosity_distance=700.0,
        inclination=2.0,
        phi_ref=1.3,
        f_min=15.0,
        f_ref=25.0,
        f_max=0.3 / (52.0 * lal.MTSUN_SI),
    )
    flags = lambda version: dict(PhenomXPrecVersion=version, **NO_MULTIBANDING)
    return [
        dict(base, lal_flags=flags(320), mode_array=[[2, 2], [2, -2], [3, 3], [3, -3]]),
        dict(base, lal_flags=flags(320), mode_array=[[2, 1], [2, -1], [4, 4], [4, -4]]),
        dict(base, lal_flags=flags(102), mode_array=[[3, 2], [3, -2], [2, 1], [2, -1]]),
        dict(base, lal_flags=flags(320), mode_array=[[2, 2], [2, -2]]),
        dict(
            base,
            mass_1=30.0,
            mass_2=30.0,
            spin_1=[0.3, 0.0, 0.2],
            spin_2=[0.0, -0.4, 0.2],
            f_max=0.3 / (60.0 * lal.MTSUN_SI),
            lal_flags=flags(320),
        ),
        dict(
            base,
            mass_1=30.0,
            mass_2=30.0,
            spin_1=[0.3, 0.0, 0.2],
            spin_2=[0.0, -0.4, 0.2],
            f_max=0.3 / (60.0 * lal.MTSUN_SI),
            lal_flags=flags(102),
        ),
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


def msa_cases():
    """IMRPhenomXP with the MSA angles (PrecVersion 223), on non-uniform and
    uniform grids, including a binary whose MSA system fails to initialise
    (|Omega_z5| > 1000), for which lalsimulation uses the NNLO angles."""
    cases = precessing_cases(
        "IMRPhenomXP", f_max_imr(0.3), lal_flags=dict(PhenomXPrecVersion=223)
    )
    msa_error = dict(
        cases[0],
        mass_1=5.223612204949848,
        mass_2=3.7902144714569808,
        spin_1=[0.28093860818235766, 0.48620119720341687, -0.09813236836152804],
        spin_2=[-0.02328864377636709, -0.684967014257048, -0.6656211288407142],
        f_min=20.0,
        f_ref=59.33970557086933,
    )
    msa_error["f_max"] = f_max_imr(0.3)(msa_error["mass_1"] + msa_error["mass_2"])
    uniform = [
        dict(case, lal_flags=dict(PhenomXPrecVersion=223))
        for case in uniform_precessing_cases()[1:3]
    ]
    return cases + [msa_error] + uniform


def spintaylor_fallback_cases():
    """IMRPhenomXP with the SpinTaylor angles (PrecVersion 320) for heavy
    binaries where the PN integration does not reach f_min, so that
    lalsimulation falls back to the MSA angles (223)."""
    configs = [
        # mass_1, mass_2, spin_1, spin_2, inclination, phi_ref, f_min, f_ref, delta_f
        (180.0, 65.0, [0.2, 0.1, -0.1], [-0.1, -0.2, 0.1], 0.4, 0.3, 15.5, 15.5, 0.0),
        (150.0, 100.0, [0.5, -0.3, 0.2], [-0.4, 0.4, 0.3], 1.3, 2.1, 20.0, 20.0, 0.0),
        (165.0, 95.0, [-0.6, 0.2, 0.5], [0.3, 0.6, -0.2], 2.2, 4.4, 20.0, 22.0, 0.0),
        (250.0, 30.0, [0.7, 0.3, -0.4], [0.0, -0.5, 0.5], 0.9, 5.9, 10.0, 25.0, 0.0),
        (140.0, 135.0, [0.1, 0.8, 0.3], [0.6, -0.1, -0.6], 1.7, 0.8, 15.0, 15.0, 0.0),
        (100.0, 60.0, [-0.3, -0.7, -0.5], [0.2, 0.2, 0.7], 2.8, 3.6, 20.0, 50.0, 0.0),
        (200.0, 70.0, [0.4, -0.6, 0.1], [-0.5, 0.0, -0.4], 0.6, 1.4, 20.0, 20.0, 0.25),
        (130.0, 90.0, [0.0, 0.5, -0.6], [0.4, 0.4, 0.2], 1.1, 2.7, 15.0, 30.0, 0.125),
    ]
    cases = list()
    for m1, m2, s1, s2, iota, phi, fmin, fref, df in configs:
        case = dict(
            approximant="IMRPhenomXP",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=1500.0,
            inclination=iota,
            phi_ref=phi,
            f_min=fmin,
            f_ref=fref,
            f_max=f_max_imr(0.3)(m1 + m2),
            lal_flags=dict(PhenomXPrecVersion=320),
        )
        if df > 0:
            case.update(delta_f=df, f_max=0.0)
        cases.append(case)
    return cases


def xphm_msa_cases():
    """IMRPhenomXPHM with the MSA angles (PrecVersion 223), multibanding off,
    on non-uniform and uniform grids, including the binary whose MSA system
    fails to initialise (lalsimulation then uses the NNLO angles)."""
    flags = dict(PhenomXPrecVersion=223, **NO_MULTIBANDING)
    cases = [
        dict(case, approximant="IMRPhenomXPHM", lal_flags=flags)
        for case in msa_cases()
    ]
    return cases


def xphm_spintaylor_fallback_cases():
    """IMRPhenomXPHM with the SpinTaylor angles (PrecVersion 320, multibanding
    off) for heavy binaries where lalsimulation falls back to the MSA angles."""
    configs = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, delta_f
        (163.5, 132.1, [0.32, -0.49, -0.14], [-0.46, -0.34, -0.31], 20.0, 25.1, 0.0),
        (128.3, 102.2, [-0.12, -0.53, 0.03], [-0.32, 0.27, -0.12], 20.0, 35.7, 0.0),
        (124.8, 87.5, [0.42, 0.55, -0.39], [0.04, 0.42, -0.49], 20.0, 27.0, 0.0),
        (136.4, 88.5, [0.07, -0.34, 0.0], [0.47, 0.45, -0.53], 20.0, 33.1, 0.0),
        (143.8, 43.2, [0.17, -0.5, -0.28], [-0.22, -0.27, 0.25], 20.0, 56.2, 0.0),
        (186.9, 68.2, [-0.19, -0.41, -0.28], [0.12, 0.25, -0.22], 15.0, 43.8, 0.0),
        (139.8, 22.9, [0.29, 0.0, -0.08], [0.15, 0.29, -0.47], 15.0, 43.7, 1.0),
        (81.2, 77.2, [0.31, -0.14, -0.34], [0.4, 0.32, 0.19], 20.0, 59.1, 0.5),
    ]
    cases = list()
    for m1, m2, s1, s2, fmin, fref, df in configs:
        case = dict(
            approximant="IMRPhenomXPHM",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=1000.0,
            inclination=1.0,
            phi_ref=0.5,
            f_min=fmin,
            f_ref=fref,
            f_max=f_max_imr(0.3)(m1 + m2),
            lal_flags=dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
        )
        if df > 0:
            case.update(delta_f=df, f_max=0.0)
        cases.append(case)
    return cases


def f_max_isco(total_mass):
    """Schwarzschild ISCO GW frequency, TaylorF2's default cutoff."""
    return 1.0 / (6.0**1.5 * np.pi * total_mass * lal.MTSUN_SI)


def f_max_imr(mf_max):
    def func(total_mass):
        return min(mf_max / (total_mass * lal.MTSUN_SI), 4096.0)

    return func


NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0.0, PhenomXPHMThresholdMband=0.0)


def xhm_mode_cases():
    """Each mode of IMRPhenomXHM separately, including equal masses and
    spins (vanishing odd modes), extreme mass ratios and high spins."""
    negative = [(2, -2), (2, -1), (3, -3), (3, -2), (4, -4)]
    configs = [
        # mass_1, mass_2, chi_1, chi_2, phi_ref, f_min, f_ref, modes
        (36.0, 29.0, 0.3, -0.4, 1.2, 20.0, 30.0, negative + [(2, 1), (3, 2)]),
        (30.0, 30.0, 0.5, 0.5, 0.4, 20.0, 20.0, negative),
        (60.0, 15.0, 0.8, 0.5, 4.0, 15.0, 20.0, negative),
        (80.0, 10.0, -0.7, 0.2, 0.3, 10.0, 25.0, negative),
        (150.0, 5.0, 0.9, -0.3, 2.5, 5.0, 8.0, negative),
        (240.0, 2.0, -0.6, 0.4, 5.0, 3.0, 4.0, negative),
        (8.0, 1.4, 0.9, 0.0, 5.5, 20.0, 50.0, negative + [(4, 4)]),
        (150.0, 100.0, 0.95, 0.95, 3.0, 10.0, 15.0, negative),
    ]
    cases = list()
    for m1, m2, c1, c2, phi, fmin, fref, modes in configs:
        cases.append(
            dict(
                approximant="IMRPhenomXHM",
                mass_1=m1,
                mass_2=m2,
                chi_1=c1,
                chi_2=c2,
                luminosity_distance=500.0,
                phi_ref=phi,
                f_min=fmin,
                f_ref=fref,
                f_max=f_max_imr(0.3)(m1 + m2),
                mode_list=modes,
            )
        )
    return cases


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
    imrphenomxhm_modes=xhm_mode_cases,
    imrphenomxp_nnlo=lambda: precessing_cases(
        "IMRPhenomXP", f_max_imr(0.3), lal_flags=dict(PhenomXPrecVersion=102)
    ),
    imrphenomxp_spintaylor=lambda: precessing_cases(
        "IMRPhenomXP", f_max_imr(0.3), lal_flags=dict(PhenomXPrecVersion=320)
    ),
    imrphenomxp_uniform=uniform_precessing_cases,
    imrphenomxp_msa=msa_cases,
    imrphenomxp_spintaylor_fallback=spintaylor_fallback_cases,
    imrphenomxp_spintaylor_angles=spintaylor_angle_cases,
    imrphenomxphm_spintaylor=lambda: precessing_cases(
        "IMRPhenomXPHM",
        f_max_imr(0.3),
        lal_flags=dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
    ),
    imrphenomxphm_nnlo=lambda: precessing_cases(
        "IMRPhenomXPHM",
        f_max_imr(0.3),
        lal_flags=dict(PhenomXPrecVersion=102, **NO_MULTIBANDING),
    ),
    imrphenomxphm_uniform=lambda: uniform_precessing_cases(
        "IMRPhenomXPHM",
        NO_MULTIBANDING,
        # equal masses: the uniform driver keeps the input order of the bodies
        extra=[(30.0, 30.0, [0.5, 0.5, 0.5], [0.5, -0.5, -0.5], 20.0, 30.0, 1.0, 512.0, 320)],
    ),
    imrphenomxphm_modes=xphm_mode_array_cases,
    imrphenomxphm_msa=xphm_msa_cases,
    imrphenomxphm_spintaylor_fallback=xphm_spintaylor_fallback_cases,
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
