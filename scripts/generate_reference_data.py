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
``SimIMRPhenomXHMFrequencySequenceOneMode`` (``imrphenomxhm_122019_modes``
and ``imrphenomxhm_122019`` for the original release of the higher-mode
amplitudes, ``PhenomXHMReleaseVersion = 122019``). The ``choose_fd_modes``
suite stores the modes of ``SimInspiralChooseFDModes`` (IMRPhenomXHM and
IMRPhenomXPHM, on uniform grids), every ``fd_modes_stride``-th sample of the
positive- and negative-frequency halves of each mode.

The ``imrphenomxpnr*`` suites cover IMRPhenomXPNR: polarizations on
frequency sequences (``imrphenomxpnr``, including the binaries for which
LAL switches parts of the PNR machinery off or changes the Euler angles)
and uniform grids (``imrphenomxpnr_uniform``), the J-frame modes of
``SimInspiralChooseFDModes`` (``imrphenomxpnr_modes``) and of
``SimIMRPhenomXPHMFrequencySequenceOneMode`` with the flags of
IMRPhenomXPNR (``imrphenomxpnr_sequence_modes``).

The time-domain ``nrsur7dq4`` suites store the polarizations of
``SimInspiralChooseTDWaveform`` (with the epoch, or the error message for
the cases where lalsimulation fails), the modes of
``SimInspiralChooseTDModes`` (``nrsur7dq4_modes``) and the output of
``PrecessingNRSurDynamics`` (``nrsur7dq4_dynamics``). They need the data
file ``NRSur7dq4_v1.0.h5`` in ``LAL_DATA_PATH``.

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


def nrsur_start_frequency(case):
    """The lowest f_min of NRSur7dq4 for a binary
    (PrecessingNRSur_StartFrequency), by bisection on f_ref with f_min = 0:
    lalsimulation fails for reference frequencies below it."""

    def valid(f_ref):
        try:
            lalsim.PrecessingNRSurDynamics(*nrsur_dynamics_args(case, f_ref))
            return True
        except RuntimeError:
            return False

    lo, hi = 1.0, 200.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if valid(mid) else (mid, hi)
    return hi


def nrsur_dynamics_args(case, f_ref):
    """Arguments of PrecessingNRSurDynamics for a waveform case (heavier body
    first) with LAL's waveform conventions at f_ref."""
    m1 = case["mass_1"] * lal.MSUN_SI
    m2 = case["mass_2"] * lal.MSUN_SI
    s1, s2 = case["spin_1"], case["spin_2"]
    if m1 < m2:
        m1, m2 = m2, m1
        s1, s2 = [-s2[0], -s2[1], s2[2]], [-s1[0], -s1[1], s1[2]]
    mass_a, mass_b = m1 / lal.MSUN_SI, m2 / lal.MSUN_SI
    m_tot_sec = (mass_a + mass_b) * lal.MTSUN_SI
    omega_ref = (m_tot_sec * f_ref) * np.pi
    return (mass_a / mass_b, *s1, *s2, omega_ref, 1.0, 0.0, 0.0, 0.0, 0.0, None,
            lalsim.NRSur7dq4)


def nrsur_case_arguments(case):
    if case.get("f_ref") == "start":
        # just above the surrogate's first node: the RK4 steps start from a
        # node other than 0 without the half nodes
        case["f_ref"] = nrsur_start_frequency(case) * case.pop("f_ref_factor")
    return (
        case["mass_1"] * lal.MSUN_SI,
        case["mass_2"] * lal.MSUN_SI,
        *case["spin_1"],
        *case["spin_2"],
    )


def nrsur_dict(case):
    """The LALDict of a NRSur7dq4 case: mode array and extrapolation."""
    params = lal_dict(case.get("lal_flags"), case.get("mode_array"))
    if case.get("unlimited_extrapolation"):
        lal.DictInsertUINT4Value(params, "unlimited_extrapolation", 1)
    return params


def generate_td(case):
    """Polarizations of SimInspiralChooseTDWaveform (or the error)."""
    case = dict(case)
    masses_spins = nrsur_case_arguments(case)
    try:
        hplus, hcross = lalsim.SimInspiralChooseTDWaveform(
            *masses_spins,
            case["luminosity_distance"] * 1e6 * lal.PC_SI,
            case["inclination"],
            case["phi_ref"],
            0.0,
            0.0,
            0.0,
            case["delta_t"],
            case["f_min"],
            case["f_ref"],
            nrsur_dict(case),
            lalsim.GetApproximantFromString(case["approximant"]),
        )
    except RuntimeError as error:
        case["lal_error"] = str(error)
        return case
    case["epoch"] = float(hplus.epoch)
    case["hplus"] = hplus.data.data.tolist()
    case["hcross"] = hcross.data.data.tolist()
    return case


def generate_td_modes(case):
    """Modes of SimInspiralChooseTDModes, every stride-th sample."""
    case = dict(case)
    masses_spins = nrsur_case_arguments(case)
    hlms = lalsim.SimInspiralChooseTDModes(
        0.0,
        case["delta_t"],
        *masses_spins,
        case["f_min"],
        case["f_ref"],
        case["luminosity_distance"] * 1e6 * lal.PC_SI,
        nrsur_dict(case),
        4,
        lalsim.GetApproximantFromString(case["approximant"]),
    )
    stride = case["sample_stride"]
    case["modes"] = list()
    for ell in range(2, 5):
        for emm in range(-ell, ell + 1):
            hlm = lalsim.SphHarmTimeSeriesGetMode(hlms, ell, emm)
            if hlm is None:
                # all modes of this ell are inactive in the mode array
                continue
            case["epoch"] = float(hlm.epoch)
            case["length"] = hlm.data.length
            case["modes"].append(
                dict(
                    ell=ell,
                    emm=emm,
                    real=hlm.data.data.real[::stride].tolist(),
                    imag=hlm.data.data.imag[::stride].tolist(),
                )
            )
    return case


def generate_dynamics(case):
    """XLALPrecessingNRSurDynamics on the dynamics nodes."""
    case = dict(case)
    out = lalsim.PrecessingNRSurDynamics(
        case["q"],
        *case["chi_a"],
        *case["chi_b"],
        case["omega_ref"],
        *case["init_quat"],
        case["init_orbphase"],
        None,
        lalsim.GetApproximantFromString(case["approximant"]),
    )
    names = ["t", "q0", "qx", "qy", "qz", "orbphase", "chiAx", "chiAy", "chiAz",
             "chiBx", "chiBy", "chiBz"]
    for name, values in zip(names, out):
        case[name] = values.data.tolist()
    return case


def sincos_sensitive_angles(n_per_kind, rng):
    """Polar angles in [0, pi] for which Darwin's __sincos_stret returns a
    different sine than sin() for theta/2, theta, 2 theta, 3 theta or
    4 theta (empty elsewhere): clang combines sin and cos of the same
    argument into __sincos_stret, so these angles test which calls the
    harmonics make."""
    import ctypes
    import platform

    if platform.system() != "Darwin":
        return []

    class SinCos(ctypes.Structure):
        _fields_ = [("sin", ctypes.c_double), ("cos", ctypes.c_double)]

    libm = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
    libm.__sincos_stret.restype = SinCos
    libm.__sincos_stret.argtypes = [ctypes.c_double]
    libm.sin.restype = ctypes.c_double
    libm.sin.argtypes = [ctypes.c_double]
    found = [0] * 5
    angles = list()
    while min(found) < n_per_kind:
        theta = float(rng.uniform(0, np.pi))
        for k, arg in enumerate([theta / 2, theta, 2 * theta, 3 * theta, 4 * theta]):
            if found[k] < n_per_kind and libm.__sincos_stret(arg).sin != libm.sin(arg):
                found[k] += 1
                angles.append(theta)
    return angles


def spin_weighted_harmonic_cases():
    """XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m) for l <= 4."""
    rng = np.random.default_rng(1729)
    thetas = sincos_sensitive_angles(20, rng)
    thetas += [0.0, np.pi / 2, np.pi, 1.5321585045470605, 1.5168235141421997, 1.497194552179717]
    thetas += rng.uniform(0, np.pi, 40).tolist()
    phis = rng.uniform(-np.pi, 2 * np.pi, len(thetas)).tolist()
    return [dict(theta=thetas, phi=phis, harmonics=True)]


def generate_fd_modes(case):
    """Modes of SimInspiralChooseFDModes. Each mode is a two-sided series
    on (-L..L) * delta_f; store h_lm(k delta_f) and h_lm(-k delta_f) for
    k = 0, stride, 2 stride, ... <= L."""
    params = lal_dict(case.get("lal_flags"), case.get("mode_array"))
    if case.get("mode_array_jframe") is not None:
        modes = lalsim.SimInspiralCreateModeArray()
        for ell, emm in case["mode_array_jframe"]:
            lalsim.SimInspiralModeArrayActivateMode(modes, ell, emm)
        lalsim.SimInspiralWaveformParamsInsertModeArrayJframe(params, modes)
    hlms = lalsim.SimInspiralChooseFDModes(
        case["mass_1"] * lal.MSUN_SI,
        case["mass_2"] * lal.MSUN_SI,
        *case["spin_1"],
        *case["spin_2"],
        case["delta_f"],
        case["f_min"],
        case["f_max"],
        case["f_ref"],
        case["phi_ref"],
        case["luminosity_distance"] * 1e6 * lal.PC_SI,
        case["inclination"],
        params,
        lalsim.GetApproximantFromString(case["approximant"]),
    )
    stride = case["fd_modes_stride"]
    case = dict(case)
    case["modes"] = list()
    node = hlms
    while node is not None:
        data = node.mode.data.data
        half = (len(data) - 1) // 2
        positive = data[half::stride]
        negative = data[half::-stride]
        case["length"] = half + 1
        case["modes"].append(
            dict(
                ell=node.l,
                emm=node.m,
                positive_real=positive.real.tolist(),
                positive_imag=positive.imag.tolist(),
                negative_real=negative.real.tolist(),
                negative_imag=negative.imag.tolist(),
            )
        )
        node = node.next
    case["modes"].sort(key=lambda mode: (mode["ell"], mode["emm"]))
    return case


def generate_sequence_modes(case, n_frequencies):
    """J-frame modes of SimIMRPhenomXPHMFrequencySequenceOneMode (with the
    flags of the case) for every l of the co-precessing modes, |m| <= l."""
    frequencies = np.geomspace(case["f_min"], case["f_max"], n_frequencies)
    freq_vector = lal.CreateREAL8Vector(len(frequencies))
    freq_vector.data = frequencies
    case = dict(case)
    case.pop("jframe_sequence")
    case["frequencies"] = frequencies.tolist()
    case["modes"] = list()
    ells = sorted({ell for ell, _ in case.get("mode_array") or [(2, 2), (3, 3), (4, 4)]})
    for ell in ells:
        for emm in range(-ell, ell + 1):
            positive, negative = lalsim.SimIMRPhenomXPHMFrequencySequenceOneMode(
                freq_vector,
                ell,
                emm,
                case["mass_1"] * lal.MSUN_SI,
                case["mass_2"] * lal.MSUN_SI,
                *case["spin_1"],
                *case["spin_2"],
                case["luminosity_distance"] * 1e6 * lal.PC_SI,
                case["inclination"],
                case["phi_ref"],
                case["f_ref"],
                lal_dict(case.get("lal_flags"), case.get("mode_array")),
            )
            case["modes"].append(
                dict(
                    ell=ell,
                    emm=emm,
                    positive_real=positive.data.data.real.tolist(),
                    positive_imag=positive.data.data.imag.tolist(),
                    negative_real=negative.data.data.real.tolist(),
                    negative_imag=negative.data.data.imag.tolist(),
                )
            )
    return case


def generate_harmonics(case):
    case = dict(case)
    case.pop("harmonics")
    case["values"] = list()
    for theta, phi in zip(case["theta"], case["phi"]):
        for ell in range(2, 5):
            for emm in range(-ell, ell + 1):
                y = lal.SpinWeightedSphericalHarmonic(theta, phi, -2, ell, emm)
                case["values"].append([ell, emm, y.real, y.imag])
    return case


def generate(case, n_frequencies):
    if "harmonics" in case:
        return generate_harmonics(case)
    if "fd_modes_stride" in case:
        return generate_fd_modes(case)
    if "jframe_sequence" in case:
        return generate_sequence_modes(case, n_frequencies // 2)
    if "sample_stride" in case:
        return generate_td_modes(case)
    if "delta_t" in case:
        return generate_td(case)
    if "omega_ref" in case:
        return generate_dynamics(case)
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


XHM_RELEASE_122019 = dict(PhenomXHMReleaseVersion=122019)

# Binaries whose higher modes exercise the branches of the 122019 amplitude:
# inspiral vetoes (q < 8 for the (2,1)) and wavy collocation points, the
# TaylorT4 (2,1) inspiral (q < 40) and its power series (q > 40), the
# (3,2) and (3,3) corners of the inspiral, the intermediate vetoes and the
# lowering of the polynomial order when it crosses zero (5th to 2nd order),
# two intermediate regions (q > 70, chi1 <= 0.9), the inspiral ansatz in the
# ringdown of the (3,3), (4,4) (q > 7, chi1 > 0.95), the (3,2) ringdown
# cutoffs, negative amplitudes and vanishing odd modes.
XHM_122019_BINARIES = [
    # mass_1, mass_2, chi_1, chi_2, f_ref
    (36.0, 29.0, 0.3, -0.4, 30.0),
    (30.0, 30.0, 0.5, -0.5, 20.0),
    (150.2452, 2.957, 0.945, -1.0, 30.0),
    (79.0458, 0.2689, -0.89, 0.959, 0.0),
    (10.7049, 10.1951, 0.63, 0.61, 0.0),
    (100.9494, 12.6503, 0.988, -0.02, 30.0),
    (148.5018, 3.4991, -0.914, -0.93, 20.0),
    (109.3352, 91.2648, -0.966, 0.991, 0.0),
    (66.0796, 65.6203, 0.77, 0.99, 0.0),
    (109.45, 109.45, -0.947, -0.947, 20.0),
    (91.0032, 63.1967, 0.989, 0.31, 30.0),
    (118.8765, 113.3236, -0.991, 0.63, 30.0),
    (28.1353, 9.1646, 0.907, -0.97, 30.0),
    (200.0, 2.0, 0.95, 0.3, 15.0),
]


def xhm_122019_mode_cases():
    """The higher modes of IMRPhenomXHM with PhenomXHMReleaseVersion = 122019."""
    negative = [(2, -1), (3, -3), (3, -2), (4, -4)]
    cases = list()
    for i, (m1, m2, c1, c2, fref) in enumerate(XHM_122019_BINARIES):
        modes = negative + [(2, 1), (3, 2)] if i == 0 else negative
        cases.append(
            dict(
                approximant="IMRPhenomXHM",
                mass_1=m1,
                mass_2=m2,
                chi_1=c1,
                chi_2=c2,
                luminosity_distance=500.0,
                phi_ref=0.1 + 0.4 * i,
                f_min=10.0,
                f_ref=fref,
                f_max=f_max_imr(0.3)(m1 + m2),
                mode_list=modes,
                lal_flags=XHM_RELEASE_122019,
            )
        )
    return cases


def xhm_122019_cases():
    """Polarizations of IMRPhenomXHM with PhenomXHMReleaseVersion = 122019
    (the (3,2) recycles the (2,2) in its ringdown)."""
    cases = list()
    for i in (0, 3, 5, 6, 8, 9, 11, 13):
        m1, m2, c1, c2, fref = XHM_122019_BINARIES[i]
        cases.append(
            dict(
                approximant="IMRPhenomXHM",
                mass_1=m1,
                mass_2=m2,
                chi_1=c1,
                chi_2=c2,
                luminosity_distance=500.0,
                inclination=0.3 + 0.35 * i,
                phi_ref=0.2 + 0.5 * i,
                f_min=10.0,
                f_ref=fref,
                f_max=f_max_imr(0.3)(m1 + m2),
                lal_flags=dict(XHM_RELEASE_122019, **NO_MULTIBANDING),
            )
        )
    return cases


def unit(vector, magnitude):
    vector = np.asarray(vector, dtype=float)
    return (magnitude * vector / np.linalg.norm(vector)).tolist()


def choose_fd_modes_cases():
    """SimInspiralChooseFDModes for IMRPhenomXHM (mode arrays, equal masses,
    f_max above the cutoff, f_ref = 0) and IMRPhenomXPHM (SpinTaylor, MSA
    and NNLO angles, a co-precessing mode array with the matching J-frame
    mode array, equal masses, and a heavy binary whose SpinTaylor angles
    fall back to the MSA angles for some J-frame modes), with
    multibanding off."""
    xhm = [
        # mass_1, mass_2, chi_1, chi_2, f_min, f_ref, delta_f, f_max, mode_array
        (36.0, 29.0, 0.3, -0.4, 20.0, 30.0, 0.5, 512.0, None),
        (60.0, 15.0, 0.8, 0.5, 15.0, 0.0, 1.0, 0.0, None),
        (30.0, 30.0, 0.5, 0.5, 20.0, 25.0, 1.0, 1024.0, [[2, 2], [2, -2], [3, -3], [4, 4]]),
        (8.0, 1.4, 0.9, 0.0, 20.0, 50.0, 0.25, 256.0, [[2, -2], [2, -1], [3, -3], [3, -2], [4, -4]]),
    ]
    cases = [
        dict(
            approximant="IMRPhenomXHM",
            mass_1=m1,
            mass_2=m2,
            spin_1=[0.0, 0.0, c1],
            spin_2=[0.0, 0.0, c2],
            luminosity_distance=500.0,
            inclination=0.7,
            phi_ref=1.2,
            f_min=fmin,
            f_ref=fref,
            delta_f=df,
            f_max=fmax,
            lal_flags=NO_MULTIBANDING,
            fd_modes_stride=8,
            **(dict(mode_array=modes) if modes else dict()),
        )
        for m1, m2, c1, c2, fmin, fref, df, fmax, modes in xhm
    ]
    l23 = [[ell, emm] for ell in (2, 3) for emm in range(-ell, ell + 1)]
    xphm = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, delta_f, f_max,
        # version, inclination, extra
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 20.0, 1.0, 512.0, 320, 0.7, {}),
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 25.0, 1.0, 512.0, 223, 0.7, {}),
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 0.0, 1.0, 512.0, 102, 0.7, {}),
        (60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0, 0.5, 0.0, 320, 2.3, {}),
        (40.0, 12.0, [0.5, -0.2, 0.3], [0.1, 0.6, -0.2], 15.0, 25.0, 1.0, 0.0, 320, 2.0,
         dict(mode_array=[[2, 2], [2, -2], [2, 1], [2, -1], [3, 3], [3, -3]],
              mode_array_jframe=l23)),
        (30.0, 30.0, [0.5, 0.5, 0.5], [0.5, -0.5, -0.5], 20.0, 30.0, 1.0, 512.0, 320, 1.6, {}),
        (180.0, 65.0, [0.2, 0.1, -0.1], [-0.1, -0.2, 0.1], 15.5, 15.5, 0.5, 0.0, 320, 0.4, {}),
    ]
    for m1, m2, s1, s2, fmin, fref, df, fmax, version, iota, extra in xphm:
        cases.append(
            dict(
                approximant="IMRPhenomXPHM",
                mass_1=m1,
                mass_2=m2,
                spin_1=s1,
                spin_2=s2,
                luminosity_distance=800.0,
                inclination=iota,
                phi_ref=0.9,
                f_min=fmin,
                f_ref=fref,
                delta_f=df,
                f_max=fmax,
                lal_flags=dict(PhenomXPrecVersion=version, **NO_MULTIBANDING),
                fd_modes_stride=8,
                **extra,
            )
        )
    return cases


# The flags SimInspiralChooseFDWaveform(Sequence) sets for IMRPhenomXPNR,
# for the IMRPhenomXPHM functions called directly.
XPNR_FLAGS = dict(
    PhenomXHMReleaseVersion=122019,
    PhenomXPNRUseTunedAngles=1,
    PhenomXPNRUseTunedCoprec=1,
    PhenomXPNRUseTunedCoprec33=0,
    PhenomXPNRForceXHMAlignment=0,
    PhenomXAntisymmetricWaveform=1,
    PhenomXPrecVersion=330,
    PhenomXPFinalSpinMod=7,
    PhenomXPHMMBandVersion=2,
)

# mass_1, mass_2, spin_1, spin_2, f_min, f_ref: binaries for which LAL's
# PNR machinery takes its less common branches
XPNR_EDGE_BINARIES = [
    # in-plane spins below 1e-7: plain IMRPhenomXPHM (NNLO angles, 122019
    # release, FinalSpinMod 7), with a negative final spin
    (30.0, 20.0, [0.0, 0.0, -0.8], [0.0, 0.0, -0.7], 20.0, 20.0),
    (30.0, 20.0, [5e-8, 0.0, 0.5], [0.0, 0.0, -0.3], 20.0, 25.0),
    # in-plane spins between 1e-7 and 1e-6: NNLO angles inside PNR
    (30.0, 20.0, [5e-7, 0.0, 0.5], [0.0, 0.0, -0.3], 20.0, 20.0),
    # just above 1e-6: SpinTaylor angles with a tiny in-plane spin
    (30.0, 20.0, [2e-6, 0.0, 0.5], [0.0, 0.0, -0.3], 20.0, 20.0),
    # transition region in q and in chi
    (90.0, 10.0, [0.5, 0.2, 0.1], [0.1, 0.0, 0.2], 15.0, 15.0),
    (110.0, 10.0, [0.5, 0.2, 0.1], [0.1, 0.0, 0.2], 15.0, 20.0),
    (30.0, 10.0, [0.9, 0.1, 0.3], [0.1, 0.0, 0.2], 20.0, 20.0),
    # outside the calibration region (PN angles, inspiral scaling of the
    # frequency map, co-precessing deviations tapered for 10 < q <= 20)
    (130.0, 10.0, [0.5, 0.2, 0.1], [0.1, 0.0, 0.2], 10.0, 10.0),
    (190.0, 10.0, [0.7, 0.2, -0.1], [0.1, 0.4, 0.2], 10.0, 12.0),
    (30.0, 5.0, [0.98, 0.1, 0.1], [0.5, 0.0, 0.2], 20.0, 20.0),
    # equal masses and spins: the odd co-precessing modes vanish
    (30.0, 30.0, [0.3, 0.2, 0.1], [0.3, 0.2, 0.1], 20.0, 20.0),
    # the SpinTaylor integration fails: MSA angles (223) inside PNR
    (250.0, 30.0, [0.7, 0.3, -0.4], [0.0, -0.5, 0.5], 10.0, 25.0),
    (180.0, 65.0, [0.2, 0.1, -0.1], [-0.1, -0.2, 0.1], 15.5, 15.5),
    # the MSA system fails as well: NNLO angles (102) inside PNR
    (324.97464539595285, 156.1317450178247, [0.033473504803287535, 0.3586421803050725, -0.8463397636474128],
     [-0.0181461364021327, -0.7837686724310844, 0.22548210669736396], 20.0, 20.0),
    # the SpinTaylor angles are integrated again from a lower frequency
    (65.39462569822626, 5.48546881140544, [0.015777294748968645, -0.01129185185344695, 0.0015599118587599989],
     [0.26347849219372566, 0.5361809956686221, -0.26720621117427656], 30.0, 52.60747784129517),
    # no merger-ringdown beta attached
    (238.87854968975887, 13.840439770449478, [-0.006661588927623341, 0.012223592277567882, -0.02867622012464532],
     [0.21532597836655717, -0.9596190530452179, 0.06428272454610019], 30.0, 36.54489151780514),
    # the lower connection frequency of beta is not found (Mf = 100)
    (141.98456902640373, 34.571276664666726, [0.5558143552228454, -0.08580576236737993, -0.781422167551453],
     [-0.4951057605866491, 0.15954331438844666, -0.5788066659662768], 10.0, 34.6828483981249),
    # a single precessing spin
    (44.776077533636894, 14.83112049039989, [-0.05546660880373749, -0.06115453893757115, 0.36144221300916063],
     [0.0, 0.0, 0.0], 30.0, 30.0),
    # light binary, long inspiral
    (5.0, 3.0, [0.4, -0.3, 0.2], [0.2, 0.3, -0.5], 20.0, 20.0),
    # high anti-aligned spins
    (40.0, 10.0, [0.1, 0.1, -0.98], [0.0, 0.2, -0.9], 20.0, 30.0),
    # f_ref far above f_min, and the lighter body first
    (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 200.0),
    (20.0, 45.0, [0.3, 0.4, -0.2], [-0.5, 0.1, 0.6], 15.0, 22.0),
]


def xpnr_cases():
    """IMRPhenomXPNR polarizations from SimInspiralChooseFDWaveformSequence:
    generic precessing binaries and the edge cases of the PNR machinery,
    and two restricted mode arrays."""
    cases = precessing_cases("IMRPhenomXPNR", f_max_imr(0.3))
    for i, (m1, m2, s1, s2, fmin, fref) in enumerate(XPNR_EDGE_BINARIES):
        cases.append(
            dict(
                approximant="IMRPhenomXPNR",
                mass_1=m1,
                mass_2=m2,
                spin_1=s1,
                spin_2=s2,
                luminosity_distance=400.0 + 50.0 * i,
                inclination=0.2 + 0.13 * i,
                phi_ref=0.1 + 0.27 * i,
                f_min=fmin,
                f_ref=fref,
                f_max=f_max_imr(0.3)(m1 + m2),
                lal_flags=dict(),
            )
        )
    for m1, m2, modes in [
        (36.0, 29.0, [[2, 1], [2, -1], [3, 3], [3, -3]]),
        (60.0, 12.0, [[3, 2], [3, -2], [4, 4], [4, -4]]),
    ]:
        cases.append(
            dict(
                approximant="IMRPhenomXPNR",
                mass_1=m1,
                mass_2=m2,
                spin_1=[0.4, 0.2, 0.3],
                spin_2=[-0.2, 0.3, -0.1],
                luminosity_distance=500.0,
                inclination=1.1,
                phi_ref=0.4,
                f_min=20.0,
                f_ref=25.0,
                f_max=f_max_imr(0.3)(m1 + m2),
                lal_flags=dict(),
                mode_array=modes,
            )
        )
    return cases


def xpnr_uniform_cases():
    """IMRPhenomXPNR on uniform grids (SimInspiralChooseFDWaveform, which
    keeps the input order of the bodies in the waveform struct)."""
    configs = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, delta_f, f_max
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 20.0, 1.0, 0.0),
        (60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0, 0.5, 512.0),
        (12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 20.0, 0.0, 1.0, 1024.0),
        (20.0, 45.0, [0.3, 0.4, -0.2], [-0.5, 0.1, 0.6], 20.0, 25.0, 0.5, 0.0),
        (30.0, 20.0, [0.0, 0.0, 0.5], [0.0, 0.0, -0.3], 20.0, 20.0, 1.0, 0.0),
        (250.0, 30.0, [0.7, 0.3, -0.4], [0.0, -0.5, 0.5], 10.0, 25.0, 0.25, 0.0),
        (130.0, 10.0, [0.5, 0.2, 0.1], [0.1, 0.0, 0.2], 10.0, 12.0, 0.25, 0.0),
    ]
    return [
        dict(
            approximant="IMRPhenomXPNR",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=600.0,
            inclination=0.9,
            phi_ref=1.7,
            f_min=fmin,
            f_ref=fref,
            delta_f=df,
            f_max=fmax,
            lal_flags=dict(NO_MULTIBANDING),
        )
        for m1, m2, s1, s2, fmin, fref, df, fmax in configs
    ]


def xpnr_mode_cases():
    """The J-frame modes of IMRPhenomXPNR from SimInspiralChooseFDModes
    (one set-up per J-frame mode, or a co-precessing mode array with the
    matching J-frame mode array), with multibanding off."""
    l23 = [[ell, emm] for ell in (2, 3) for emm in range(-ell, ell + 1)]
    configs = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, delta_f, f_max, extra
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 25.0, 1.0, 512.0, {}),
        (60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0, 0.5, 0.0, {}),
        (40.0, 12.0, [0.5, -0.2, 0.3], [0.1, 0.6, -0.2], 15.0, 25.0, 1.0, 0.0,
         dict(mode_array=[[2, 2], [2, -2], [2, 1], [2, -1], [3, 3], [3, -3]], mode_array_jframe=l23)),
        (30.0, 30.0, [0.5, 0.5, 0.5], [0.5, -0.5, -0.5], 20.0, 30.0, 1.0, 512.0, {}),
        (30.0, 20.0, [0.0, 0.0, 0.5], [0.0, 0.0, -0.3], 20.0, 20.0, 1.0, 0.0, {}),
        (180.0, 65.0, [0.2, 0.1, -0.1], [-0.1, -0.2, 0.1], 15.5, 15.5, 0.5, 0.0, {}),
    ]
    return [
        dict(
            approximant="IMRPhenomXPNR",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=800.0,
            inclination=0.7,
            phi_ref=0.9,
            f_min=fmin,
            f_ref=fref,
            delta_f=df,
            f_max=fmax,
            lal_flags=dict(NO_MULTIBANDING),
            fd_modes_stride=8,
            **extra,
        )
        for m1, m2, s1, s2, fmin, fref, df, fmax, extra in configs
    ]


def xpnr_sequence_mode_cases():
    """The J-frame modes of SimIMRPhenomXPHMFrequencySequenceOneMode with
    the flags of IMRPhenomXPNR (FrequencyDomainModel::modes)."""
    configs = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, mode_array
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 25.0, None),
        (30.0, 20.0, [5e-8, 0.0, 0.5], [0.0, 0.0, -0.3], 20.0, 20.0, None),
        (60.0, 12.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0, [[2, 1], [2, -1], [3, 3], [3, -3]]),
    ]
    return [
        dict(
            approximant="IMRPhenomXPNR",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=700.0,
            inclination=1.3,
            phi_ref=0.5,
            f_min=fmin,
            f_ref=fref,
            f_max=f_max_imr(0.3)(m1 + m2),
            lal_flags=dict(XPNR_FLAGS, **NO_MULTIBANDING),
            jframe_sequence=True,
            **(dict(mode_array=modes) if modes else dict()),
        )
        for m1, m2, s1, s2, fmin, fref, modes in configs
    ]


def nrsur7dq4_cases():
    """NRSur7dq4 polarizations from SimInspiralChooseTDWaveform: mass ratios
    1 to 6, spins up to 0.99, f_ref != f_min, f_min = 0 (the whole
    surrogate), m1 < m2 (label switch), and LAL's error cases."""
    generic_99 = unit([0.6, -0.6, 0.53], 0.99)
    configs = [
        # mass_1, mass_2, spin_1, spin_2, distance, inclination, phi_ref,
        # f_min, f_ref, delta_t
        (40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 100.0, 0.7, 0.1, 20.0, 20.0, 1 / 4096),
        (30.0, 30.0, [0.0, 0.0, 0.5], [0.0, 0.0, 0.5], 400.0, 0.4, 1.3, 20.0, 20.0, 1 / 2048),
        (30.0, 30.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 400.0, 1.0, 0.0, 20.0, 0.0, 1 / 2048),
        (75.0, 25.0, [0.6, 0.3, -0.5], [-0.2, 0.7, 0.1], 800.0, 2.1, 4.0, 15.0, 25.0, 1 / 2048),
        (96.0, 24.0, [0.7, -0.5, 0.4], [0.3, 0.3, -0.8], 1500.0, 1.2, 2.5, 12.0, 20.0, 1 / 1024),
        (120.0, 20.0, generic_99, unit([-0.7, 0.3, 0.6], 0.99), 2000.0, 0.9, 5.2, 10.0, 12.0, 1 / 1024),
        (150.0, 50.0, [0.0, 0.0, 0.99], [0.0, 0.0, -0.99], 1000.0, 0.3, 0.8, 10.0, 10.0, 1 / 1024),
        (20.0, 45.0, [0.3, 0.4, -0.2], [-0.5, 0.1, 0.6], 500.0, 1.4, 3.3, 20.0, 30.0, 1 / 4096),
        (50.0, 30.0, [0.2, -0.6, 0.3], [0.5, 0.2, 0.1], 600.0, 0.0, 0.0, 0.0, 0.0, 1 / 2048),
        (50.0, 30.0, [0.2, -0.6, 0.3], [0.5, 0.2, 0.1], 600.0, 2.0, 1.0, 0.0, 30.0, 1 / 2048),
        (200.0, 100.0, [-0.5, 0.5, 0.1], [0.0, -0.8, 0.3], 2500.0, 2.5, 1.0, 0.0, 0.0, 1 / 512),
        (14.0, 6.0, [0.1, 0.8, 0.3], [0.4, -0.4, -0.4], 300.0, 3.1, 0.5, 60.0, 60.0, 1 / 8192),
        (36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 500.0, 0.7, 1.2, 0.0, "start", 1 / 2048),
        (27.0, 30.0, [0.0, 0.0, 0.0], generic_99, 700.0, 1.9, 2.2, 0.0, "start", 1 / 2048),
        (45.0, 15.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 300.0, np.pi, 0.3, 25.0, 25.0, 1 / 2048),
    ]
    start_factors = iter([1.002, 1.004])
    cases = list()
    for m1, m2, s1, s2, dl, iota, phi, fmin, fref, dt in configs:
        case = dict(
            approximant="NRSur7dq4",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=dl,
            inclination=iota,
            phi_ref=phi,
            f_min=fmin,
            f_ref=fref,
            delta_t=dt,
        )
        if fref == "start":
            case["f_ref_factor"] = next(start_factors)
        cases.append(case)
    errors = [
        # f_min below the surrogate start, q > 6.01, f_ref < f_min,
        # omega_ref > 0.201, f_ref below the surrogate start, |chi| > 1
        (40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 5.0, 20.0),
        (70.0, 10.0, [0.1, 0.0, 0.0], [0.0, 0.0, 0.1], 20.0, 20.0),
        (40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 30.0, 20.0),
        (40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 20.0, 400.0),
        (40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 0.0, 5.0),
        (40.0, 20.0, [0.6, -0.6, 0.6], [0.1, 0.5, -0.3], 20.0, 20.0),
    ]
    # coorbital mode arrays and unlimited extrapolation
    options = [
        dict(mode_array=[(2, 2), (2, -2)]),
        dict(mode_array=[(ell, m) for ell in (2, 3) for m in range(-ell, ell + 1)]),
        dict(mode_array=[(2, 0), (3, 1), (4, -4)]),
    ]
    for option in options:
        cases.append(
            dict(
                approximant="NRSur7dq4",
                mass_1=40.0,
                mass_2=20.0,
                spin_1=[0.4, -0.3, 0.2],
                spin_2=[0.1, 0.5, -0.3],
                luminosity_distance=100.0,
                inclination=0.7,
                phi_ref=0.1,
                f_min=20.0,
                f_ref=20.0,
                delta_t=1 / 2048,
                **option,
            )
        )
    cases.append(
        dict(
            approximant="NRSur7dq4",
            mass_1=70.0,
            mass_2=10.0,
            spin_1=[0.1, 0.0, 0.0],
            spin_2=[0.0, 0.0, 0.1],
            luminosity_distance=100.0,
            inclination=0.5,
            phi_ref=0.0,
            f_min=20.0,
            f_ref=20.0,
            delta_t=1 / 2048,
            unlimited_extrapolation=True,
        )
    )
    for m1, m2, s1, s2, fmin, fref in errors:
        cases.append(
            dict(
                approximant="NRSur7dq4",
                mass_1=m1,
                mass_2=m2,
                spin_1=s1,
                spin_2=s2,
                luminosity_distance=100.0,
                inclination=0.5,
                phi_ref=0.0,
                f_min=fmin,
                f_ref=fref,
                delta_t=1 / 1024,
            )
        )
    return cases


def nrsur7dq4_mode_cases():
    """NRSur7dq4 modes from SimInspiralChooseTDModes (every other sample)."""
    configs = [
        # mass_1, mass_2, spin_1, spin_2, f_min, f_ref, delta_t
        (40.0, 20.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 20.0, 20.0, 1 / 1024),
        (20.0, 45.0, [0.3, 0.4, -0.2], [-0.5, 0.1, 0.6], 20.0, 30.0, 1 / 1024),
        (50.0, 30.0, [0.2, -0.6, 0.3], [0.5, 0.2, 0.1], 0.0, 0.0, 1 / 1024),
        (100.0, 25.0, [0.5, 0.1, 0.8], [-0.3, -0.3, 0.2], 12.0, 15.0, 1 / 512),
    ]
    return [
        dict(
            approximant="NRSur7dq4",
            mass_1=m1,
            mass_2=m2,
            spin_1=s1,
            spin_2=s2,
            luminosity_distance=400.0,
            f_min=fmin,
            f_ref=fref,
            delta_t=dt,
            sample_stride=2,
        )
        for m1, m2, s1, s2, fmin, fref, dt in configs
    ] + [
        dict(
            approximant="NRSur7dq4",
            mass_1=40.0,
            mass_2=20.0,
            spin_1=[0.4, -0.3, 0.2],
            spin_2=[0.1, 0.5, -0.3],
            luminosity_distance=400.0,
            f_min=20.0,
            f_ref=20.0,
            delta_t=1 / 1024,
            sample_stride=2,
            mode_array=[(2, 2), (2, -2), (3, 3), (3, -3)],
        )
    ]


def nrsur7dq4_dynamics_cases():
    """XLALPrecessingNRSurDynamics with general quaternions and phases."""
    norm = np.linalg.norm([0.9, 0.1, -0.3, 0.2])
    configs = [
        # q, chi_a, chi_b, omega_ref, init_quat, init_orbphase
        (2.0, [0.4, -0.3, 0.2], [0.1, 0.5, -0.3], 0.02, [1.0, 0.0, 0.0, 0.0], 0.0),
        (4.5, [0.7, 0.1, -0.6], [-0.2, 0.2, 0.9], 0.0, [x / norm for x in [0.9, 0.1, -0.3, 0.2]], 1.2),
        (1.0, [0.0, 0.0, 0.0], [0.3, 0.3, 0.3], 0.05, [1.0, 0.0, 0.0, 0.0], 2.5),
        (5.9, unit([0.2, 0.9, -0.1], 0.99), [0.0, 0.0, -0.5], 0.03, [1.0, 0.0, 0.0, 0.0], 0.4),
    ]
    return [
        dict(
            approximant="NRSur7dq4",
            q=q,
            chi_a=chi_a,
            chi_b=chi_b,
            omega_ref=omega,
            init_quat=quat,
            init_orbphase=phase,
        )
        for q, chi_a, chi_b, omega, quat, phase in configs
    ]


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
    imrphenomxhm_122019_modes=xhm_122019_mode_cases,
    imrphenomxhm_122019=xhm_122019_cases,
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
    choose_fd_modes=choose_fd_modes_cases,
    imrphenomxpnr=xpnr_cases,
    imrphenomxpnr_uniform=xpnr_uniform_cases,
    imrphenomxpnr_modes=xpnr_mode_cases,
    imrphenomxpnr_sequence_modes=xpnr_sequence_mode_cases,
    nrsur7dq4=nrsur7dq4_cases,
    nrsur7dq4_modes=nrsur7dq4_mode_cases,
    nrsur7dq4_dynamics=nrsur7dq4_dynamics_cases,
    spin_weighted_harmonics=spin_weighted_harmonic_cases,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-frequencies", type=int, default=256)
    parser.add_argument("--suites", nargs="*", default=list(SUITES))
    parser.add_argument("--outdir", default=str(DATA_DIR))
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name in args.suites:
        cases = [generate(case, args.n_frequencies) for case in SUITES[name]()]
        output = dict(
            lalsimulation_version=lalsim.__version__,
            lal_version=lal.__version__,
            cases=cases,
        )
        path = outdir / f"{name}.json"
        with open(path, "w") as ff:
            json.dump(output, ff)
        print(f"Wrote {len(cases)} cases to {path}")


if __name__ == "__main__":
    main()
