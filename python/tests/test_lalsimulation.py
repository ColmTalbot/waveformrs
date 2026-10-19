"""Agreement with lalsimulation's Python API at the crate's tolerances.

The crate's tests compare against stored LALSimulation 6.2.0 output to a
relative error (max |difference| / max |LAL|) of 1e-10 (TaylorF2, IMRPhenomD)
or 1e-9 (IMRPhenomX) and a mismatch of 1e-12 on the conda-forge macOS arm64
build that it mirrors. Elsewhere the higher-mode models are compared with the
crate's cross-build tolerances (docs/linux.md).
"""

import numpy as np
import pytest

import waveformrs as wf
from conftest import MIRRORED_BUILD, compare, requires_nrsur

BINARY = dict(
    mass_1=36.0,
    mass_2=29.0,
    luminosity_distance=500.0,
    inclination=0.7,
    phi_ref=1.2,
    f_ref=20.0,
)
ALIGNED = dict(BINARY, spin_1=(0.0, 0.0, 0.4), spin_2=(0.0, 0.0, -0.3))
PRECESSING = dict(BINARY, spin_1=(0.4, 0.2, 0.3), spin_2=(-0.2, 0.3, -0.1))
UNEQUAL = dict(
    mass_1=12.0,
    mass_2=48.0,
    spin_1=(-0.3, 0.5, 0.2),
    spin_2=(0.6, -0.1, 0.4),
    luminosity_distance=800.0,
    inclination=2.1,
    phi_ref=-0.4,
    f_ref=25.0,
)
FREQUENCIES = np.geomspace(20.0, 1024.0, 400)

# (approximant, binary, options, (relative error, mismatch))
FD_CASES = [
    ("TaylorF2", dict(ALIGNED, mass_1=1.6, mass_2=1.3, f_ref=0.0), {}, (1e-10, 1e-12)),
    ("TaylorF2", dict(ALIGNED, mass_1=10.0, mass_2=8.0), {}, (1e-10, 1e-12)),
    ("IMRPhenomD", ALIGNED, {}, (1e-10, 1e-12)),
    ("IMRPhenomD", dict(ALIGNED, mass_1=80.0, f_ref=0.0), {}, (1e-10, 1e-12)),
    ("IMRPhenomXAS", ALIGNED, {}, (1e-9, 1e-12)),
    ("IMRPhenomXHM", ALIGNED, {}, (1e-9, 1e-12)),
    ("IMRPhenomXHM", dict(ALIGNED, mass_1=10.0, mass_2=50.0), {}, (1e-9, 1e-12)),
    ("IMRPhenomXHM", ALIGNED, dict(mode_array=[(2, 2), (2, -2), (3, 3), (3, -3)]), (1e-9, 1e-12)),
    ("IMRPhenomXP", PRECESSING, dict(precession_version=320), (1e-9, 1e-12)),
    ("IMRPhenomXP", PRECESSING, dict(precession_version=223), (1e-9, 1e-12)),
    ("IMRPhenomXPHM", PRECESSING, dict(precession_version=102), (1e-9, 1e-12)),
    ("IMRPhenomXPHM", PRECESSING, dict(precession_version=223), (1e-9, 1e-12)),
    ("IMRPhenomXPHM", PRECESSING, dict(precession_version=320), (1e-9, 1e-12)),
    ("IMRPhenomXPHM", UNEQUAL, dict(precession_version=102), (1e-9, 1e-12)),
    ("IMRPhenomXPHM", UNEQUAL, dict(precession_version=223), (1e-9, 1e-12)),
    ("IMRPhenomXPHM", UNEQUAL, dict(precession_version=320), (1e-9, 1e-12)),
    (
        "IMRPhenomXPHM",
        PRECESSING,
        dict(precession_version=320, mode_array=[(2, 2), (3, 3)]),
        (1e-9, 1e-12),
    ),
    # LALSimulation 6.2.0 has no parity fix (lalsuite !2593, test_xpnr.py)
    ("IMRPhenomXPNR", PRECESSING, dict(pnr_parity_fix=False), (1e-9, 1e-12)),
    ("IMRPhenomXPNR", UNEQUAL, dict(pnr_parity_fix=False), (1e-9, 1e-12)),
]

# the crate's tolerances where it does not mirror the reference build
CROSS_BUILD = {
    "IMRPhenomXHM": (1e-3, 1e-6),
    "IMRPhenomXPHM": (1e-3, 1e-6),
    "IMRPhenomXP": (1e-7, 1e-12),
    "IMRPhenomXPNR": (1e-3, 1e-6),
}


def tolerances(approximant, strict):
    if MIRRORED_BUILD:
        return strict
    return CROSS_BUILD.get(approximant, strict)


def case_id(case):
    approximant, binary, options, _ = case
    label = f"{approximant}-m{binary['mass_1']:g}+{binary['mass_2']:g}"
    return label + "".join(f"-{key}={value}" for key, value in options.items())


def lal_args(lal, binary):
    """The binary in the units and order of SimInspiralChooseFDWaveform."""
    s1 = binary.get("spin_1", (0.0, 0.0, 0.0))
    s2 = binary.get("spin_2", (0.0, 0.0, 0.0))
    return (
        binary["mass_1"] * lal.MSUN_SI,
        binary["mass_2"] * lal.MSUN_SI,
        *s1,
        *s2,
        binary["luminosity_distance"] * 1e6 * lal.PC_SI,
        binary["inclination"],
    )


def lal_dict(lal, lalsim, options):
    """A LALDict with the options and multibanding off."""
    params = lal.CreateDict()
    lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(params, 0)
    lalsim.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(params, 0)
    if "precession_version" in options:
        lalsim.SimInspiralWaveformParamsInsertPhenomXPrecVersion(
            params, options["precession_version"]
        )
    if options.get("mode_array") is not None:
        modes = lalsim.SimInspiralCreateModeArray()
        for ell, emm in options["mode_array"]:
            lalsim.SimInspiralModeArrayActivateMode(modes, ell, emm)
        lalsim.SimInspiralWaveformParamsInsertModeArray(params, modes)
    return params


def lal_sequence(lal, lalsim, approximant, binary, options, frequencies):
    """SimInspiralChooseFDWaveformSequence."""
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, inclination = lal_args(lal, binary)
    vector = lal.CreateREAL8Vector(len(frequencies))
    vector.data = frequencies
    hp, hc = lalsim.SimInspiralChooseFDWaveformSequence(
        binary["phi_ref"],
        m1,
        m2,
        s1x,
        s1y,
        s1z,
        s2x,
        s2y,
        s2z,
        binary["f_ref"],
        distance,
        inclination,
        lal_dict(lal, lalsim, options),
        lalsim.GetApproximantFromString(approximant),
        vector,
    )
    return hp.data.data, hc.data.data


@pytest.mark.parametrize("case", FD_CASES, ids=[case_id(case) for case in FD_CASES])
def test_polarizations(lal, case):
    approximant, binary, options, strict = case
    lal_mod, lalsim = lal
    rel, mismatch = tolerances(approximant, strict)
    hp, hc = wf.polarizations(approximant, FREQUENCIES, **binary, **options)
    lp, lc = lal_sequence(lal_mod, lalsim, approximant, binary, options, FREQUENCIES)
    for test, reference in ((hp, lp), (hc, lc)):
        error, mm = compare(test, reference)
        assert error <= rel and mm <= mismatch, (error, mm)


@pytest.mark.parametrize(
    "case",
    [case for case in FD_CASES if case[0] in ("IMRPhenomXAS", "IMRPhenomXHM", "IMRPhenomXPHM")],
    ids=case_id,
)
def test_polarizations_on_grid(lal, case):
    """The fast grid paths, on a uniform grid, against the sequence of LAL."""
    approximant, binary, options, strict = case
    lal_mod, lalsim = lal
    rel, mismatch = tolerances(approximant, strict)
    grid = wf.FrequencyGrid.uniform(20.0, 1024.0, 0.5)
    hp, hc = wf.polarizations_on_grid(approximant, grid, **binary, **options)
    lp, lc = lal_sequence(lal_mod, lalsim, approximant, binary, options, grid.frequencies)
    for test, reference in ((hp, lp), (hc, lc)):
        error, mm = compare(test, reference)
        assert error <= rel and mm <= mismatch, (error, mm)


def test_options_object_matches_lal_defaults(lal):
    """WaveformOptions.lalsimulation_defaults() is LAL with an empty LALDict
    (apart from multibanding)."""
    lal_mod, lalsim = lal
    options = wf.WaveformOptions.lalsimulation_defaults()
    hp, _ = wf.polarizations("IMRPhenomXP", FREQUENCIES, options=options, **PRECESSING)
    lp, _ = lal_sequence(lal_mod, lalsim, "IMRPhenomXP", PRECESSING, {}, FREQUENCIES)
    error, mm = compare(hp, lp)
    assert error <= tolerances("IMRPhenomXP", (1e-9, 1e-12))[0], error


CHOOSE_FD_MODES_CASES = [
    ("IMRPhenomXHM", ALIGNED, {}),
    ("IMRPhenomXPHM", PRECESSING, dict(precession_version=223)),
    ("IMRPhenomXPHM", PRECESSING, dict(precession_version=320)),
]


@pytest.mark.parametrize(
    "case", CHOOSE_FD_MODES_CASES, ids=[case_id((*case, None)) for case in CHOOSE_FD_MODES_CASES]
)
def test_choose_fd_modes(lal, case):
    approximant, binary, options = case
    lal_mod, lalsim = lal
    delta_f, f_min, f_max = 0.25, 20.0, 512.0
    result = wf.choose_fd_modes(approximant, delta_f, f_min, f_max, **binary, **options)
    frequencies, modes = result
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, inclination = lal_args(lal_mod, binary)
    node = lalsim.SimInspiralChooseFDModes(
        m1,
        m2,
        s1x,
        s1y,
        s1z,
        s2x,
        s2y,
        s2z,
        delta_f,
        f_min,
        f_max,
        binary["f_ref"],
        binary["phi_ref"],
        distance,
        inclination,
        lal_dict(lal_mod, lalsim, options),
        lalsim.GetApproximantFromString(approximant),
    )
    lal_modes = {}
    while node is not None:
        lal_modes[(node.l, node.m)] = node.mode
        node = node.next
    assert sorted(lal_modes) == sorted(modes)
    for key, series in lal_modes.items():
        reference = series.data.data
        assert len(reference) == len(frequencies)
        # LAL's two-sided series: sample i at (i - (n - 1) / 2) delta_f
        half = (len(reference) - 1) // 2
        np.testing.assert_array_equal(
            frequencies, (np.arange(len(reference)) - half) * series.deltaF
        )
        # the l = 3 modes of the higher-mode models amplify rounding off the
        # mirrored build (the crate's choose_fd_modes_l3 tolerances)
        rel, mismatch = (1e-9, 1e-12) if MIRRORED_BUILD else ((2e-2, 4e-5) if key[0] == 3 else (1e-8, 1e-12))
        error, mm = compare(modes[key], reference)
        assert error <= rel and mm <= mismatch, (key, error, mm)


def test_modes_polarizations_round_trip():
    """FrequencyModes.polarizations() rebuilds the polarizations."""
    for approximant, binary, options, rtol in (
        ("IMRPhenomXAS", ALIGNED, {}, 1e-13),
        ("IMRPhenomXPHM", PRECESSING, dict(precession_version=102), 1e-12),
    ):
        hp, hc = wf.polarizations(approximant, FREQUENCIES, **binary, **options)
        result = wf.modes(approximant, FREQUENCIES, **binary, **options)
        mp, mc = result.polarizations()
        assert compare(mp, hp)[0] < rtol
        assert compare(mc, hc)[0] < rtol


def lal_td(lal, lalsim, approximant, binary, delta_t, f_min, options=None):
    """SimInspiralChooseTDWaveform."""
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, inclination = lal_args(lal, binary)
    return lalsim.SimInspiralChooseTDWaveform(
        m1,
        m2,
        s1x,
        s1y,
        s1z,
        s2x,
        s2y,
        s2z,
        distance,
        inclination,
        binary["phi_ref"],
        0.0,
        0.0,
        0.0,
        delta_t,
        f_min,
        binary["f_ref"],
        lal_dict(lal, lalsim, options or {}),
        lalsim.GetApproximantFromString(approximant),
    )


def check_series(result, lal_hp, lal_hc, tolerance):
    """Same epoch, sampling and length; values to `tolerance` of the peak."""
    assert result.epoch == float(lal_hp.epoch)
    assert len(result.hplus) == lal_hp.data.length
    for test, reference in ((result.hplus, lal_hp.data.data), (result.hcross, lal_hc.data.data)):
        error = np.max(np.abs(test - reference)) / np.max(np.abs(reference))
        assert error <= tolerance, error


NRSUR_BINARY = dict(
    mass_1=40.0,
    mass_2=20.0,
    spin_1=(0.4, -0.3, 0.2),
    spin_2=(0.1, 0.5, -0.3),
    luminosity_distance=100.0,
    inclination=0.7,
    phi_ref=0.1,
    f_ref=20.0,
)


@requires_nrsur
@pytest.mark.parametrize("f_min", [20.0, 0.0])
def test_nrsur7dq4_td_polarizations(lal, f_min):
    lal_mod, lalsim = lal
    delta_t = 1.0 / 4096
    result = wf.td_polarizations("NRSur7dq4", delta_t, f_min, **NRSUR_BINARY)
    lp, lc = lal_td(lal_mod, lalsim, "NRSur7dq4", NRSUR_BINARY, delta_t, f_min)
    assert result.delta_t == lp.deltaT
    check_series(result, lp, lc, 1e-12)
    # again, from the cached modes of the thread's workspace
    again = wf.td_polarizations("NRSur7dq4", delta_t, f_min, **dict(NRSUR_BINARY, inclination=1.3))
    lp, lc = lal_td(
        lal_mod, lalsim, "NRSur7dq4", dict(NRSUR_BINARY, inclination=1.3), delta_t, f_min
    )
    check_series(again, lp, lc, 1e-12)


@requires_nrsur
def test_nrsur7dq4_td_modes(lal):
    lal_mod, lalsim = lal
    delta_t = 1.0 / 2048
    b = NRSUR_BINARY
    result = wf.td_modes("NRSur7dq4", delta_t, 20.0, **b)
    times, modes = result
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, _ = lal_args(lal_mod, b)
    node = lalsim.SimInspiralChooseTDModes(
        b["phi_ref"],
        delta_t,
        m1,
        m2,
        s1x,
        s1y,
        s1z,
        s2x,
        s2y,
        s2z,
        20.0,
        b["f_ref"],
        distance,
        None,
        4,
        lalsim.NRSur7dq4,
    )
    count = 0
    while node is not None:
        reference = node.mode.data.data
        assert result.epoch == float(node.mode.epoch)
        error = np.max(np.abs(modes[node.l, node.m] - reference)) / np.max(np.abs(reference))
        assert error <= 1e-12, ((node.l, node.m), error)
        count += 1
        node = node.next
    assert count == len(modes)
    assert len(times) == len(modes[2, 2])


def test_imrphenomxas_td_polarizations(lal):
    """SimInspiralChooseTDWaveform of a frequency-domain model."""
    lal_mod, lalsim = lal
    delta_t = 1.0 / 2048
    result = wf.td_polarizations("IMRPhenomXAS", delta_t, 20.0, **ALIGNED)
    lp, lc = lal_td(lal_mod, lalsim, "IMRPhenomXAS", ALIGNED, delta_t, 20.0)
    check_series(result, lp, lc, 1e-10)


def lal_simulate_td(lal, lalsim, approximant, binary, delta_t, f_min, options):
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, inclination = lal_args(lal, binary)
    return lalsim.SimInspiralTD(
        m1,
        m2,
        s1x,
        s1y,
        s1z,
        s2x,
        s2y,
        s2z,
        distance,
        inclination,
        binary["phi_ref"],
        0.0,
        0.0,
        0.0,
        delta_t,
        f_min,
        binary["f_ref"],
        lal_dict(lal, lalsim, options),
        lalsim.GetApproximantFromString(approximant),
    )


@pytest.mark.parametrize(
    "approximant, binary, options",
    [
        ("IMRPhenomXAS", ALIGNED, {}),
        ("IMRPhenomXPHM", PRECESSING, dict(precession_version=320)),
        ("TaylorF2", dict(ALIGNED, mass_1=10.0, mass_2=8.0), {}),
    ],
)
def test_simulate_td(lal, approximant, binary, options):
    lal_mod, lalsim = lal
    delta_t = 1.0 / 2048
    result = wf.simulate_td(approximant, delta_t, 20.0, **binary, **options)
    lp, lc = lal_simulate_td(lal_mod, lalsim, approximant, binary, delta_t, 20.0, options)
    tolerance = 1e-10 if MIRRORED_BUILD or approximant in ("IMRPhenomXAS", "TaylorF2") else 1e-3
    check_series(result, lp, lc, tolerance)


@requires_nrsur
@pytest.mark.parametrize("delta_f", [0.0, 0.125])
def test_simulate_fd_nrsur7dq4(lal, delta_f):
    lal_mod, lalsim = lal
    b = NRSUR_BINARY
    result = wf.simulate_fd("NRSur7dq4", delta_f, 20.0, 2048.0, **b)
    m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, distance, inclination = lal_args(lal_mod, b)
    lp, lc = lalsim.SimInspiralFD(
        m1,
        m2,
        s1x,
        s1y,
        s1z,
        s2x,
        s2y,
        s2z,
        distance,
        inclination,
        b["phi_ref"],
        0.0,
        0.0,
        0.0,
        delta_f,
        20.0,
        2048.0,
        b["f_ref"],
        None,
        lalsim.NRSur7dq4,
    )
    assert result.delta_f == lp.deltaF
    assert result.epoch == float(lp.epoch)
    assert len(result.hplus) == lp.data.length
    np.testing.assert_array_equal(result.frequencies, np.arange(lp.data.length) * lp.deltaF)
    for test, reference in ((result.hplus, lp.data.data), (result.hcross, lc.data.data)):
        error = np.max(np.abs(test - reference)) / np.max(np.abs(reference))
        assert error <= 1e-10, error
