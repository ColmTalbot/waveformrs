"""Batches: independent of the number of threads, equal to single calls, with
per-binary failures reported (or raised with ``strict=True``)."""

import numpy as np
import pytest

import waveformrs as wf
from conftest import requires_nrsur


def draw(n, precessing=True, seed=1):
    rng = np.random.default_rng(seed)
    mass_1 = rng.uniform(10.0, 60.0, n)
    mass_2 = rng.uniform(10.0, 60.0, n)
    spins = rng.uniform(-0.5, 0.5, (2, n, 3))
    if not precessing:
        spins[..., :2] = 0.0
    return dict(
        mass_1=mass_1,
        mass_2=mass_2,
        spin_1=spins[0],
        spin_2=spins[1],
        luminosity_distance=rng.uniform(100.0, 1000.0, n),
        inclination=rng.uniform(0.0, np.pi, n),
        phi_ref=rng.uniform(0.0, 2 * np.pi, n),
        f_ref=20.0,
    )


def binary(parameters, i):
    """Binary i of a dict of arrays as keyword arguments."""
    return {
        key: (value[i] if np.ndim(value) > 0 else value) for key, value in parameters.items()
    }


GRID = wf.FrequencyGrid.uniform(20.0, 512.0, 0.5)
CASES = [
    ("IMRPhenomXAS", False, {}),
    ("IMRPhenomXHM", False, {}),
    ("IMRPhenomXPHM", True, dict(precession_version=320)),
    ("IMRPhenomXPHM", True, dict(precession_version=223)),
]


@pytest.mark.parametrize("approximant, precessing, options", CASES)
def test_threads_and_single_calls(approximant, precessing, options):
    parameters = draw(12, precessing)
    serial = wf.batch_polarizations(approximant, GRID, parameters, threads=1, **options)
    assert serial.ok.all() and not serial.errors
    assert serial.hplus.shape == (12, len(GRID))
    assert serial.hplus.dtype == np.complex128
    for threads in (2, 3, 8, 0):
        result = wf.batch_polarizations(approximant, GRID, parameters, threads=threads, **options)
        np.testing.assert_array_equal(result.hplus, serial.hplus)
        np.testing.assert_array_equal(result.hcross, serial.hcross)
    for i in range(12):
        hp, hc = wf.polarizations_on_grid(approximant, GRID, **binary(parameters, i), **options)
        np.testing.assert_array_equal(serial.hplus[i], hp)
        np.testing.assert_array_equal(serial.hcross[i], hc)


def test_exact_matches_polarizations():
    """exact=True evaluates as polarizations() (the frequency sequence)."""
    parameters = draw(6)
    result = wf.batch_polarizations(
        "IMRPhenomXPHM", GRID, parameters, threads=3, exact=True, precession_version=102
    )
    for i in range(6):
        hp, hc = wf.polarizations(
            "IMRPhenomXPHM", GRID.frequencies, **binary(parameters, i), precession_version=102
        )
        np.testing.assert_array_equal(result.hplus[i], hp)
        np.testing.assert_array_equal(result.hcross[i], hc)


def test_parameter_forms():
    """A dict with spin vectors or components, a structured array and a
    table give the same binaries."""
    parameters = draw(5)
    table = wf.parameter_table(parameters)
    assert table.shape == (5, 12) and table.flags.c_contiguous
    components = {k: v for k, v in parameters.items() if not k.startswith("spin")}
    for i, axis in enumerate("xyz"):
        components[f"spin_1{axis}"] = parameters["spin_1"][:, i]
        components[f"spin_2{axis}"] = parameters["spin_2"][:, i]
    np.testing.assert_array_equal(wf.parameter_table(components), table)
    structured = np.zeros(
        5,
        dtype=[
            ("mass_1", "f8"),
            ("mass_2", "f8"),
            ("spin_1", "f8", 3),
            ("spin_2", "f8", 3),
            ("luminosity_distance", "f8"),
            ("inclination", "f8"),
            ("phi_ref", "f8"),
            ("f_ref", "f8"),
        ],
    )
    for key, value in parameters.items():
        structured[key] = value
    np.testing.assert_array_equal(wf.parameter_table(structured), table)
    np.testing.assert_array_equal(wf.parameter_table(table), table)
    # defaults and broadcasting
    defaults = wf.parameter_table(dict(mass_1=[30.0, 40.0], mass_2=20.0))
    np.testing.assert_array_equal(
        defaults,
        [
            [30.0, 20.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0],
            [40.0, 20.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0],
        ],
    )
    a = wf.batch_polarizations("IMRPhenomXP", GRID, parameters)
    b = wf.batch_polarizations("IMRPhenomXP", GRID, structured)
    np.testing.assert_array_equal(a.hplus, b.hplus)
    with pytest.raises(TypeError, match="unknown binary parameters"):
        wf.parameter_table(dict(mass_1=1.0, mass_2=1.0, mass1=2.0))
    with pytest.raises(TypeError, match="missing binary parameter"):
        wf.parameter_table(dict(mass_1=1.0))
    with pytest.raises(TypeError, match="not both"):
        wf.parameter_table(dict(mass_1=1.0, mass_2=1.0, spin_1=(0, 0, 0), spin_1z=0.1))


def test_failures_are_reported():
    parameters = draw(6, precessing=False)
    parameters["mass_1"][[1, 4]] = -5.0
    parameters["spin_1"][2] = [0.3, 0.0, 0.1]  # in-plane spin: rejected by XAS
    result = wf.batch_polarizations("IMRPhenomXAS", GRID, parameters, threads=2)
    np.testing.assert_array_equal(result.ok, [True, False, False, True, False, True])
    assert sorted(result.errors) == [1, 2, 4]
    assert isinstance(result.errors[1], wf.InvalidParameterError)
    assert "transverse" in str(result.errors[2])
    assert not result.hplus[[1, 2, 4]].any() and not result.hcross[[1, 2, 4]].any()
    hp, _ = wf.polarizations_on_grid("IMRPhenomXAS", GRID, **binary(parameters, 0))
    np.testing.assert_array_equal(result.hplus[0], hp)
    with pytest.raises(wf.InvalidParameterError, match="^binary 1: "):
        wf.batch_polarizations("IMRPhenomXAS", GRID, parameters, threads=2, strict=True)


def test_batch_rejects_time_domain_models():
    with pytest.raises(wf.UnsupportedError, match="batch_td_polarizations"):
        wf.batch_polarizations("NRSur7dq4", GRID, draw(2))


def test_batch_accepts_frequencies():
    parameters = draw(3, precessing=False)
    a = wf.batch_polarizations("IMRPhenomD", GRID, parameters)
    b = wf.batch_polarizations("IMRPhenomD", GRID.frequencies, parameters)
    np.testing.assert_array_equal(a.hplus, b.hplus)
    hplus, hcross = a
    assert hplus is a.hplus and hcross is a.hcross


@requires_nrsur
def test_nrsur7dq4_batch():
    rng = np.random.default_rng(3)
    n = 6
    parameters = dict(
        mass_1=rng.uniform(40.0, 60.0, n),
        mass_2=rng.uniform(15.0, 40.0, n),
        spin_1=rng.uniform(-0.4, 0.4, (n, 3)),
        spin_2=rng.uniform(-0.4, 0.4, (n, 3)),
        luminosity_distance=400.0,
        inclination=rng.uniform(0.0, np.pi, n),
        phi_ref=rng.uniform(0.0, 2 * np.pi, n),
        f_ref=20.0,
    )
    parameters["mass_2"][3] = 1.0  # mass ratio 50: too much extrapolation
    delta_t = 1.0 / 2048
    serial = wf.batch_td_polarizations("NRSur7dq4", delta_t, 20.0, parameters)
    assert len(serial) == n
    assert serial[3] is None and not serial.ok[3]
    assert isinstance(serial.errors[3], wf.InvalidParameterError)
    for threads in (2, 4):
        result = wf.batch_td_polarizations("NRSur7dq4", delta_t, 20.0, parameters, threads=threads)
        for a, b in zip(result, serial):
            if b is None:
                assert a is None
                continue
            assert (a.epoch, a.delta_t) == (b.epoch, b.delta_t)
            np.testing.assert_array_equal(a.hplus, b.hplus)
            np.testing.assert_array_equal(a.hcross, b.hcross)
    for i, waveform in enumerate(serial):
        if waveform is None:
            continue
        single = wf.td_polarizations("NRSur7dq4", delta_t, 20.0, **binary(parameters, i))
        assert single.epoch == waveform.epoch
        np.testing.assert_array_equal(single.hplus, waveform.hplus)
        np.testing.assert_array_equal(single.hcross, waveform.hcross)
    with pytest.raises(wf.InvalidParameterError, match="^binary 3: "):
        wf.batch_td_polarizations("NRSur7dq4", delta_t, 20.0, parameters, strict=True)


def test_td_batch_of_frequency_domain_model():
    parameters = draw(3, precessing=False)
    result = wf.batch_td_polarizations("IMRPhenomXAS", 1.0 / 2048, 20.0, parameters, threads=2)
    for i, waveform in enumerate(result):
        single = wf.td_polarizations("IMRPhenomXAS", 1.0 / 2048, 20.0, **binary(parameters, i))
        np.testing.assert_array_equal(single.hplus, waveform.hplus)
        assert single.epoch == waveform.epoch
