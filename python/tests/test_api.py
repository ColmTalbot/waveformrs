"""The Python interface: approximants, options, inputs and result classes."""

import inspect
import pickle

import numpy as np
import pytest

import waveformrs as wf
from conftest import requires_nrsur

F = np.linspace(20.0, 512.0, 200)
BINARY = dict(mass_1=30.0, mass_2=20.0, spin_1=(0.0, 0.0, 0.3), inclination=0.4, f_ref=20.0)
PRECESSING = dict(BINARY, spin_1=(0.3, -0.2, 0.3), spin_2=(0.1, 0.4, -0.2))


def test_approximants():
    names = wf.approximants()
    assert names[:3] == ["TaylorF2", "IMRPhenomD", "IMRPhenomXAS"]
    assert "IMRPhenomXPHM" in names and "NRSur7dq4" in names
    assert wf.approximants("time") == ["NRSur7dq4"]
    assert set(wf.approximants("frequency")) == set(names) - {"NRSur7dq4"}
    assert wf.is_precessing("IMRPhenomXPHM") and not wf.is_precessing("IMRPhenomXHM")
    assert wf.is_frequency_domain("IMRPhenomD") and not wf.is_frequency_domain("NRSur7dq4")
    with pytest.raises(ValueError):
        wf.approximants("both")


def test_waveform_options():
    default = wf.WaveformOptions()
    assert default.precession_version == 320
    assert default.mode_array is None
    assert default.version_fallback is True
    assert wf.WaveformOptions.lalsimulation_defaults().precession_version == 223
    options = wf.WaveformOptions(precession_version="msa", mode_array=[[2, 2], (3, -3)])
    assert options.precession_version == 223
    assert options.mode_array == [(2, 2), (3, -3)]
    assert options == wf.WaveformOptions(precession_version=223, mode_array=[(2, 2), (3, -3)])
    assert options != default
    changed = options.replace(mode_array=None, version_fallback=False)
    assert changed.mode_array is None and not changed.version_fallback
    assert options.mode_array is not None  # unchanged
    assert "precession_version=223" in repr(options)
    assert wf.WaveformOptions(precession_version=np.int64(102)).precession_version == 102
    modes = np.array([[2, 2], [2, -2]])
    assert wf.WaveformOptions(mode_array=modes).mode_array == [(2, 2), (2, -2)]


def test_options_object_and_keywords():
    """options= and keywords give the same waveform; keywords override."""
    a = wf.polarizations("IMRPhenomXPHM", F, precession_version=102, **PRECESSING)
    options = wf.WaveformOptions(precession_version=102)
    b = wf.polarizations("IMRPhenomXPHM", F, options=options, **PRECESSING)
    c = wf.polarizations(
        "IMRPhenomXPHM", F, options=wf.WaveformOptions(), precession_version=102, **PRECESSING
    )
    d = wf.polarizations("IMRPhenomXPHM", F, **PRECESSING)  # 320
    for x in (b, c):
        np.testing.assert_array_equal(a[0], x[0])
    assert not np.array_equal(a[0], d[0])


def test_inputs():
    """Lists, non-contiguous arrays, other dtypes and numpy scalars."""
    hp, hc = wf.polarizations("IMRPhenomXAS", F, **BINARY)
    assert hp.dtype == np.complex128 and hp.shape == F.shape
    for frequencies in (list(F), tuple(F), np.repeat(F, 2)[::2], F[::-1][::-1]):
        hp2, _ = wf.polarizations("IMRPhenomXAS", frequencies, **BINARY)
        np.testing.assert_array_equal(hp2, hp)
    # other dtypes are converted to float64
    single = F.astype(np.float32)
    np.testing.assert_array_equal(
        wf.polarizations("IMRPhenomXAS", single, **BINARY)[0],
        wf.polarizations("IMRPhenomXAS", single.astype(np.float64), **BINARY)[0],
    )
    with pytest.raises(ValueError):
        wf.polarizations("IMRPhenomXAS", np.ones((2, 3)), **BINARY)
    scalars = dict(BINARY, mass_1=np.float64(30.0), mass_2=20, spin_1=np.array([0.0, 0.0, 0.3]))
    np.testing.assert_array_equal(wf.polarizations("IMRPhenomXAS", F, **scalars)[0], hp)


def test_defaults():
    """Distance 1 Mpc, angles and f_ref zero, spins zero."""
    a = wf.polarizations("IMRPhenomD", F, mass_1=30.0, mass_2=20.0)
    b = wf.polarizations(
        "IMRPhenomD",
        F,
        mass_1=30.0,
        mass_2=20.0,
        spin_1=(0, 0, 0),
        spin_2=(0, 0, 0),
        luminosity_distance=1.0,
        inclination=0.0,
        phi_ref=0.0,
        f_ref=0.0,
    )
    np.testing.assert_array_equal(a[0], b[0])


def test_frequency_grid():
    grid = wf.FrequencyGrid.uniform(20.0, 512.0, 0.25)
    assert len(grid) == 1969
    np.testing.assert_array_equal(grid.frequencies, 20.0 + 0.25 * np.arange(1969))
    assert "1969" in repr(grid)
    custom = wf.FrequencyGrid(F)
    np.testing.assert_array_equal(custom.frequencies, F)
    a = wf.polarizations_on_grid("IMRPhenomXHM", custom, **BINARY)
    b = wf.polarizations_on_grid("IMRPhenomXHM", F, **BINARY)
    np.testing.assert_array_equal(a[0], b[0])
    c = wf.polarizations("IMRPhenomXHM", F, **BINARY)
    assert np.max(np.abs(a[0] - c[0])) <= 1e-12 * np.max(np.abs(c[0]))
    with pytest.raises(ValueError):
        wf.FrequencyGrid.uniform(20.0, 10.0, 0.25)


def test_modes_result():
    result = wf.modes("IMRPhenomXHM", F, **BINARY)
    frequencies, modes = result
    assert frequencies is not None and len(frequencies) == 2 * len(F)
    np.testing.assert_array_equal(frequencies[len(F):], F)
    np.testing.assert_array_equal(frequencies[: len(F)], -F[::-1])
    assert (2, -2) in modes and (3, 3) in modes
    np.testing.assert_array_equal(result[2, -2], modes[2, -2])
    # aligned spins: m < 0 at positive and m > 0 at negative frequencies
    assert not result[2, -2][: len(F)].any() and not result[2, 2][len(F):].any()
    assert result.theta == BINARY["inclination"] and result.phi == pytest.approx(np.pi / 2)
    with pytest.raises(KeyError):
        result[5, 5]
    # leaving modes out of the sum
    modes_22 = wf.modes("IMRPhenomXAS", F, **BINARY)
    hp, _ = wf.polarizations("IMRPhenomXAS", F, **BINARY)
    del result.modes[3, 3], result.modes[3, -3], result.modes[3, 2], result.modes[3, -2]
    del result.modes[2, 1], result.modes[2, -1], result.modes[4, 4], result.modes[4, -4]
    np.testing.assert_allclose(modes_22.polarizations()[0], hp, rtol=0, atol=1e-13 * abs(hp).max())


def test_choose_fd_modes_grid():
    result = wf.choose_fd_modes("IMRPhenomXHM", 0.25, 20.0, 512.0, **BINARY)
    n = len(result.frequencies)
    assert n % 2 == 1 and result.frequencies[n // 2] == 0.0
    assert result.frequencies[n // 2 + 1] == 0.25
    hp, hc = result.polarizations()
    assert len(hp) == n // 2 + 1


@requires_nrsur
def test_time_series_result():
    binary = dict(PRECESSING, mass_1=40.0, mass_2=20.0)
    result = wf.td_polarizations("NRSur7dq4", 1.0 / 2048, 20.0, **binary)
    hplus, hcross = result
    assert hplus is result.hplus and hplus.dtype == np.float64
    times = result.times
    assert len(times) == len(hplus)
    assert times[0] == result.epoch and result.epoch < 0 < times[-1]
    assert np.diff(times) == pytest.approx(1.0 / 2048)
    assert "TimeSeriesPolarizations" in repr(result)
    # the cached modes of the thread's workspace give the same as a fresh one
    other = wf.td_polarizations("NRSur7dq4", 1.0 / 2048, 20.0, **dict(binary, mass_1=45.0))
    again = wf.td_polarizations("NRSur7dq4", 1.0 / 2048, 20.0, **binary)
    np.testing.assert_array_equal(again.hplus, result.hplus)
    assert len(other.hplus) != len(hplus)


def test_frequency_series_result():
    result = wf.simulate_fd("IMRPhenomXAS", 0.0, 20.0, 1024.0, **BINARY)
    hplus, hcross = result
    assert hplus.dtype == np.complex128
    assert result.frequencies[1] == result.delta_f
    assert len(result.frequencies) == len(hplus)
    assert result.epoch < 0


def test_signatures_and_docs():
    signature = inspect.signature(wf.polarizations)
    assert "mass_1" in signature.parameters and "options" in signature.parameters
    assert signature.parameters["spin_1"].default == (0.0, 0.0, 0.0)
    for name in wf.__all__:
        obj = getattr(wf, name)
        if callable(obj):
            assert obj.__doc__, name


def test_results_are_not_shared():
    """Each call returns new arrays (the GIL is released while they are
    written, before they are returned)."""
    a = wf.polarizations("IMRPhenomXAS", F, **BINARY)
    b = wf.polarizations("IMRPhenomXAS", F, **BINARY)
    assert a[0] is not b[0]
    a[0][:] = 0
    assert b[0].any()


def test_pickle():
    """Options and grids can be sent to other processes."""
    options = wf.WaveformOptions(precession_version=102, mode_array=[(2, 2)])
    assert pickle.loads(pickle.dumps(options)) == options
    assert wf.WaveformOptions(**options.as_dict()) == options
    # the continuous variants survive the round trip
    continuous = wf.WaveformOptions(precession_version="SpinTaylorContinuous", nrsur_continuous=True)
    assert continuous.precession_version == "SpinTaylorContinuous"
    assert continuous.nrsur_continuous
    assert pickle.loads(pickle.dumps(continuous)) == continuous
    assert wf.WaveformOptions(**continuous.as_dict()) == continuous
    assert repr(options) == (
        "WaveformOptions(precession_version=102, mode_array=[(2, 2)], version_fallback=True, "
        "lambda_1=0.0, lambda_2=0.0, nrsur_continuous=False)"
    )
    grid = wf.FrequencyGrid.uniform(20.0, 64.0, 0.5)
    np.testing.assert_array_equal(pickle.loads(pickle.dumps(grid)).frequencies, grid.frequencies)
