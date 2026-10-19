"""The crate's errors as Python exceptions."""

import os
import subprocess
import sys

import numpy as np
import pytest

import waveformrs as wf
from conftest import requires_nrsur

F = np.linspace(20.0, 512.0, 100)
BINARY = dict(mass_1=30.0, mass_2=20.0)


def test_hierarchy():
    for name in wf.exceptions.__all__:
        assert issubclass(getattr(wf, name), wf.WaveformError)
    assert issubclass(wf.InvalidParameterError, ValueError)
    assert issubclass(wf.UnsupportedError, ValueError)
    assert issubclass(wf.UnknownApproximantError, ValueError)
    assert issubclass(wf.ModelError, RuntimeError)
    assert issubclass(wf.PrecessionVersionUnavailableError, wf.ModelError)
    assert issubclass(wf.DataFileError, OSError)
    assert issubclass(wf.DataFileNotFoundError, FileNotFoundError)


def test_unknown_approximant():
    with pytest.raises(wf.UnknownApproximantError, match="known: TaylorF2"):
        wf.polarizations("IMRPhenomXYZ", F, **BINARY)
    with pytest.raises(ValueError):
        wf.is_precessing("imrphenomxphm")  # case-sensitive, as LAL


def test_keyword_errors():
    with pytest.raises(TypeError, match="mass_2"):
        wf.polarizations("IMRPhenomD", F, mass_1=30.0)
    with pytest.raises(TypeError, match="unexpected keyword argument 'mass1'"):
        wf.polarizations("IMRPhenomD", F, mass1=30.0, mass_2=20.0)
    with pytest.raises(TypeError, match="spin_1 must be a sequence of three"):
        wf.polarizations("IMRPhenomD", F, spin_1=(0.0, 0.1), **BINARY)
    with pytest.raises(TypeError, match="mass_1 must be a real number"):
        wf.polarizations("IMRPhenomD", F, mass_1="30", mass_2=20.0)
    with pytest.raises(TypeError, match="options must be a WaveformOptions"):
        wf.polarizations("IMRPhenomD", F, options=dict(mode_array=None), **BINARY)
    with pytest.raises(TypeError, match="unexpected keyword"):
        wf.WaveformOptions(precesion_version=223)
    with pytest.raises(TypeError, match="go in `parameters`"):
        wf.batch_polarizations("IMRPhenomD", F, dict(BINARY), inclination=0.3)


def test_invalid_parameters():
    with pytest.raises(wf.InvalidParameterError, match="masses must be positive"):
        wf.polarizations("IMRPhenomXAS", F, mass_1=-1.0, mass_2=20.0)
    with pytest.raises(wf.InvalidParameterError, match="transverse spins"):
        wf.polarizations("IMRPhenomXHM", F, spin_1=(0.1, 0.0, 0.0), **BINARY)
    with pytest.raises(wf.InvalidParameterError):
        wf.polarizations("IMRPhenomXPHM", F, spin_1=(0.9, 0.0, 0.9), **BINARY)
    with pytest.raises(wf.InvalidParameterError, match="precession version 101"):
        wf.WaveformOptions(precession_version=101)
    with pytest.raises(wf.InvalidParameterError, match="unknown precession version"):
        wf.WaveformOptions(precession_version="MSA2")
    with pytest.raises(wf.InvalidParameterError):
        wf.polarizations("IMRPhenomXHM", F, mode_array=[(5, 5)], **BINARY)
    for frequencies in (-F, np.r_[0.0, F], np.r_[F, np.nan]):
        with pytest.raises(wf.InvalidParameterError, match="positive and finite"):
            wf.polarizations("IMRPhenomXHM", frequencies, **BINARY)
        with pytest.raises(wf.InvalidParameterError, match="positive and finite"):
            wf.modes("IMRPhenomXAS", frequencies, **BINARY)
        with pytest.raises(wf.InvalidParameterError, match="positive and finite"):
            wf.FrequencyGrid(frequencies)
    with pytest.raises(wf.InvalidParameterError):
        wf.FrequencyGrid.uniform(0.0, 512.0, 0.25)


def test_unsupported():
    with pytest.raises(wf.UnsupportedError, match="mode array"):
        wf.polarizations("TaylorF2", F, mode_array=[(2, 2)], **BINARY)
    with pytest.raises(wf.UnsupportedError, match="time-domain model"):
        wf.polarizations("NRSur7dq4", F, **BINARY)
    with pytest.raises(wf.UnsupportedError, match="SimInspiralChooseFDModes"):
        wf.choose_fd_modes("IMRPhenomXAS", 0.25, 20.0, 512.0, **BINARY)
    with pytest.raises(wf.UnsupportedError, match="SimInspiralChooseTDWaveform"):
        wf.td_polarizations("TaylorF2", 1.0 / 4096, 20.0, **BINARY)
    with pytest.raises(wf.UnsupportedError, match="time-domain modes"):
        wf.td_modes("IMRPhenomXAS", 1.0 / 4096, 20.0, **BINARY)


def test_precession_version_unavailable():
    tiny = dict(BINARY, spin_1=(1e-9, 0.0, 0.3), f_ref=20.0)
    # LAL switches to the NNLO angles for in-plane spins below 1e-7 ...
    wf.polarizations("IMRPhenomXP", F, precession_version=320, **tiny)
    # ... which can be switched off
    with pytest.raises(wf.PrecessionVersionUnavailableError, match="320"):
        wf.polarizations("IMRPhenomXP", F, precession_version=320, version_fallback=False, **tiny)
    options = wf.WaveformOptions(precession_version=320, version_fallback=False)
    with pytest.raises(wf.ModelError):
        wf.polarizations("IMRPhenomXPHM", F, options=options, **tiny)


@requires_nrsur
def test_nrsur7dq4_errors():
    with pytest.raises(wf.InvalidParameterError, match="too much extrapolation"):
        wf.td_polarizations("NRSur7dq4", 1.0 / 4096, 20.0, mass_1=70.0, mass_2=5.0)
    with pytest.raises(wf.InvalidParameterError, match="sampling interval"):
        wf.td_polarizations("NRSur7dq4", -1.0, 20.0, mass_1=40.0, mass_2=20.0)


def test_conversion_arguments():
    with pytest.raises(wf.InvalidParameterError, match="positive"):
        wf.simulate_td("IMRPhenomXAS", 0.0, 20.0, **BINARY)


def test_data_file_not_found():
    """Without NRSUR7DQ4_DATA and LAL_DATA_PATH NRSur7dq4 raises
    DataFileNotFoundError (in a fresh process: the data are loaded once)."""
    env = {k: v for k, v in os.environ.items() if k not in ("NRSUR7DQ4_DATA", "LAL_DATA_PATH")}
    code = (
        "import waveformrs as wf\n"
        "assert wf.nrsur7dq4_data_file() is None\n"
        "try:\n"
        "    wf.td_polarizations('NRSur7dq4', 1 / 4096, 20.0, mass_1=40.0, mass_2=20.0)\n"
        "except FileNotFoundError as error:\n"
        "    print(type(error).__name__)\n"
    )
    output = subprocess.run(
        [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
    )
    assert output.stdout.strip() == "DataFileNotFoundError"
