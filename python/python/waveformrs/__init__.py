# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
# from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
"""Gravitational-wave models that reproduce LALSimulation 6.2.0 to rounding.

Python bindings of the Rust crate ``waveformrs``: TaylorF2, IMRPhenomD,
IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, IMRPhenomXPNR (frequency
domain) and NRSur7dq4 (time domain), named as in LALSimulation.

The binary is given by keyword arguments with LAL's parameters in the units of
``SimInspiralChooseFDWaveform`` except for the masses and the distance:

* ``mass_1``, ``mass_2``: component masses in solar masses (required);
* ``spin_1``, ``spin_2``: dimensionless spin vectors ``(x, y, z)`` at the
  reference frequency, with the orbital angular momentum along ``z``
  (default zero);
* ``luminosity_distance``: in Mpc (default 1);
* ``inclination``, ``phi_ref``: in rad (default zero);
* ``f_ref``: reference frequency in Hz (default zero: LAL's default of each
  model).

The settings of LAL's ``LALDict`` are :class:`WaveformOptions`, passed as
``options=`` and/or as keyword arguments (``precession_version``,
``mode_array``, ``version_fallback``, ``pnr_parity_fix``).

>>> import numpy as np
>>> import waveformrs
>>> frequencies = np.linspace(20.0, 1024.0, 4017)
>>> hplus, hcross = waveformrs.polarizations(
...     "IMRPhenomXPHM", frequencies, mass_1=36.0, mass_2=29.0,
...     spin_1=(0.4, 0.2, 0.3), spin_2=(-0.2, 0.3, -0.1),
...     luminosity_distance=500.0, inclination=0.7, phi_ref=1.2, f_ref=20.0,
...     precession_version=223,
... )

Every function releases the GIL while it computes; the batch functions
(:func:`batch_polarizations`, :func:`batch_td_polarizations`) also spread
many binaries over threads.
"""

from . import exceptions
from ._batch import batch_polarizations, batch_td_polarizations, parameter_table
from ._waveformrs import (
    BatchPolarizations,
    BatchTimeSeries,
    FrequencyGrid,
    FrequencyModes,
    FrequencySeriesPolarizations,
    TimeSeriesModes,
    TimeSeriesPolarizations,
    WaveformOptions,
    __version__,
    approximants,
    choose_fd_modes,
    is_frequency_domain,
    is_precessing,
    modes,
    nrsur7dq4_data_file,
    polarizations,
    polarizations_on_grid,
    simulate_fd,
    simulate_td,
    td_modes,
    td_polarizations,
)
from .exceptions import (
    DataFileError,
    DataFileNotFoundError,
    InvalidParameterError,
    ModelError,
    PrecessionVersionUnavailableError,
    UnknownApproximantError,
    UnsupportedError,
    WaveformError,
)

__all__ = [
    "__version__",
    # approximants and options
    "approximants",
    "is_frequency_domain",
    "is_precessing",
    "WaveformOptions",
    # frequency domain
    "polarizations",
    "FrequencyGrid",
    "polarizations_on_grid",
    "modes",
    "choose_fd_modes",
    "FrequencyModes",
    # time domain and conversions
    "td_polarizations",
    "td_modes",
    "simulate_td",
    "simulate_fd",
    "TimeSeriesPolarizations",
    "TimeSeriesModes",
    "FrequencySeriesPolarizations",
    "nrsur7dq4_data_file",
    # batches
    "batch_polarizations",
    "batch_td_polarizations",
    "parameter_table",
    "BatchPolarizations",
    "BatchTimeSeries",
    # errors
    "exceptions",
    "WaveformError",
    "InvalidParameterError",
    "UnsupportedError",
    "UnknownApproximantError",
    "ModelError",
    "PrecessionVersionUnavailableError",
    "DataFileError",
    "DataFileNotFoundError",
]
