# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
# from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
"""Batched evaluation: tables of binaries from dicts of arrays or structured
arrays, passed to the batch functions of the extension module."""

from __future__ import annotations

from typing import Any, Mapping, Optional, Union

import numpy as np

from . import _waveformrs
from ._waveformrs import (
    BatchPolarizations,
    BatchTimeSeries,
    FrequencyGrid,
    WaveformOptions,
)

__all__ = ["batch_polarizations", "batch_td_polarizations", "parameter_table"]

#: The columns of a table of binaries, in order.
COLUMNS = (
    "mass_1",
    "mass_2",
    "spin_1x",
    "spin_1y",
    "spin_1z",
    "spin_2x",
    "spin_2y",
    "spin_2z",
    "luminosity_distance",
    "inclination",
    "phi_ref",
    "f_ref",
)

_DEFAULTS = {
    "spin_1x": 0.0,
    "spin_1y": 0.0,
    "spin_1z": 0.0,
    "spin_2x": 0.0,
    "spin_2y": 0.0,
    "spin_2z": 0.0,
    "luminosity_distance": 1.0,
    "inclination": 0.0,
    "phi_ref": 0.0,
    "f_ref": 0.0,
}

_SPIN_VECTORS = ("spin_1", "spin_2")
_ALLOWED = frozenset(COLUMNS) | frozenset(_SPIN_VECTORS)

Parameters = Union[Mapping[str, Any], np.ndarray]


def parameter_table(parameters: Parameters) -> np.ndarray:
    """The binaries of ``parameters`` as a C-contiguous ``(n, 12)`` float64
    table with the columns ``COLUMNS``.

    ``parameters`` is a mapping (a ``dict``, a pandas ``DataFrame``, ...) or
    a numpy structured array with the names of the keyword arguments of
    :func:`waveformrs.polarizations`: ``mass_1`` and ``mass_2`` (required),
    ``spin_1`` and ``spin_2`` as ``(n, 3)`` arrays or their components
    ``spin_1x``, ``spin_1y``, ``spin_1z``, ... (default zero),
    ``luminosity_distance`` (default 1 Mpc), ``inclination``, ``phi_ref`` and
    ``f_ref`` (default zero). Scalars are broadcast to the other columns.
    An ``(n, 12)`` array is taken as the table itself.
    """
    if isinstance(parameters, np.ndarray) and parameters.dtype.names is None:
        table = np.ascontiguousarray(parameters, dtype=np.float64)
        if table.ndim != 2 or table.shape[1] != len(COLUMNS):
            raise ValueError(
                f"an unstructured parameters array must have shape (n, {len(COLUMNS)}) "
                f"with the columns {', '.join(COLUMNS)}"
            )
        return table
    if isinstance(parameters, np.ndarray):
        names = list(parameters.dtype.names)
    elif hasattr(parameters, "keys") and hasattr(parameters, "__getitem__"):
        names = list(parameters.keys())
    else:
        raise TypeError(
            "parameters must be a mapping of arrays (e.g. a dict), a structured array "
            f"or an (n, {len(COLUMNS)}) array, not {type(parameters).__name__}"
        )
    unknown = [name for name in names if name not in _ALLOWED]
    if unknown:
        raise TypeError(
            f"unknown binary parameters {unknown}: use {', '.join(sorted(_ALLOWED))}"
        )
    for name in ("mass_1", "mass_2"):
        if name not in names:
            raise TypeError(f"missing binary parameter {name!r}")

    values = {}
    for name in names:
        if name in _SPIN_VECTORS:
            components = [f"{name}{axis}" for axis in "xyz"]
            if any(component in names for component in components):
                raise TypeError(f"give {name} or its components, not both")
            spin = np.asarray(parameters[name], dtype=np.float64)
            if spin.shape[-1:] != (3,):
                raise ValueError(f"{name} must have shape (3,) or (n, 3), not {spin.shape}")
            for axis, component in enumerate(components):
                values[component] = spin[..., axis]
        else:
            values[name] = np.asarray(parameters[name], dtype=np.float64)
    columns = [values.get(name, _DEFAULTS.get(name)) for name in COLUMNS]
    columns = np.broadcast_arrays(*[np.asarray(column, dtype=np.float64) for column in columns])
    shape = columns[0].shape
    if len(shape) > 1:
        raise ValueError(f"the binary parameters must be one-dimensional, not of shape {shape}")
    table = np.empty((shape[0] if shape else 1, len(COLUMNS)), dtype=np.float64)
    for index, column in enumerate(columns):
        table[:, index] = column
    return table


def batch_polarizations(
    approximant: str,
    grid: Union[FrequencyGrid, Any],
    parameters: Parameters,
    *,
    threads: int = 1,
    strict: bool = False,
    exact: bool = False,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> BatchPolarizations:
    """The frequency-domain polarizations of many binaries on one grid.

    Parameters
    ----------
    approximant:
        A frequency-domain approximant (LAL's name).
    grid:
        A :class:`FrequencyGrid`, or the frequencies (Hz) to build one from.
    parameters:
        The binaries: a dict of arrays (or scalars), a structured array or an
        ``(n, 12)`` table (see :func:`parameter_table`).
    threads:
        The number of threads (0 for all cores). The results do not depend on
        it, and the GIL is released while they are computed.
    strict:
        Raise the error of the first binary that fails instead of reporting
        the failures in ``ok`` and ``errors``.
    exact:
        Evaluate each binary as :func:`polarizations` does
        (``SimInspiralChooseFDWaveformSequence``, bit for bit) instead of with
        the fast paths of :func:`polarizations_on_grid` (which agree to about
        1e-13 relative).
    options, **option_kwargs:
        The :class:`WaveformOptions`, as for :func:`polarizations`.

    Returns
    -------
    BatchPolarizations
        ``hplus`` and ``hcross`` as ``(n_binaries, n_frequencies)`` complex128
        arrays (zero rows for failed binaries), the mask ``ok`` and the
        exceptions ``errors`` (``{index: exception}``). Unpacks as
        ``hplus, hcross = result``.
    """
    return _waveformrs.batch_polarizations(
        approximant,
        grid,
        parameter_table(parameters),
        threads=threads,
        strict=strict,
        exact=exact,
        options=options,
        **option_kwargs,
    )


def batch_td_polarizations(
    approximant: str,
    delta_t: float,
    f_min: float,
    parameters: Parameters,
    *,
    threads: int = 1,
    strict: bool = False,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> BatchTimeSeries:
    """The time-domain polarizations (``SimInspiralChooseTDWaveform``, as
    :func:`td_polarizations`) of many binaries, for NRSur7dq4 with one
    workspace per thread.

    Parameters are as for :func:`batch_polarizations`, with the sampling
    interval ``delta_t`` (s) and the starting frequency ``f_min`` (Hz).

    Returns
    -------
    BatchTimeSeries
        ``waveforms[i]``, a :class:`TimeSeriesPolarizations` or ``None`` for
        the binaries that failed, with ``ok`` and ``errors`` as for
        :func:`batch_polarizations`.
    """
    return _waveformrs.batch_td_polarizations(
        approximant,
        delta_t,
        f_min,
        parameter_table(parameters),
        threads=threads,
        strict=strict,
        options=options,
        **option_kwargs,
    )
