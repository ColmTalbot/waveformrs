# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
# from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
"""Type stubs of the extension module ``waveformrs._waveformrs``.

The binary is given by keyword arguments (``BinaryKwargs``) followed by the
options of :class:`WaveformOptions` (``options=`` and/or its keywords).
"""

from typing import Any, Dict, Iterator, List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt

__version__: str

ComplexArray = npt.NDArray[np.complex128]
RealArray = npt.NDArray[np.float64]
Mode = Tuple[int, int]
Spin = Union[Sequence[float], npt.ArrayLike]
PrecessionVersion = Union[Literal[102, 223, 320], Literal["NNLO", "MSA", "SpinTaylor"], int, str]

class WaveformOptions:
    """The settings of LALSimulation's ``LALDict`` that the models take.

    * ``precession_version``: LAL's ``PhenomXPrecVersion`` of IMRPhenomXP/XPHM:
      102 (``"NNLO"``), 223 (``"MSA"``) or 320 (``"SpinTaylor"``, the default
      here; LALSimulation 6.2.0 uses 223, see :meth:`lalsimulation_defaults`).
    * ``mode_array``: LAL's ``ModeArray`` (``(l, m)`` pairs) or ``None`` for
      all modes.
    * ``version_fallback``: whether IMRPhenomXP/XPHM/XPNR may change the
      precession version as LALSimulation does (default ``True``).
    * ``lambda_1``, ``lambda_2``: LAL's ``lambda1``, ``lambda2``, the tidal
      deformabilities of the bodies of ``mass_1`` and ``mass_2`` (default 0)
      of the NRTidal approximants; the others reject non-zero values.
    """

    def __init__(
        self,
        *,
        precession_version: PrecessionVersion = 320,
        mode_array: Optional[Sequence[Sequence[int]]] = None,
        version_fallback: bool = True,
        lambda_1: float = 0.0,
        lambda_2: float = 0.0,
        **options: Any,
    ) -> None: ...
    @staticmethod
    def lalsimulation_defaults() -> WaveformOptions:
        """LALSimulation 6.2.0's defaults: MSA angles (223) and all modes."""
    def replace(self, **options: Any) -> WaveformOptions:
        """A copy with the given options changed."""
    @property
    def precession_version(self) -> int: ...
    @property
    def mode_array(self) -> Optional[List[Mode]]: ...
    @property
    def version_fallback(self) -> bool: ...
    @property
    def lambda_1(self) -> float: ...
    @property
    def lambda_2(self) -> float: ...
    def __eq__(self, other: object) -> bool: ...

class FrequencyGrid:
    """Frequencies (Hz) with the powers the models need precomputed, for
    many waveforms on the same frequencies."""

    def __init__(self, frequencies: npt.ArrayLike) -> None: ...
    @staticmethod
    def uniform(f_min: float, f_max: float, delta_f: float) -> FrequencyGrid:
        """``f_min, f_min + delta_f, ...`` up to and including ``f_max``."""
    @property
    def frequencies(self) -> RealArray: ...
    def __len__(self) -> int: ...

class TimeSeriesPolarizations:
    """``h_+(t)``, ``h_x(t)`` sampled every ``delta_t`` s from ``epoch`` (s).
    Unpacks as ``hplus, hcross``."""

    @property
    def hplus(self) -> RealArray: ...
    @property
    def hcross(self) -> RealArray: ...
    @property
    def epoch(self) -> float: ...
    @property
    def delta_t(self) -> float: ...
    @property
    def times(self) -> RealArray: ...
    def __iter__(self) -> Iterator[RealArray]: ...

class FrequencySeriesPolarizations:
    """``h_+(f)``, ``h_x(f)`` on ``k * delta_f`` with the ``epoch`` (s) of the
    time series. Unpacks as ``hplus, hcross``."""

    @property
    def hplus(self) -> ComplexArray: ...
    @property
    def hcross(self) -> ComplexArray: ...
    @property
    def epoch(self) -> float: ...
    @property
    def delta_f(self) -> float: ...
    @property
    def frequencies(self) -> RealArray: ...
    def __iter__(self) -> Iterator[ComplexArray]: ...

class FrequencyModes:
    """Modes ``h_lm(f)`` on two-sided frequencies (LAL's
    ``SphHarmFrequencySeries`` layout). Unpacks as ``frequencies, modes``."""

    @property
    def frequencies(self) -> RealArray: ...
    @property
    def modes(self) -> Dict[Mode, ComplexArray]: ...
    @property
    def theta(self) -> float: ...
    @property
    def phi(self) -> float: ...
    @property
    def polarization_angle(self) -> float: ...
    def polarizations(self) -> Tuple[ComplexArray, ComplexArray]:
        """``(h_+, h_x)`` at the non-negative frequencies from ``modes``."""
    def __getitem__(self, mode: Mode) -> ComplexArray: ...
    def __iter__(self) -> Iterator[Union[RealArray, Dict[Mode, ComplexArray]]]: ...

class TimeSeriesModes:
    """Modes ``h_lm(t)`` sampled every ``delta_t`` s from ``epoch`` (s).
    Unpacks as ``times, modes``."""

    @property
    def modes(self) -> Dict[Mode, ComplexArray]: ...
    @property
    def epoch(self) -> float: ...
    @property
    def delta_t(self) -> float: ...
    @property
    def times(self) -> RealArray: ...
    def __getitem__(self, mode: Mode) -> ComplexArray: ...
    def __iter__(self) -> Iterator[Union[RealArray, Dict[Mode, ComplexArray]]]: ...

class BatchPolarizations:
    """``(n_binaries, n_frequencies)`` polarizations with the success mask
    ``ok`` and ``errors`` (``{index: exception}``). Unpacks as ``hplus,
    hcross``."""

    @property
    def hplus(self) -> ComplexArray: ...
    @property
    def hcross(self) -> ComplexArray: ...
    @property
    def ok(self) -> npt.NDArray[np.bool_]: ...
    @property
    def errors(self) -> Dict[int, Exception]: ...
    def __iter__(self) -> Iterator[ComplexArray]: ...

class BatchTimeSeries:
    """One :class:`TimeSeriesPolarizations` (or ``None``) per binary, with
    ``ok`` and ``errors``."""

    @property
    def waveforms(self) -> List[Optional[TimeSeriesPolarizations]]: ...
    @property
    def ok(self) -> npt.NDArray[np.bool_]: ...
    @property
    def errors(self) -> Dict[int, Exception]: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Optional[TimeSeriesPolarizations]: ...
    def __iter__(self) -> Iterator[Optional[TimeSeriesPolarizations]]: ...

def approximants(domain: Optional[Literal["frequency", "time"]] = None) -> List[str]:
    """LAL's names of the approximants (of one domain)."""

def is_precessing(approximant: str) -> bool: ...
def is_frequency_domain(approximant: str) -> bool: ...
def nrsur7dq4_data_file() -> Optional[str]:
    """The NRSur7dq4 data file from ``NRSUR7DQ4_DATA`` or ``LAL_DATA_PATH``."""

def polarizations(
    approximant: str,
    frequencies: npt.ArrayLike,
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> Tuple[ComplexArray, ComplexArray]:
    """``(h_+, h_x)`` at ``frequencies`` (Hz), as
    ``SimInspiralChooseFDWaveformSequence``."""

def polarizations_on_grid(
    approximant: str,
    grid: Union[FrequencyGrid, npt.ArrayLike],
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> Tuple[ComplexArray, ComplexArray]:
    """``(h_+, h_x)`` on a :class:`FrequencyGrid` with the fast grid paths."""

def modes(
    approximant: str,
    frequencies: npt.ArrayLike,
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> FrequencyModes:
    """The modes ``h_lm`` at ``+-frequencies`` (``SimInspiralChooseFDModes``
    conventions)."""

def choose_fd_modes(
    approximant: str,
    delta_f: float,
    f_min: float,
    f_max: float = 0.0,
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> FrequencyModes:
    """``SimInspiralChooseFDModes`` (IMRPhenomXHM, IMRPhenomXPHM,
    IMRPhenomXPNR)."""

def td_polarizations(
    approximant: str,
    delta_t: float,
    f_min: float,
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> TimeSeriesPolarizations:
    """``SimInspiralChooseTDWaveform``."""

def td_modes(
    approximant: str,
    delta_t: float,
    f_min: float,
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> TimeSeriesModes:
    """``SimInspiralChooseTDModes`` (NRSur7dq4)."""

def simulate_td(
    approximant: str,
    delta_t: float,
    f_min: float,
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> TimeSeriesPolarizations:
    """``SimInspiralTD``: the conditioned time-domain waveform."""

def simulate_fd(
    approximant: str,
    delta_f: float,
    f_min: float,
    f_max: float,
    *,
    mass_1: float,
    mass_2: float,
    spin_1: Spin = (0.0, 0.0, 0.0),
    spin_2: Spin = (0.0, 0.0, 0.0),
    luminosity_distance: float = 1.0,
    inclination: float = 0.0,
    phi_ref: float = 0.0,
    f_ref: float = 0.0,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> FrequencySeriesPolarizations:
    """``SimInspiralFD``: the conditioned frequency-domain waveform."""

def batch_polarizations(
    approximant: str,
    grid: Union[FrequencyGrid, npt.ArrayLike],
    parameters: npt.ArrayLike,
    *,
    threads: int = 1,
    strict: bool = False,
    exact: bool = False,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> BatchPolarizations:
    """Polarizations of the ``(n, 12)`` table of binaries ``parameters``
    (see ``waveformrs.parameter_table``) on one grid."""

def batch_td_polarizations(
    approximant: str,
    delta_t: float,
    f_min: float,
    parameters: npt.ArrayLike,
    *,
    threads: int = 1,
    strict: bool = False,
    options: Optional[WaveformOptions] = None,
    **option_kwargs: Any,
) -> BatchTimeSeries:
    """Time-domain polarizations of the ``(n, 12)`` table of binaries."""
