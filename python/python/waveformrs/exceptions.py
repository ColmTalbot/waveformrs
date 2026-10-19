"""Exceptions raised by waveformrs.

Every error of the Rust crate (``waveformrs::Error``) is raised as a subclass
of :class:`WaveformError`, which also derives from the matching built-in
exception, so ``except ValueError`` catches invalid parameters as it would for
most Python code.

=====================================  ==========================================
exception                              raised for
=====================================  ==========================================
:class:`InvalidParameterError`         parameters LALSimulation rejects (masses,
                                       spins, frequencies, in-plane spins for an
                                       aligned-spin model, mode arrays, ...)
:class:`UnsupportedError`              requests the approximant does not support
                                       (a mode array for a (2,2)-mode model, a
                                       time-domain request for a frequency-domain
                                       model, ...)
:class:`UnknownApproximantError`       names that are not approximants of the
                                       crate
:class:`ModelError`                    failures of the model itself (SpinTaylor
                                       or PNR angles, NRSur7dq4 evaluation)
:class:`PrecessionVersionUnavailableError`  the requested Euler angles cannot be
                                       used and ``version_fallback=False``
:class:`DataFileError`                 the NRSur7dq4 data file cannot be read
:class:`DataFileNotFoundError`         the NRSur7dq4 data file was not found
=====================================  ==========================================
"""

__all__ = [
    "WaveformError",
    "InvalidParameterError",
    "UnsupportedError",
    "UnknownApproximantError",
    "ModelError",
    "PrecessionVersionUnavailableError",
    "DataFileError",
    "DataFileNotFoundError",
]


class WaveformError(Exception):
    """Base class of the errors of waveformrs."""


class InvalidParameterError(WaveformError, ValueError):
    """A parameter that LALSimulation rejects for this approximant."""


class UnsupportedError(WaveformError, ValueError):
    """The approximant does not support the request."""


class UnknownApproximantError(WaveformError, ValueError):
    """The name is not one of :func:`waveformrs.approximants`."""


class ModelError(WaveformError, RuntimeError):
    """The model could not be evaluated for this binary."""


class PrecessionVersionUnavailableError(ModelError):
    """The requested precession version (``PhenomXPrecVersion``) cannot be
    used for this binary and the automatic fallback that LALSimulation would
    make is switched off (``version_fallback=False``)."""


class DataFileError(WaveformError, OSError):
    """The NRSur7dq4 data file could not be read."""


class DataFileNotFoundError(DataFileError, FileNotFoundError):
    """The NRSur7dq4 data file was not found: set ``NRSUR7DQ4_DATA`` to the
    file or ``LAL_DATA_PATH`` to its directory."""
