// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! The crate's [`waveformrs::Error`] as Python exceptions.
//!
//! The exception classes are defined in Python (`waveformrs/exceptions.py`),
//! so that they can derive from both `WaveformError` and the matching
//! built-in exception (`ValueError`, `RuntimeError`, `OSError`, ...).

use pyo3::prelude::*;
use waveformrs::imrphenomx::IMRPhenomXError;
use waveformrs::nrsur::NRSurError;
use waveformrs::Error;

pyo3::import_exception!(waveformrs.exceptions, WaveformError);
pyo3::import_exception!(waveformrs.exceptions, InvalidParameterError);
pyo3::import_exception!(waveformrs.exceptions, UnsupportedError);
pyo3::import_exception!(waveformrs.exceptions, UnknownApproximantError);
pyo3::import_exception!(waveformrs.exceptions, ModelError);
pyo3::import_exception!(waveformrs.exceptions, PrecessionVersionUnavailableError);
pyo3::import_exception!(waveformrs.exceptions, DataFileError);
pyo3::import_exception!(waveformrs.exceptions, DataFileNotFoundError);

/// The kind of Python exception for an error of the crate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    InvalidParameter,
    Unsupported,
    UnknownApproximant,
    Model,
    PrecessionVersionUnavailable,
    DataFile,
    DataFileNotFound,
    Other,
}

fn imrphenomx_kind(err: &IMRPhenomXError) -> Kind {
    // IMRPhenomXError may gain variants: anything not listed is a parameter
    // LAL rejects.
    #[allow(unreachable_patterns)]
    match err {
        IMRPhenomXError::PrecessionVersionUnavailable { .. } => Kind::PrecessionVersionUnavailable,
        IMRPhenomXError::SpinTaylorAngles(_)
        | IMRPhenomXError::SpinTaylorFallbackToMsa(_)
        | IMRPhenomXError::PnrAngles(_) => Kind::Model,
        _ => Kind::InvalidParameter,
    }
}

fn nrsur_kind(err: &NRSurError) -> Kind {
    #[allow(unreachable_patterns)]
    match err {
        NRSurError::DataFileNotFound(_) => Kind::DataFileNotFound,
        NRSurError::Hdf5(_) | NRSurError::InvalidData(_) => Kind::DataFile,
        NRSurError::Failure(_) => Kind::Model,
        _ => Kind::InvalidParameter,
    }
}

fn kind(err: &Error) -> Kind {
    match err {
        Error::IMRPhenomX(err) => imrphenomx_kind(err),
        Error::NRSur(err) => nrsur_kind(err),
        Error::InvalidRequest(why) if *why == crate::models::PANICKED => Kind::Model,
        Error::InvalidParameter { .. } | Error::InPlaneSpins(_) | Error::InvalidRequest(_) => {
            Kind::InvalidParameter
        }
        Error::Unsupported { .. } => Kind::Unsupported,
        Error::UnknownApproximant(_) => Kind::UnknownApproximant,
        _ => Kind::Other,
    }
}

/// The Python exception for `err`, with `prefix` before its message.
pub fn to_pyerr_with_prefix(err: &Error, prefix: &str) -> PyErr {
    let message = format!("{prefix}{err}");
    match kind(err) {
        Kind::InvalidParameter => InvalidParameterError::new_err(message),
        Kind::Unsupported => UnsupportedError::new_err(message),
        Kind::UnknownApproximant => UnknownApproximantError::new_err(message),
        Kind::Model => ModelError::new_err(message),
        Kind::PrecessionVersionUnavailable => PrecessionVersionUnavailableError::new_err(message),
        Kind::DataFile => DataFileError::new_err(message),
        Kind::DataFileNotFound => DataFileNotFoundError::new_err(message),
        Kind::Other => WaveformError::new_err(message),
    }
}

/// The Python exception for `err`.
pub fn to_pyerr(err: Error) -> PyErr {
    to_pyerr_with_prefix(&err, "")
}

/// `.map_err(to_pyerr)` for results of the crate.
pub trait IntoPyResult<T> {
    fn py_result(self) -> PyResult<T>;
}

impl<T> IntoPyResult<T> for waveformrs::Result<T> {
    fn py_result(self) -> PyResult<T> {
        self.map_err(to_pyerr)
    }
}

/// An `InvalidParameterError` with `message`.
pub fn invalid_parameter(message: impl Into<String>) -> PyErr {
    InvalidParameterError::new_err(message.into())
}
