// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Conversions between array-likes and contiguous `f64` slices.

use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;

use crate::errors::invalid_parameter;

/// `numpy.ascontiguousarray(value, dtype=numpy.float64)`.
fn ascontiguous_f64<'py>(value: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    static ASCONTIGUOUSARRAY: PyOnceLock<Py<PyAny>> = PyOnceLock::new();
    let py = value.py();
    ASCONTIGUOUSARRAY
        .import(py, "numpy", "ascontiguousarray")?
        .call1((value, "float64"))
}

/// `value` as a contiguous one-dimensional `float64` array, without a copy
/// when it is one already.
pub fn f64_array1<'py>(
    value: &Bound<'py, PyAny>,
    name: &str,
) -> PyResult<PyReadonlyArray1<'py, f64>> {
    if let Ok(array) = value.extract::<PyReadonlyArray1<'py, f64>>() {
        if array.is_c_contiguous() {
            return Ok(array);
        }
    }
    let converted = ascontiguous_f64(value)?;
    converted
        .extract::<PyReadonlyArray1<'py, f64>>()
        .map_err(|_| {
            PyValueError::new_err(format!(
                "{name} must be a one-dimensional array of real numbers"
            ))
        })
}

/// `value` as a C-contiguous two-dimensional `float64` array.
pub fn f64_array2<'py>(
    value: &Bound<'py, PyAny>,
    name: &str,
) -> PyResult<PyReadonlyArray2<'py, f64>> {
    if let Ok(array) = value.extract::<PyReadonlyArray2<'py, f64>>() {
        if array.is_c_contiguous() {
            return Ok(array);
        }
    }
    let converted = ascontiguous_f64(value)?;
    converted
        .extract::<PyReadonlyArray2<'py, f64>>()
        .map_err(|_| {
            PyValueError::new_err(format!(
                "{name} must be a two-dimensional array of real numbers"
            ))
        })
}

/// Frequencies must be positive and finite (the models are not defined
/// elsewhere, and a zero or negative first frequency would become the
/// reference frequency of some).
pub fn check_frequencies(frequencies: &[f64]) -> PyResult<()> {
    match frequencies.iter().find(|&&f| !(f > 0.0 && f.is_finite())) {
        Some(f) => Err(invalid_parameter(format!(
            "frequencies must be positive and finite, got {f}"
        ))),
        None => Ok(()),
    }
}
