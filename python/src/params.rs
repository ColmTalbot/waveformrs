// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Approximant names, `WaveformOptions` and the keyword arguments of a
//! binary.
//!
//! Every waveform function takes the parameters of the binary as keyword
//! arguments with the names of [`waveformrs::Binary`] and, in the same
//! `**kwargs`, the options of [`waveformrs::WaveformOptions`] (optionally on
//! top of a `WaveformOptions` object passed as `options=`).
//!
//! To bind a new option of `WaveformOptions`, add an arm to
//! [`apply_option`] and an entry to [`option_items`] (and its name to
//! [`OPTION_NAMES`], a getter to [`PyWaveformOptions`] and a line to the
//! stub `_waveformrs.pyi`).

use std::str::FromStr;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyString, PyTuple};
use waveformrs::imrphenomx::PrecessionVersion;
use waveformrs::{Approximant, Binary, WaveformOptions};

use crate::errors::{invalid_parameter, IntoPyResult};

/// The approximant named `name` (LAL's name, case-sensitive).
pub fn parse_approximant(name: &str) -> PyResult<Approximant> {
    Approximant::from_str(name).py_result()
}

/// The keyword arguments of a binary, in the order of `Binary`.
pub const BINARY_NAMES: [&str; 8] = [
    "mass_1",
    "mass_2",
    "spin_1",
    "spin_2",
    "luminosity_distance",
    "inclination",
    "phi_ref",
    "f_ref",
];

/// The options that can be given as keyword arguments.
pub const OPTION_NAMES: [&str; 5] = [
    "precession_version",
    "mode_array",
    "version_fallback",
    "lambda_1",
    "lambda_2",
];

/// Set the option `name` of `options` to `value`; `Ok(false)` if `name` is
/// not an option.
///
/// Each arm maps a keyword argument onto a builder method of
/// `WaveformOptions` (or a field, to reset it).
fn apply_option(
    options: &mut WaveformOptions,
    name: &str,
    value: &Bound<'_, PyAny>,
) -> PyResult<bool> {
    let current = std::mem::take(options);
    *options = match name {
        "precession_version" => current.with_precession_version(precession_version(value)?),
        "mode_array" => {
            if value.is_none() {
                let mut options = current;
                options.mode_array = None;
                options
            } else {
                current.with_mode_array(&mode_array(value)?)
            }
        }
        "version_fallback" => current.with_version_fallback(extract_bool(value, name)?),
        "lambda_1" => {
            let lambda_2 = current.lambda_2;
            current.with_tidal_deformabilities(extract_f64(value, name)?, lambda_2)
        }
        "lambda_2" => {
            let lambda_1 = current.lambda_1;
            current.with_tidal_deformabilities(lambda_1, extract_f64(value, name)?)
        }
        _ => {
            *options = current;
            return Ok(false);
        }
    };
    Ok(true)
}

/// The options as `(keyword, value)`, for `as_dict`, `repr` and pickling:
/// `WaveformOptions(**dict(option_items(options)))` is `options`.
fn option_items<'py>(
    py: Python<'py>,
    options: &WaveformOptions,
) -> PyResult<Vec<(&'static str, Bound<'py, PyAny>)>> {
    Ok(vec![
        (
            "precession_version",
            options
                .precession_version
                .lal_flag()
                .into_pyobject(py)?
                .into_any(),
        ),
        (
            "mode_array",
            options.mode_array.clone().into_pyobject(py)?.into_any(),
        ),
        (
            "version_fallback",
            options
                .version_fallback
                .into_pyobject(py)?
                .to_owned()
                .into_any(),
        ),
        ("lambda_1", options.lambda_1.into_pyobject(py)?.into_any()),
        ("lambda_2", options.lambda_2.into_pyobject(py)?.into_any()),
    ])
}

/// A precession version from LAL's `PhenomXPrecVersion` (102, 223, 320) or
/// a name (`"NNLO"`, `"MSA"`, `"SpinTaylor"`, `"SpinTaylorContinuous"`,
/// case-insensitive).
fn precession_version(value: &Bound<'_, PyAny>) -> PyResult<PrecessionVersion> {
    if let Ok(name) = value.cast::<PyString>() {
        let name = name.to_cow()?;
        return match name.to_ascii_lowercase().as_str() {
            "nnlo" => Ok(PrecessionVersion::Nnlo),
            "msa" => Ok(PrecessionVersion::Msa),
            "spintaylor" | "spin_taylor" => Ok(PrecessionVersion::SpinTaylor),
            "spintaylorcontinuous" | "spintaylor_continuous" | "spin_taylor_continuous" => {
                Ok(PrecessionVersion::SpinTaylorContinuous)
            }
            _ => Err(invalid_parameter(format!(
                "unknown precession version {name:?}: use 102 (\"NNLO\"), 223 (\"MSA\"), \
                 320 (\"SpinTaylor\") or \"SpinTaylorContinuous\""
            ))),
        };
    }
    let flag: i64 = value.extract().map_err(|_| {
        PyTypeError::new_err(
            "precession_version must be an int (LAL's PhenomXPrecVersion) or a str",
        )
    })?;
    match flag {
        102 => Ok(PrecessionVersion::Nnlo),
        223 => Ok(PrecessionVersion::Msa),
        320 => Ok(PrecessionVersion::SpinTaylor),
        _ => Err(invalid_parameter(format!(
            "precession version {flag} is not implemented: use 102 (NNLO), 223 (MSA) or 320 \
             (SpinTaylor)"
        ))),
    }
}

/// A mode array: a sequence of `(l, m)` pairs.
fn mode_array(value: &Bound<'_, PyAny>) -> PyResult<Vec<(u32, i32)>> {
    let error = || PyTypeError::new_err("mode_array must be a sequence of (l, m) pairs or None");
    if value.is_instance_of::<PyString>() {
        return Err(error());
    }
    let mut modes = Vec::new();
    for item in value.try_iter().map_err(|_| error())? {
        let [ell, emm]: [i64; 2] = item?.extract().map_err(|_| error())?;
        let ell =
            u32::try_from(ell).map_err(|_| invalid_parameter(format!("invalid l = {ell}")))?;
        let emm =
            i32::try_from(emm).map_err(|_| invalid_parameter(format!("invalid m = {emm}")))?;
        modes.push((ell, emm));
    }
    Ok(modes)
}

fn extract_bool(value: &Bound<'_, PyAny>, name: &str) -> PyResult<bool> {
    value
        .extract()
        .map_err(|_| PyTypeError::new_err(format!("{name} must be a bool")))
}

fn extract_f64(value: &Bound<'_, PyAny>, name: &str) -> PyResult<f64> {
    value.extract().map_err(|_| {
        PyTypeError::new_err(format!(
            "{name} must be a real number, got {}",
            type_name(value)
        ))
    })
}

fn extract_spin(value: &Bound<'_, PyAny>, name: &str) -> PyResult<[f64; 3]> {
    value.extract().map_err(|_| {
        PyTypeError::new_err(format!(
            "{name} must be a sequence of three real numbers (x, y, z), got {}",
            type_name(value)
        ))
    })
}

fn type_name(value: &Bound<'_, PyAny>) -> String {
    value
        .get_type()
        .name()
        .map(|name| name.to_string())
        .unwrap_or_else(|_| "?".to_string())
}

fn unknown_keyword(name: &str) -> PyErr {
    PyTypeError::new_err(format!(
        "unexpected keyword argument '{name}': the binary parameters are {} and the options \
         {} (or options=WaveformOptions(...))",
        BINARY_NAMES.join(", "),
        OPTION_NAMES.join(", ")
    ))
}

/// The options of a `WaveformOptions` object (if any) changed by the
/// option keywords.
fn options_with<'py>(
    base: Option<&Bound<'py, PyAny>>,
    keywords: impl IntoIterator<Item = (Bound<'py, PyString>, Bound<'py, PyAny>)>,
) -> PyResult<WaveformOptions> {
    let mut options = match base {
        Some(base) if !base.is_none() => base
            .cast::<PyWaveformOptions>()
            .map_err(|_| PyTypeError::new_err("options must be a WaveformOptions or None"))?
            .get()
            .inner
            .clone(),
        _ => WaveformOptions::default(),
    };
    for (key, value) in keywords {
        let name = key.to_cow()?;
        if !apply_option(&mut options, &name, &value)? {
            return Err(unknown_keyword(&name));
        }
    }
    Ok(options)
}

/// The options given as keyword arguments (`options=` and the option
/// names), rejecting anything else.
pub fn parse_options(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<WaveformOptions> {
    let mut base = None;
    let mut keywords = Vec::new();
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs.iter() {
            let key = key.cast_into::<PyString>()?;
            let name = key.to_cow()?;
            if name == "options" {
                base = Some(value);
            } else if BINARY_NAMES.contains(&name.as_ref()) {
                return Err(PyTypeError::new_err(format!(
                    "unexpected keyword argument '{name}': the binary parameters of a batch \
                     go in `parameters`"
                )));
            } else {
                keywords.push((key, value));
            }
        }
    }
    options_with(base.as_ref(), keywords)
}

/// The binary and the options given as keyword arguments.
///
/// `mass_1` and `mass_2` are required; the spins default to zero, the
/// distance to 1 Mpc, the angles to zero and `f_ref` to zero (LAL's default
/// reference frequency of each model).
pub fn parse_request(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<(Binary, WaveformOptions)> {
    let mut mass_1 = None;
    let mut mass_2 = None;
    let mut binary = Binary {
        mass_1: f64::NAN,
        mass_2: f64::NAN,
        spin_1: [0.0; 3],
        spin_2: [0.0; 3],
        luminosity_distance: 1.0,
        inclination: 0.0,
        phi_ref: 0.0,
        f_ref: 0.0,
    };
    let mut base = None;
    let mut keywords = Vec::new();
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs.iter() {
            let key = key.cast_into::<PyString>()?;
            let name = key.to_cow()?;
            let name = name.as_ref();
            match name {
                "mass_1" => mass_1 = Some(extract_f64(&value, name)?),
                "mass_2" => mass_2 = Some(extract_f64(&value, name)?),
                "spin_1" => binary.spin_1 = extract_spin(&value, name)?,
                "spin_2" => binary.spin_2 = extract_spin(&value, name)?,
                "luminosity_distance" => binary.luminosity_distance = extract_f64(&value, name)?,
                "inclination" => binary.inclination = extract_f64(&value, name)?,
                "phi_ref" => binary.phi_ref = extract_f64(&value, name)?,
                "f_ref" => binary.f_ref = extract_f64(&value, name)?,
                "options" => base = Some(value),
                _ => keywords.push((key, value)),
            }
        }
    }
    // unknown keywords (typos) first
    let options = options_with(base.as_ref(), keywords)?;
    let missing =
        |name: &str| PyTypeError::new_err(format!("missing required keyword argument '{name}'"));
    binary.mass_1 = mass_1.ok_or_else(|| missing("mass_1"))?;
    binary.mass_2 = mass_2.ok_or_else(|| missing("mass_2"))?;
    Ok((binary, options))
}

/// The settings of LALSimulation's `LALDict` that the models take
/// (`waveformrs::WaveformOptions`).
///
/// Every option is a keyword argument, here and in the waveform functions:
///
/// * ``precession_version``: the Euler angles of IMRPhenomXP/XPHM, LAL's
///   ``PhenomXPrecVersion``: 102 (``"NNLO"``), 223 (``"MSA"``) or 320
///   (``"SpinTaylor"``, the default here; LALSimulation 6.2.0 defaults to
///   223, see ``WaveformOptions.lalsimulation_defaults()``), or
///   ``"SpinTaylorContinuous"``, the SpinTaylor angles without their jumps
///   in the masses and spins (reported as 320).
/// * ``mode_array``: LAL's ``ModeArray``, a sequence of ``(l, m)`` pairs, or
///   ``None`` for all modes of the model.
/// * ``version_fallback``: whether IMRPhenomXP/XPHM/XPNR may change the
///   precession version automatically as LALSimulation does (default
///   ``True``); with ``False`` they raise ``PrecessionVersionUnavailableError``
///   instead.
/// * ``lambda_1``, ``lambda_2``: LAL's ``lambda1``, ``lambda2``, the
///   dimensionless tidal deformabilities of the bodies of ``mass_1`` and
///   ``mass_2`` (default 0), used by the NRTidal approximants
///   (``IMRPhenomD_NRTidalv2``, ``IMRPhenomXAS_NRTidalv2``,
///   ``IMRPhenomXAS_NRTidalv3``, ``IMRPhenomXP_NRTidalv2``,
///   ``IMRPhenomXP_NRTidalv3``); the other approximants reject non-zero
///   values.
///
/// Multibanding is not implemented: the models reproduce LAL with
/// ``PhenomXHMThresholdMband = PhenomXPHMThresholdMband = 0``.
#[pyclass(name = "WaveformOptions", module = "waveformrs", frozen, eq)]
#[derive(Debug, PartialEq)]
pub struct PyWaveformOptions {
    pub inner: WaveformOptions,
}

#[pymethods]
impl PyWaveformOptions {
    #[new]
    #[pyo3(
        signature = (**kwargs),
        text_signature = "(*, precession_version=320, mode_array=None, version_fallback=True, \
                          lambda_1=0.0, lambda_2=0.0)"
    )]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let keywords = match kwargs {
            Some(kwargs) => kwargs
                .iter()
                .map(|(key, value)| Ok((key.cast_into::<PyString>()?, value)))
                .collect::<PyResult<Vec<_>>>()?,
            None => Vec::new(),
        };
        Ok(Self {
            inner: options_with(None, keywords)?,
        })
    }

    /// The defaults of LALSimulation 6.2.0 for an empty ``LALDict``: the MSA
    /// Euler angles (223) and all modes.
    #[staticmethod]
    fn lalsimulation_defaults() -> Self {
        Self {
            inner: WaveformOptions::lalsimulation_defaults(),
        }
    }

    /// A copy with the options given as keyword arguments changed.
    #[pyo3(signature = (**kwargs), text_signature = "($self, **options)")]
    fn replace(slf: &Bound<'_, Self>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let keywords = match kwargs {
            Some(kwargs) => kwargs
                .iter()
                .map(|(key, value)| Ok((key.cast_into::<PyString>()?, value)))
                .collect::<PyResult<Vec<_>>>()?,
            None => Vec::new(),
        };
        Ok(Self {
            inner: options_with(Some(slf.as_any()), keywords)?,
        })
    }

    /// LAL's ``PhenomXPrecVersion`` of the Euler angles (102, 223 or 320).
    #[getter]
    fn precession_version(&self) -> u32 {
        self.inner.precession_version.lal_flag()
    }

    /// The mode array as a list of ``(l, m)``, or ``None`` for all modes.
    #[getter]
    fn mode_array(&self) -> Option<Vec<(u32, i32)>> {
        self.inner.mode_array.clone()
    }

    /// Whether the precession version may change automatically.
    #[getter]
    fn version_fallback(&self) -> bool {
        self.inner.version_fallback
    }

    /// The tidal deformability of the body of ``mass_1`` (LAL's ``lambda1``).
    #[getter]
    fn lambda_1(&self) -> f64 {
        self.inner.lambda_1
    }

    /// The tidal deformability of the body of ``mass_2`` (LAL's ``lambda2``).
    #[getter]
    fn lambda_2(&self) -> f64 {
        self.inner.lambda_2
    }

    /// The options as a dict of keyword arguments.
    fn as_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (key, value) in option_items(py, &self.inner)? {
            dict.set_item(key, value)?;
        }
        Ok(dict)
    }

    /// Pickling: rebuilt from its keyword arguments.
    fn __getnewargs_ex__<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<(Bound<'py, PyTuple>, Bound<'py, PyDict>)> {
        Ok((PyTuple::empty(py), self.as_dict(py)?))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let items = option_items(py, &self.inner)?
            .into_iter()
            .map(|(key, value)| Ok(format!("{key}={}", value.repr()?)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(format!("WaveformOptions({})", items.join(", ")))
    }
}
