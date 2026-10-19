//! The result classes: numpy arrays with their sampling metadata.
//!
//! Each class unpacks like a tuple of its main arrays (`hplus, hcross =
//! result`, `frequencies, modes = result`) and keeps the rest as attributes.

use num_complex::Complex64;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyIterator, PyList, PyTuple};
use waveformrs::{FrequencyMode, FrequencyModes};

use crate::errors::to_pyerr_with_prefix;

/// `(h_+, h_x)` as complex numpy arrays.
pub type Pair<'py> = (
    Bound<'py, PyArray1<Complex64>>,
    Bound<'py, PyArray1<Complex64>>,
);

/// An iterator over `items`, for `__iter__`.
fn iterate<'py>(py: Python<'py>, items: &[&Py<PyAny>]) -> PyResult<Bound<'py, PyIterator>> {
    PyTuple::new(py, items.iter().map(|item| item.bind(py)))?
        .into_any()
        .try_iter()
}

/// `epoch + j delta_t` for `j < length`.
fn sample_times(
    py: Python<'_>,
    epoch: f64,
    delta_t: f64,
    length: usize,
) -> Bound<'_, PyArray1<f64>> {
    (0..length)
        .map(|j| epoch + j as f64 * delta_t)
        .collect::<Vec<f64>>()
        .into_pyarray(py)
}

/// Time-domain polarizations ``h_+(t)``, ``h_x(t)`` sampled every
/// ``delta_t`` seconds from ``epoch`` (seconds, relative to the reference
/// time of the model; rounded to nanoseconds as LAL's ``LIGOTimeGPS``).
///
/// Unpacks as ``hplus, hcross = result``.
#[pyclass(name = "TimeSeriesPolarizations", module = "waveformrs", frozen)]
pub struct PyTimeSeriesPolarizations {
    /// ``h_+`` (float64).
    #[pyo3(get)]
    hplus: Py<PyAny>,
    /// ``h_x`` (float64).
    #[pyo3(get)]
    hcross: Py<PyAny>,
    /// Time of the first sample in seconds.
    #[pyo3(get)]
    epoch: f64,
    /// Sampling interval in seconds.
    #[pyo3(get)]
    delta_t: f64,
    length: usize,
}

impl PyTimeSeriesPolarizations {
    pub fn new(
        py: Python<'_>,
        hplus: Vec<f64>,
        hcross: Vec<f64>,
        epoch: f64,
        delta_t: f64,
    ) -> Self {
        let length = hplus.len();
        Self {
            hplus: hplus.into_pyarray(py).into_any().unbind(),
            hcross: hcross.into_pyarray(py).into_any().unbind(),
            epoch,
            delta_t,
            length,
        }
    }
}

#[pymethods]
impl PyTimeSeriesPolarizations {
    /// The sample times ``epoch + j * delta_t`` in seconds.
    #[getter]
    fn times<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        sample_times(py, self.epoch, self.delta_t, self.length)
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        iterate(py, &[&self.hplus, &self.hcross])
    }

    fn __repr__(&self) -> String {
        format!(
            "TimeSeriesPolarizations(length={}, epoch={}, delta_t={})",
            self.length, self.epoch, self.delta_t
        )
    }
}

/// Frequency-domain polarizations ``h_+(f)``, ``h_x(f)`` on the grid ``k *
/// delta_f``, ``k = 0, 1, ...``, with the ``epoch`` (seconds) of the
/// corresponding time series, as ``SimInspiralFD`` returns them.
///
/// Unpacks as ``hplus, hcross = result``.
#[pyclass(name = "FrequencySeriesPolarizations", module = "waveformrs", frozen)]
pub struct PyFrequencySeriesPolarizations {
    /// ``h_+`` (complex128).
    #[pyo3(get)]
    hplus: Py<PyAny>,
    /// ``h_x`` (complex128).
    #[pyo3(get)]
    hcross: Py<PyAny>,
    /// Epoch of the series in seconds.
    #[pyo3(get)]
    epoch: f64,
    /// Frequency spacing in Hz.
    #[pyo3(get)]
    delta_f: f64,
    length: usize,
}

impl PyFrequencySeriesPolarizations {
    pub fn new(
        py: Python<'_>,
        hplus: Vec<Complex64>,
        hcross: Vec<Complex64>,
        epoch: f64,
        delta_f: f64,
    ) -> Self {
        let length = hplus.len();
        Self {
            hplus: hplus.into_pyarray(py).into_any().unbind(),
            hcross: hcross.into_pyarray(py).into_any().unbind(),
            epoch,
            delta_f,
            length,
        }
    }
}

#[pymethods]
impl PyFrequencySeriesPolarizations {
    /// The frequencies ``k * delta_f`` in Hz.
    #[getter]
    fn frequencies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        sample_times(py, 0.0, self.delta_f, self.length)
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        iterate(py, &[&self.hplus, &self.hcross])
    }

    fn __repr__(&self) -> String {
        format!(
            "FrequencySeriesPolarizations(length={}, epoch={}, delta_f={})",
            self.length, self.epoch, self.delta_f
        )
    }
}

/// Frequency-domain modes ``h_lm(f)`` without the spherical harmonics, in
/// the layout of LAL's ``SphHarmFrequencySeries`` (``SimInspiralChooseFDModes``):
/// each mode on the two-sided frequencies ``-f_{n-1}, ..., -f_0, f_0, ...,
/// f_{n-1}`` (with ``f = 0`` once when ``f_0 = 0``).
///
/// Unpacks as ``frequencies, modes = result``; ``result[l, m]`` is one mode.
///
/// The polarizations follow from the modes with the spin-weighted spherical
/// harmonics at ``(theta, phi)`` and a rotation of the polarization basis by
/// ``polarization_angle`` (``polarizations()``): for the aligned-spin models
/// ``theta`` is the inclination and ``phi = pi / 2``; the precessing models
/// return their modes in LAL's inertial J frame.
#[pyclass(name = "FrequencyModes", module = "waveformrs", frozen)]
pub struct PyFrequencyModes {
    /// The two-sided frequencies in Hz (ascending).
    #[pyo3(get)]
    frequencies: Py<PyAny>,
    /// ``{(l, m): h_lm}`` (complex128 arrays on ``frequencies``).
    #[pyo3(get)]
    modes: Py<PyDict>,
    /// Polar angle of the line of sight in the frame of the modes.
    #[pyo3(get)]
    theta: f64,
    /// Azimuth of the line of sight in the frame of the modes.
    #[pyo3(get)]
    phi: f64,
    /// Rotation of the polarization basis (zero for the aligned-spin models).
    #[pyo3(get)]
    polarization_angle: f64,
    /// The number of non-negative frequencies.
    n_positive: usize,
}

impl PyFrequencyModes {
    pub fn new(py: Python<'_>, modes: FrequencyModes) -> PyResult<Self> {
        let n = modes.frequencies.len();
        // a zero frequency appears once, as in LAL's two-sided series
        let skip = usize::from(n > 0 && modes.frequencies[0] == 0.0);
        let mut frequencies = Vec::with_capacity(2 * n - skip);
        frequencies.extend(modes.frequencies[skip..].iter().rev().map(|f| -f));
        frequencies.extend_from_slice(&modes.frequencies);
        let dict = PyDict::new(py);
        for mode in modes.modes {
            let mut h = Vec::with_capacity(2 * n - skip);
            h.extend(mode.negative[skip..].iter().rev());
            h.extend_from_slice(&mode.positive);
            dict.set_item((mode.ell, mode.m), h.into_pyarray(py))?;
        }
        Ok(Self {
            frequencies: frequencies.into_pyarray(py).into_any().unbind(),
            modes: dict.unbind(),
            theta: modes.theta,
            phi: modes.phi,
            polarization_angle: modes.polarization_angle,
            n_positive: n,
        })
    }
}

#[pymethods]
impl PyFrequencyModes {
    /// The polarizations ``(h_+, h_x)`` at the non-negative frequencies
    /// (``frequencies[-n:]``) from the modes in ``modes`` (which may be
    /// edited to leave modes out).
    fn polarizations<'py>(&self, py: Python<'py>) -> PyResult<Pair<'py>> {
        let frequencies = self
            .frequencies
            .bind(py)
            .extract::<PyReadonlyArray1<'py, f64>>()?;
        let frequencies = frequencies.as_slice()?;
        let n2 = frequencies.len();
        let n = self.n_positive;
        let mut modes = Vec::new();
        for (key, value) in self.modes.bind(py).iter() {
            let (ell, m): (u32, i32) = key.extract()?;
            let h = value
                .extract::<PyReadonlyArray1<'py, Complex64>>()
                .map_err(|_| {
                    PyValueError::new_err(format!(
                        "mode ({ell}, {m}) must be a complex128 array on the frequencies"
                    ))
                })?;
            let h = h.as_slice()?;
            if h.len() != n2 {
                return Err(PyValueError::new_err(format!(
                    "mode ({ell}, {m}) has {} samples, not {n2}",
                    h.len()
                )));
            }
            modes.push(FrequencyMode {
                ell,
                m,
                positive: h[n2 - n..].to_vec(),
                // h(-f_k) = h[n - 1 - k], also with a single f = 0 sample
                negative: h[..n].iter().rev().copied().collect(),
            });
        }
        let modes = FrequencyModes {
            frequencies: frequencies[n2 - n..].to_vec(),
            theta: self.theta,
            phi: self.phi,
            polarization_angle: self.polarization_angle,
            modes,
        };
        let (hplus, hcross) = py.detach(|| modes.polarizations());
        Ok((hplus.into_pyarray(py), hcross.into_pyarray(py)))
    }

    fn __getitem__<'py>(&self, py: Python<'py>, key: (i64, i64)) -> PyResult<Bound<'py, PyAny>> {
        self.modes
            .bind(py)
            .get_item(key)?
            .ok_or_else(|| PyKeyError::new_err(key))
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let modes = self.modes.clone_ref(py).into_any();
        iterate(py, &[&self.frequencies, &modes])
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        let keys: Vec<String> = self
            .modes
            .bind(py)
            .keys()
            .iter()
            .map(|key| key.to_string())
            .collect();
        format!(
            "FrequencyModes(n_frequencies={}, modes=[{}])",
            self.n_positive,
            keys.join(", ")
        )
    }
}

/// Time-domain modes ``h_lm(t)`` (complex, without the spherical
/// harmonics) sampled every ``delta_t`` seconds from ``epoch``, as
/// ``SimInspiralChooseTDModes``.
///
/// Unpacks as ``times, modes = result``; ``result[l, m]`` is one mode.
#[pyclass(name = "TimeSeriesModes", module = "waveformrs", frozen)]
pub struct PyTimeSeriesModes {
    /// ``{(l, m): h_lm}`` (complex128 arrays).
    #[pyo3(get)]
    modes: Py<PyDict>,
    /// Time of the first sample in seconds.
    #[pyo3(get)]
    epoch: f64,
    /// Sampling interval in seconds.
    #[pyo3(get)]
    delta_t: f64,
    length: usize,
}

impl PyTimeSeriesModes {
    pub fn new(py: Python<'_>, modes: waveformrs::model::TimeSeriesModes) -> PyResult<Self> {
        let length = modes.modes.first().map_or(0, |mode| mode.data.len());
        let dict = PyDict::new(py);
        for mode in modes.modes {
            dict.set_item((mode.ell, mode.m), mode.data.into_pyarray(py))?;
        }
        Ok(Self {
            modes: dict.unbind(),
            epoch: modes.epoch,
            delta_t: modes.delta_t,
            length,
        })
    }
}

#[pymethods]
impl PyTimeSeriesModes {
    /// The sample times ``epoch + j * delta_t`` in seconds.
    #[getter]
    fn times<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        sample_times(py, self.epoch, self.delta_t, self.length)
    }

    fn __getitem__<'py>(&self, py: Python<'py>, key: (i64, i64)) -> PyResult<Bound<'py, PyAny>> {
        self.modes
            .bind(py)
            .get_item(key)?
            .ok_or_else(|| PyKeyError::new_err(key))
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let times = self.times(py).into_any().unbind();
        let modes = self.modes.clone_ref(py).into_any();
        iterate(py, &[&times, &modes])
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "TimeSeriesModes(length={}, epoch={}, delta_t={}, n_modes={})",
            self.length,
            self.epoch,
            self.delta_t,
            self.modes.bind(py).len()
        )
    }
}

/// The success mask and the exceptions of a batch (raising the first
/// exception if `strict`).
pub fn batch_errors<'py>(
    py: Python<'py>,
    results: &[waveformrs::Result<()>],
    strict: bool,
) -> PyResult<(Bound<'py, PyArray1<bool>>, Bound<'py, PyDict>)> {
    if strict {
        if let Some((i, Err(err))) = results.iter().enumerate().find(|(_, r)| r.is_err()) {
            return Err(to_pyerr_with_prefix(err, &format!("binary {i}: ")));
        }
    }
    let errors = PyDict::new(py);
    for (i, result) in results.iter().enumerate() {
        if let Err(err) = result {
            errors.set_item(i, to_pyerr_with_prefix(err, "").into_value(py))?;
        }
    }
    let ok: Vec<bool> = results.iter().map(Result::is_ok).collect();
    Ok((ok.into_pyarray(py), errors))
}

/// Frequency-domain polarizations of a batch of binaries: row ``i`` of
/// ``hplus`` and ``hcross`` (``(n_binaries, n_frequencies)`` complex128
/// arrays) is binary ``i``.
///
/// ``ok[i]`` is ``False`` for the binaries that failed (their rows are
/// zero) and ``errors[i]`` the exception they raised.
///
/// Unpacks as ``hplus, hcross = result``.
#[pyclass(name = "BatchPolarizations", module = "waveformrs", frozen)]
pub struct PyBatchPolarizations {
    /// ``h_+``, one row per binary.
    #[pyo3(get)]
    hplus: Py<PyAny>,
    /// ``h_x``, one row per binary.
    #[pyo3(get)]
    hcross: Py<PyAny>,
    /// Whether each binary succeeded.
    #[pyo3(get)]
    ok: Py<PyAny>,
    /// ``{index: exception}`` for the binaries that failed.
    #[pyo3(get)]
    errors: Py<PyDict>,
}

impl PyBatchPolarizations {
    pub fn new(
        hplus: Bound<'_, PyArray2<Complex64>>,
        hcross: Bound<'_, PyArray2<Complex64>>,
        ok: Bound<'_, PyArray1<bool>>,
        errors: Bound<'_, PyDict>,
    ) -> Self {
        Self {
            hplus: hplus.into_any().unbind(),
            hcross: hcross.into_any().unbind(),
            ok: ok.into_any().unbind(),
            errors: errors.unbind(),
        }
    }
}

#[pymethods]
impl PyBatchPolarizations {
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        iterate(py, &[&self.hplus, &self.hcross])
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let shape: Vec<usize> = self.hplus.bind(py).getattr("shape")?.extract()?;
        Ok(format!(
            "BatchPolarizations(n_binaries={}, n_frequencies={}, failed={})",
            shape[0],
            shape[1],
            self.errors.bind(py).len()
        ))
    }
}

/// Time-domain polarizations of a batch of binaries: ``waveforms[i]`` is the
/// ``TimeSeriesPolarizations`` of binary ``i``, or ``None`` if it failed
/// (``ok[i]`` is then ``False`` and ``errors[i]`` the exception).
///
/// Iterates over ``waveforms``.
#[pyclass(name = "BatchTimeSeries", module = "waveformrs", frozen)]
pub struct PyBatchTimeSeries {
    /// One ``TimeSeriesPolarizations`` (or ``None``) per binary.
    #[pyo3(get)]
    waveforms: Py<PyList>,
    /// Whether each binary succeeded.
    #[pyo3(get)]
    ok: Py<PyAny>,
    /// ``{index: exception}`` for the binaries that failed.
    #[pyo3(get)]
    errors: Py<PyDict>,
}

impl PyBatchTimeSeries {
    pub fn new(
        waveforms: Bound<'_, PyList>,
        ok: Bound<'_, PyArray1<bool>>,
        errors: Bound<'_, PyDict>,
    ) -> Self {
        Self {
            waveforms: waveforms.unbind(),
            ok: ok.into_any().unbind(),
            errors: errors.unbind(),
        }
    }
}

#[pymethods]
impl PyBatchTimeSeries {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.waveforms.bind(py).len()
    }

    fn __getitem__<'py>(
        &self,
        py: Python<'py>,
        index: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.waveforms.bind(py).as_any().get_item(index)
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        self.waveforms.bind(py).as_any().try_iter()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "BatchTimeSeries(n_binaries={}, failed={})",
            self.waveforms.bind(py).len(),
            self.errors.bind(py).len()
        )
    }
}
