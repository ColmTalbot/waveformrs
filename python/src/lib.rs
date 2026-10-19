//! Python bindings of `waveformrs` (the extension module
//! `waveformrs._waveformrs`; the package `waveformrs` re-exports it).
//!
//! Every function releases the GIL while the waveforms are computed and
//! writes the results straight into numpy arrays where their size is known
//! in advance.

mod arrays;
mod batch;
mod errors;
mod models;
mod params;
mod results;

use num_complex::Complex64;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use waveformrs::grid::FrequencyGrid;
use waveformrs::nrsur::TimeSeriesPolarizations;
use waveformrs::{choose_fd_model, conversions, Approximant, TimeDomainModel};

use crate::arrays::{check_frequencies, f64_array1};
use crate::errors::IntoPyResult;
use crate::models::{td_model, with_workspace};
use crate::params::{parse_approximant, parse_request, PyWaveformOptions};
use crate::results::{
    Pair, PyBatchPolarizations, PyBatchTimeSeries, PyFrequencyModes,
    PyFrequencySeriesPolarizations, PyTimeSeriesModes, PyTimeSeriesPolarizations,
};

/// The names of the approximants (LAL's names), all of them or those of one
/// domain (``"frequency"`` or ``"time"``).
#[pyfunction]
#[pyo3(signature = (domain = None))]
fn approximants(domain: Option<&str>) -> PyResult<Vec<&'static str>> {
    let keep = |approximant: &Approximant| match domain {
        None => Ok(true),
        Some("frequency") => Ok(approximant.is_frequency_domain()),
        Some("time") => Ok(!approximant.is_frequency_domain()),
        Some(other) => Err(PyValueError::new_err(format!(
            "domain must be \"frequency\", \"time\" or None, not {other:?}"
        ))),
    };
    let mut names = Vec::new();
    for approximant in Approximant::ALL {
        if keep(&approximant)? {
            names.push(approximant.name());
        }
    }
    Ok(names)
}

/// Whether the approximant accepts in-plane spins.
#[pyfunction]
fn is_precessing(approximant: &str) -> PyResult<bool> {
    Ok(parse_approximant(approximant)?.is_precessing())
}

/// Whether the approximant is a frequency-domain model (otherwise it is a
/// time-domain model).
#[pyfunction]
fn is_frequency_domain(approximant: &str) -> PyResult<bool> {
    Ok(parse_approximant(approximant)?.is_frequency_domain())
}

/// The NRSur7dq4 data file found through the environment variables
/// ``NRSUR7DQ4_DATA`` (the file) or ``LAL_DATA_PATH`` (its directory), or
/// ``None``.
#[pyfunction]
fn nrsur7dq4_data_file() -> Option<String> {
    waveformrs::nrsur::find_data_file().map(|path| path.to_string_lossy().into_owned())
}

/// Frequencies (Hz, positive and increasing) prepared for evaluating many
/// waveforms: the powers of the frequency that the models need are computed
/// once. Results on a grid agree with ``polarizations`` to rounding (about
/// 1e-13 relative) rather than bit for bit.
#[pyclass(name = "FrequencyGrid", module = "waveformrs", frozen)]
pub struct PyFrequencyGrid {
    pub inner: FrequencyGrid,
}

#[pymethods]
impl PyFrequencyGrid {
    #[new]
    fn new(frequencies: &Bound<'_, PyAny>) -> PyResult<Self> {
        let frequencies = f64_array1(frequencies, "frequencies")?;
        let frequencies = frequencies.as_slice()?;
        check_frequencies(frequencies)?;
        Ok(Self {
            inner: FrequencyGrid::new(frequencies),
        })
    }

    /// The uniform grid ``f_min, f_min + delta_f, ...`` up to and including
    /// ``f_max`` (to within half a step).
    #[staticmethod]
    fn uniform(f_min: f64, f_max: f64, delta_f: f64) -> PyResult<Self> {
        if !(delta_f > 0.0 && f_min > 0.0 && f_max >= f_min && f_max.is_finite()) {
            return Err(errors::invalid_parameter(
                "the uniform grid needs delta_f > 0 and finite 0 < f_min <= f_max",
            ));
        }
        Ok(Self {
            inner: FrequencyGrid::uniform(f_min, f_max, delta_f),
        })
    }

    /// The frequencies in Hz (a copy).
    #[getter]
    fn frequencies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, self.inner.frequencies())
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Pickling: rebuilt from its frequencies.
    fn __getnewargs__<'py>(&self, py: Python<'py>) -> (Bound<'py, PyArray1<f64>>,) {
        (self.frequencies(py),)
    }

    fn __repr__(&self) -> String {
        let f = self.inner.frequencies();
        match (f.first(), f.last()) {
            (Some(first), Some(last)) => {
                format!("FrequencyGrid(n={}, from {first} Hz to {last} Hz)", f.len())
            }
            _ => "FrequencyGrid(n=0)".to_string(),
        }
    }
}

/// Two zeroed complex arrays of length `n`.
fn zeros_pair(py: Python<'_>, n: usize) -> Pair<'_> {
    (
        PyArray1::<Complex64>::zeros(py, n, false),
        PyArray1::<Complex64>::zeros(py, n, false),
    )
}

/// ``(h_+, h_x)`` of ``approximant`` at ``frequencies`` (Hz, positive and
/// increasing), as ``SimInspiralChooseFDWaveformSequence``.
///
/// The binary is given by keyword arguments (``mass_1``, ``mass_2`` in
/// solar masses; ``spin_1``, ``spin_2`` as ``(x, y, z)``; ``luminosity_distance``
/// in Mpc; ``inclination``, ``phi_ref`` in rad; ``f_ref`` in Hz), followed by
/// the options (``options=WaveformOptions(...)`` and/or its keywords).
#[pyfunction]
#[pyo3(
    signature = (approximant, frequencies, **kwargs),
    text_signature = "(approximant, frequencies, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn polarizations<'py>(
    py: Python<'py>,
    approximant: &str,
    frequencies: &Bound<'py, PyAny>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Pair<'py>> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let frequencies = f64_array1(frequencies, "frequencies")?;
    let f = frequencies.as_slice()?;
    check_frequencies(f)?;
    let (hplus, hcross) = zeros_pair(py, f.len());
    {
        let mut hp = hplus.readwrite();
        let mut hc = hcross.readwrite();
        let (hp, hc) = (hp.as_slice_mut()?, hc.as_slice_mut()?);
        py.detach(|| {
            choose_fd_model(approximant, &binary, &options)?.polarizations_into(f, hp, hc)
        })
        .py_result()?;
    }
    Ok((hplus, hcross))
}

/// A `FrequencyGrid` from a `FrequencyGrid` or an array of frequencies.
pub(crate) fn grid_from<'py>(grid: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyFrequencyGrid>> {
    match grid.cast::<PyFrequencyGrid>() {
        Ok(grid) => Ok(grid.clone()),
        Err(_) => Bound::new(grid.py(), PyFrequencyGrid::new(grid)?),
    }
}

/// ``(h_+, h_x)`` of ``approximant`` on a ``FrequencyGrid`` (or an array of
/// frequencies), with the fast paths for repeated evaluation on the same
/// frequencies. Keyword arguments as for ``polarizations``.
#[pyfunction]
#[pyo3(
    signature = (approximant, grid, **kwargs),
    text_signature = "(approximant, grid, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn polarizations_on_grid<'py>(
    py: Python<'py>,
    approximant: &str,
    grid: &Bound<'py, PyAny>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Pair<'py>> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let grid = grid_from(grid)?;
    let grid = &grid.get().inner;
    let (hplus, hcross) = zeros_pair(py, grid.len());
    {
        let mut hp = hplus.readwrite();
        let mut hc = hcross.readwrite();
        let (hp, hc) = (hp.as_slice_mut()?, hc.as_slice_mut()?);
        py.detach(|| {
            choose_fd_model(approximant, &binary, &options)?.polarizations_on_grid(grid, hp, hc)
        })
        .py_result()?;
    }
    Ok((hplus, hcross))
}

/// The modes ``h_lm`` (without the spherical harmonics) of ``approximant``
/// at ``frequencies`` (Hz, positive and increasing) and at ``-frequencies``,
/// in the conventions of ``SimInspiralChooseFDModes``: IMRPhenomXHM its
/// modes, IMRPhenomXPHM/XPNR the J-frame modes, TaylorF2, IMRPhenomD and
/// IMRPhenomXAS their (2, +-2) modes and IMRPhenomXP its J-frame (2, m)
/// modes. Keyword arguments as for ``polarizations``.
#[pyfunction]
#[pyo3(
    signature = (approximant, frequencies, **kwargs),
    text_signature = "(approximant, frequencies, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn modes(
    py: Python<'_>,
    approximant: &str,
    frequencies: &Bound<'_, PyAny>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyFrequencyModes> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let frequencies = f64_array1(frequencies, "frequencies")?;
    let f = frequencies.as_slice()?;
    check_frequencies(f)?;
    let modes = py
        .detach(|| choose_fd_model(approximant, &binary, &options)?.modes(f))
        .py_result()?;
    PyFrequencyModes::new(py, modes)
}

/// ``SimInspiralChooseFDModes`` (IMRPhenomXHM, IMRPhenomXPHM and
/// IMRPhenomXPNR): the modes on the two-sided uniform grid ``k * delta_f``,
/// ``|k| < n``, zero below ``f_min`` and above ``min(f_max, f_cut)``
/// (``f_max = 0`` selects the model's cutoff). Keyword arguments as for
/// ``polarizations``; ``f_ref = 0`` selects ``f_min``.
#[pyfunction]
#[pyo3(
    signature = (approximant, delta_f, f_min, f_max = 0.0, **kwargs),
    text_signature = "(approximant, delta_f, f_min, f_max=0.0, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn choose_fd_modes(
    py: Python<'_>,
    approximant: &str,
    delta_f: f64,
    f_min: f64,
    f_max: f64,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyFrequencyModes> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let modes = py
        .detach(|| {
            waveformrs::choose_fd_modes(approximant, &binary, delta_f, f_min, f_max, &options)
        })
        .py_result()?;
    PyFrequencyModes::new(py, modes)
}

/// ``SimInspiralChooseTDWaveform``: the time-domain polarizations sampled
/// every ``delta_t`` seconds from the (2, 2) frequency ``f_min`` (Hz; zero
/// for the start of the model for NRSur7dq4). NRSur7dq4 is evaluated with a
/// workspace kept per thread (no allocations and no new surrogate
/// evaluation when only the extrinsic parameters change); the
/// frequency-domain approximants (except TaylorF2) are transformed as LAL
/// does (``XLALSimInspiralTDFromFD``). Keyword arguments as for
/// ``polarizations``.
#[pyfunction]
#[pyo3(
    signature = (approximant, delta_t, f_min, **kwargs),
    text_signature = "(approximant, delta_t, f_min, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn td_polarizations(
    py: Python<'_>,
    approximant: &str,
    delta_t: f64,
    f_min: f64,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyTimeSeriesPolarizations> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let h = py
        .detach(|| {
            let model = td_model(approximant, &binary, &options)?;
            let mut out = TimeSeriesPolarizations::default();
            with_workspace(|workspace| {
                model.polarizations_into(workspace, delta_t, f_min, &mut out)
            })?;
            Ok(out)
        })
        .py_result()?;
    Ok(PyTimeSeriesPolarizations::new(
        py, h.hplus, h.hcross, h.epoch, h.delta_t,
    ))
}

/// ``SimInspiralChooseTDModes`` of NRSur7dq4: the inertial-frame modes
/// sampled every ``delta_t`` seconds from ``f_min``. Keyword arguments as
/// for ``polarizations`` (``mode_array`` selects the coorbital modes).
#[pyfunction]
#[pyo3(
    signature = (approximant, delta_t, f_min, **kwargs),
    text_signature = "(approximant, delta_t, f_min, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn td_modes(
    py: Python<'_>,
    approximant: &str,
    delta_t: f64,
    f_min: f64,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyTimeSeriesModes> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let modes = py
        .detach(|| {
            let model = td_model(approximant, &binary, &options)?;
            let mut out = waveformrs::model::TimeSeriesModes::default();
            with_workspace(|workspace| model.modes_into(workspace, delta_t, f_min, &mut out))?;
            Ok(out)
        })
        .py_result()?;
    PyTimeSeriesModes::new(py, modes)
}

/// ``SimInspiralTD``: the conditioned time-domain waveform of any
/// approximant (tapered and high-pass filtered as LAL does; transformed
/// from the frequency domain for the frequency-domain approximants).
/// Keyword arguments as for ``polarizations``.
#[pyfunction]
#[pyo3(
    signature = (approximant, delta_t, f_min, **kwargs),
    text_signature = "(approximant, delta_t, f_min, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn simulate_td(
    py: Python<'_>,
    approximant: &str,
    delta_t: f64,
    f_min: f64,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyTimeSeriesPolarizations> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let h = py
        .detach(|| conversions::simulate_td(approximant, &binary, &options, delta_t, f_min))
        .py_result()?;
    let epoch = h.hplus.epoch.seconds();
    let delta_t = h.hplus.delta_t;
    Ok(PyTimeSeriesPolarizations::new(
        py,
        h.hplus.data,
        h.hcross.data,
        epoch,
        delta_t,
    ))
}

/// ``SimInspiralFD``: the conditioned frequency-domain waveform of any
/// approximant on the grid ``k * delta_f`` up to ``f_max`` (``delta_f = 0``
/// chooses the spacing; transformed from the time domain for NRSur7dq4).
/// Keyword arguments as for ``polarizations``.
#[pyfunction]
#[pyo3(
    signature = (approximant, delta_f, f_min, f_max, **kwargs),
    text_signature = "(approximant, delta_f, f_min, f_max, *, mass_1, mass_2, spin_1=(0.0, 0.0, 0.0), spin_2=(0.0, 0.0, 0.0), luminosity_distance=1.0, inclination=0.0, phi_ref=0.0, f_ref=0.0, options=None, **option_kwargs)"
)]
fn simulate_fd(
    py: Python<'_>,
    approximant: &str,
    delta_f: f64,
    f_min: f64,
    f_max: f64,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyFrequencySeriesPolarizations> {
    let approximant = parse_approximant(approximant)?;
    let (binary, options) = parse_request(kwargs)?;
    let h = py
        .detach(|| conversions::simulate_fd(approximant, &binary, &options, delta_f, f_min, f_max))
        .py_result()?;
    let epoch = h.hplus.epoch.seconds();
    let delta_f = h.hplus.delta_f;
    Ok(PyFrequencySeriesPolarizations::new(
        py,
        h.hplus.data,
        h.hcross.data,
        epoch,
        delta_f,
    ))
}

#[pymodule]
fn _waveformrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyWaveformOptions>()?;
    m.add_class::<PyFrequencyGrid>()?;
    m.add_class::<PyTimeSeriesPolarizations>()?;
    m.add_class::<PyFrequencySeriesPolarizations>()?;
    m.add_class::<PyFrequencyModes>()?;
    m.add_class::<PyTimeSeriesModes>()?;
    m.add_class::<PyBatchPolarizations>()?;
    m.add_class::<PyBatchTimeSeries>()?;
    m.add_function(wrap_pyfunction!(approximants, m)?)?;
    m.add_function(wrap_pyfunction!(is_precessing, m)?)?;
    m.add_function(wrap_pyfunction!(is_frequency_domain, m)?)?;
    m.add_function(wrap_pyfunction!(nrsur7dq4_data_file, m)?)?;
    m.add_function(wrap_pyfunction!(polarizations, m)?)?;
    m.add_function(wrap_pyfunction!(polarizations_on_grid, m)?)?;
    m.add_function(wrap_pyfunction!(modes, m)?)?;
    m.add_function(wrap_pyfunction!(choose_fd_modes, m)?)?;
    m.add_function(wrap_pyfunction!(td_polarizations, m)?)?;
    m.add_function(wrap_pyfunction!(td_modes, m)?)?;
    m.add_function(wrap_pyfunction!(simulate_td, m)?)?;
    m.add_function(wrap_pyfunction!(simulate_fd, m)?)?;
    m.add_function(wrap_pyfunction!(batch::batch_polarizations, m)?)?;
    m.add_function(wrap_pyfunction!(batch::batch_td_polarizations, m)?)?;
    Ok(())
}
