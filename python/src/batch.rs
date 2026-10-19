//! Many binaries at once with `waveformrs::batch::Batch`.
//!
//! The binaries come as a table with one row per binary and the columns of
//! [`COLUMNS`] (the package's `batch_polarizations` builds it from a dict
//! of arrays or a structured array). The frequency-domain results are
//! written by the batch threads straight into the rows of numpy arrays,
//! with the GIL released.

use num_complex::Complex64;
use numpy::{PyArray2, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use waveformrs::batch::Batch;
use waveformrs::nrsur::TimeSeriesPolarizations;
use waveformrs::{choose_fd_model, Binary, Error};

use crate::arrays::f64_array2;
use crate::errors::to_pyerr;
use crate::grid_from;
use crate::models::{guard, td_model, Guarded, Sequence};
use crate::params::{parse_approximant, parse_options};
use crate::results::{
    batch_errors, PyBatchPolarizations, PyBatchTimeSeries, PyTimeSeriesPolarizations,
};

/// The columns of a table of binaries.
const COLUMNS: [&str; 12] = [
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
];

/// The binaries of a table with the columns of [`COLUMNS`].
fn binaries(table: &PyReadonlyArray2<'_, f64>) -> PyResult<Vec<Binary>> {
    let shape = table.shape();
    if shape[1] != COLUMNS.len() {
        return Err(PyValueError::new_err(format!(
            "parameters must have {} columns ({}), not {}",
            COLUMNS.len(),
            COLUMNS.join(", "),
            shape[1]
        )));
    }
    Ok(table
        .as_slice()?
        .chunks_exact(COLUMNS.len())
        .map(|row| Binary {
            mass_1: row[0],
            mass_2: row[1],
            spin_1: [row[2], row[3], row[4]],
            spin_2: [row[5], row[6], row[7]],
            luminosity_distance: row[8],
            inclination: row[9],
            phi_ref: row[10],
            f_ref: row[11],
        })
        .collect())
}

/// The frequency-domain polarizations of many binaries on one grid.
///
/// `parameters` is a `(n_binaries, 12)` table (see the package's
/// `batch_polarizations`); the result has one row per binary.
#[pyfunction]
#[pyo3(signature = (approximant, grid, parameters, *, threads = 1, strict = false, exact = false, **kwargs))]
#[allow(clippy::too_many_arguments)]
pub fn batch_polarizations<'py>(
    py: Python<'py>,
    approximant: &str,
    grid: &Bound<'py, PyAny>,
    parameters: &Bound<'py, PyAny>,
    threads: usize,
    strict: bool,
    exact: bool,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyBatchPolarizations> {
    let approximant = parse_approximant(approximant)?;
    if !approximant.is_frequency_domain() {
        return Err(to_pyerr(Error::Unsupported {
            approximant,
            what:
                "frequency-domain waveforms (it is a time-domain model: use batch_td_polarizations)",
        }));
    }
    let options = parse_options(kwargs)?;
    let table = f64_array2(parameters, "parameters")?;
    let binaries = binaries(&table)?;
    let grid = grid_from(grid)?;
    let grid = &grid.get().inner;
    let hplus = PyArray2::<Complex64>::zeros(py, [binaries.len(), grid.len()], false);
    let hcross = PyArray2::<Complex64>::zeros(py, [binaries.len(), grid.len()], false);
    let results = {
        let mut hp = hplus.readwrite();
        let mut hc = hcross.readwrite();
        let (hp, hc) = (hp.as_slice_mut()?, hc.as_slice_mut()?);
        let batch = Batch::new(threads);
        py.detach(|| {
            if exact {
                batch.polarizations_on_grid(
                    &binaries,
                    |binary| {
                        guard(|| choose_fd_model(approximant, binary, &options))
                            .map(|model| Guarded(Sequence(model)))
                    },
                    grid,
                    hp,
                    hc,
                )
            } else {
                batch.polarizations_on_grid(
                    &binaries,
                    |binary| guard(|| choose_fd_model(approximant, binary, &options)).map(Guarded),
                    grid,
                    hp,
                    hc,
                )
            }
        })
    };
    let (ok, errors) = batch_errors(py, &results, strict)?;
    Ok(PyBatchPolarizations::new(hplus, hcross, ok, errors))
}

/// The time-domain polarizations (`SimInspiralChooseTDWaveform`) of many
/// binaries, with one NRSur7dq4 workspace per thread.
#[pyfunction]
#[pyo3(signature = (approximant, delta_t, f_min, parameters, *, threads = 1, strict = false, **kwargs))]
#[allow(clippy::too_many_arguments)]
pub fn batch_td_polarizations<'py>(
    py: Python<'py>,
    approximant: &str,
    delta_t: f64,
    f_min: f64,
    parameters: &Bound<'py, PyAny>,
    threads: usize,
    strict: bool,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyBatchTimeSeries> {
    let approximant = parse_approximant(approximant)?;
    let options = parse_options(kwargs)?;
    let table = f64_array2(parameters, "parameters")?;
    let binaries = binaries(&table)?;
    let mut out: Vec<TimeSeriesPolarizations> =
        (0..binaries.len()).map(|_| Default::default()).collect();
    let results = py.detach(|| {
        Batch::new(threads).time_domain_polarizations(
            &binaries,
            |binary| guard(|| td_model(approximant, binary, &options)).map(Guarded),
            delta_t,
            f_min,
            &mut out,
        )
    });
    let (ok, errors) = batch_errors(py, &results, strict)?;
    let waveforms = PyList::empty(py);
    for (h, result) in out.into_iter().zip(&results) {
        if result.is_ok() {
            let h = PyTimeSeriesPolarizations::new(py, h.hplus, h.hcross, h.epoch, h.delta_t);
            waveforms.append(Bound::new(py, h)?)?;
        } else {
            waveforms.append(py.None())?;
        }
    }
    Ok(PyBatchTimeSeries::new(waveforms, ok, errors))
}
