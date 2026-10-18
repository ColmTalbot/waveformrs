//! The precessing numerical-relativity surrogate NRSur7dq4 (Varma et al.
//! 2019, arXiv:1905.09300), a port of `LALSimIMRPrecessingNRSur.c` from
//! LALSimulation 6.2.0.
//!
//! | module      | LALSimulation source                                      |
//! |-------------|-----------------------------------------------------------|
//! | `data`      | `PrecessingNRSur_Init` and the loaders (HDF5 via [`crate::hdf5`]) |
//! | `fit`       | `NRSur7dq4_eval_fit`, `NRSur7dq4_effective_spins`, `ipow` |
//! | `dynamics`  | `PrecessingNRSur_IntegrateDynamics` (RK4 start, AB4 steps, `get_t_ref`) |
//! | `frames`    | `TransformModesCoorbitalToInertial`, `WignerDMatrices_Compute` (`LALSimNRSurrogateUtilities.c`) |
//! | `spline`    | the GSL natural cubic spline and 4-point polynomial interpolation |
//! | `model`     | `PrecessingNRSur_core`, `XLALSimInspiralPrecessingNRSurPolarizations`, `XLALSimInspiralPrecessingNRSurModes`, `XLALPrecessingNRSurDynamics` |
//!
//! # Data
//!
//! The model needs the data file `NRSur7dq4_v1.0.h5` (zenodo.org/records/
//! 14999310, also in `lalsuite-waveform-data`). It is read with the
//! pure-Rust reader in [`crate::hdf5`] into an [`NRSur7dq4Data`] that is
//! loaded once (a few milliseconds) and shared by any number of waveform
//! evaluations: pass `&data` to [`NRSur7dq4::new`], wrap it in an `Arc`, or
//! use the process-wide copy of [`NRSur7dq4Data::global`], which finds the
//! file with [`find_data_file`] (`NRSUR7DQ4_DATA`, then `LAL_DATA_PATH`).
//!
//! # Accuracy
//!
//! Every operation of the C code is reproduced in the same order, with the
//! fused multiply-adds of clang's arm64 build (see [`crate::spintaylor::fp`]),
//! the GSL splines and polynomial interpolation operation by operation, the
//! fused column-by-column accumulation of the OpenBLAS `dgemv` kernel for
//! the empirical interpolants, Darwin's `__sincos_stret` where clang combines
//! `sin` and `cos` (`trig`) and the C library's `cexp` in the spherical
//! harmonics. On Apple silicon the polarizations, modes and dynamics agree
//! bit for bit with the conda-forge build of lalsimulation 6.2.0 (all
//! reference cases and 200 random draws); elsewhere expect agreement to
//! rounding (`~1e-15`).
//!
//! # Repeated evaluation
//!
//! For many waveforms in sequence (e.g. in a sampler) use
//! [`NRSur7dq4::polarizations_with`] / [`NRSur7dq4::modes_with`] with a
//! [`NRSurWorkspace`] and an output struct that are reused across calls:
//! nothing is then allocated per waveform, and everything that does not
//! depend on the binary (the fit index tables, the factorised spline
//! systems of the dynamics and coorbital time grids, the interpolation
//! weights between them and the Wigner-D coefficients) is computed once with
//! the data. The workspace also keeps the inertial-frame modes of the last
//! binary, so that changing only extrinsic parameters (inclination, phase,
//! distance, sampling interval) skips the surrogate evaluation.
//!
//! ```no_run
//! use waveformrs::nrsur::{
//!     NRSur7dq4, NRSur7dq4Data, NRSurWorkspace, PrecessingBinary, TimeSeriesPolarizations,
//! };
//!
//! // found through NRSUR7DQ4_DATA or LAL_DATA_PATH; or NRSur7dq4Data::load(path)
//! let data = NRSur7dq4Data::global()?;
//! let binary = PrecessingBinary {
//!     mass_1: 40.0,
//!     mass_2: 20.0,
//!     spin_1: [0.4, -0.3, 0.2],
//!     spin_2: [0.1, 0.5, -0.3],
//!     luminosity_distance: 100.0,
//!     inclination: 0.7,
//!     phi_ref: 0.1,
//!     f_ref: 20.0,
//! };
//! let h = NRSur7dq4::new(data, &binary)?.polarizations(1.0 / 4096.0, 20.0)?;
//! println!("{} samples from t = {} s", h.hplus.len(), h.epoch);
//!
//! // many waveforms without allocations
//! let mut workspace = NRSurWorkspace::new();
//! let mut out = TimeSeriesPolarizations::default();
//! for inclination in [0.1, 0.5, 1.0] {
//!     let model = NRSur7dq4::new(data, &PrecessingBinary { inclination, ..binary })?;
//!     model.polarizations_with(&mut workspace, 1.0 / 4096.0, 20.0, &mut out)?;
//! }
//! # Ok::<(), waveformrs::nrsur::NRSurError>(())
//! ```
//!
//! NRSur7dq2 has the same structure (with different fits, vector fits and
//! limits, see [`SurrogateVersion`]) but is not implemented.

mod data;
mod dynamics;
mod fit;
mod frames;
mod model;
mod simd;
mod spline;
mod trig;

use std::fmt;
use std::path::PathBuf;

use once_cell::sync::OnceCell;

pub use crate::imrphenomx::PrecessingBinary;
pub use data::{NRSur7dq4Data, SurrogateVersion};
pub use model::{
    Mode, NRSur7dq4, NRSurDynamics, NRSurWorkspace, TimeSeriesModes, TimeSeriesPolarizations,
};

/// File name of the NRSur7dq4 data (`NRSUR7DQ4_DATAFILE`).
pub const NRSUR7DQ4_DATAFILE: &str = "NRSur7dq4_v1.0.h5";

/// Locate `NRSur7dq4_v1.0.h5`: the file named by the environment variable
/// `NRSUR7DQ4_DATA`, else the first directory of `LAL_DATA_PATH` (a
/// colon-separated list, as used by LAL) that contains it.
pub fn find_data_file() -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("NRSUR7DQ4_DATA") {
        let path = PathBuf::from(path);
        if path.is_file() {
            return Some(path);
        }
    }
    let dirs = std::env::var_os("LAL_DATA_PATH")?;
    std::env::split_paths(&dirs)
        .map(|dir| dir.join(NRSUR7DQ4_DATAFILE))
        .find(|path| path.is_file())
}

static GLOBAL_DATA: OnceCell<NRSur7dq4Data> = OnceCell::new();

impl NRSur7dq4Data {
    /// The process-wide surrogate data, loaded on first use from the file
    /// found by [`find_data_file`]. A failed load is retried on the next
    /// call.
    pub fn global() -> Result<&'static Self, NRSurError> {
        GLOBAL_DATA.get_or_try_init(|| {
            let path = find_data_file().ok_or_else(|| {
                NRSurError::DataFileNotFound(format!(
                    "unable to find {NRSUR7DQ4_DATAFILE}: set NRSUR7DQ4_DATA to the file or \
                     LAL_DATA_PATH to its directory (the data are available at \
                     https://zenodo.org/records/14999310)"
                ))
            })?;
            Self::load(path)
        })
    }
}

/// Errors from NRSur7dq4, following the failures of LALSimulation.
#[derive(Debug)]
pub enum NRSurError {
    /// The data file could not be located.
    DataFileNotFound(String),
    /// The data file could not be read.
    Hdf5(crate::hdf5::Hdf5Error),
    /// The data file does not have the expected content.
    InvalidData(String),
    /// A mass is not positive and finite.
    InvalidMass,
    /// The sampling interval is not positive and finite.
    InvalidTimeStep(f64),
    /// The mass ratio `m1 / m2` (heavier over lighter) is below 1
    /// (`XLAL_EINVAL`).
    InvalidMassRatio(f64),
    /// The mass ratio exceeds the extrapolation limit 6.01 (`XLAL_EDOM`).
    MassRatioTooLarge { q: f64, q_max: f64 },
    /// A spin magnitude exceeds 1 (`XLAL_EDOM`).
    SpinTooLarge {
        label: &'static str,
        magnitude: f64,
        chi_max: f64,
    },
    /// The reference frequency is below the starting frequency
    /// (`get_dimless_omega`, `XLAL_EINVAL`).
    ReferenceBelowMinimum { f_ref: f64, f_min: f64 },
    /// The dimensionless orbital reference frequency exceeds 0.201
    /// (`XLAL_EINVAL`).
    ReferenceFrequencyTooHigh { omega: f64 },
    /// The orbital frequency is below that of the first surrogate node for
    /// this binary (`XLAL_EDOM`).
    FrequencyBelowSurrogate { omega: f64, omega_min: f64 },
    /// `f_min` is neither zero nor above the surrogate's starting frequency
    /// for this binary (`XLAL_EDOM`).
    StartFrequency { f_min: f64, f_start: f64 },
    /// The mode array contains modes the model does not have.
    InvalidModeArray,
    /// Other failures of the surrogate evaluation (`XLAL_FAILURE`).
    Failure(String),
}

impl fmt::Display for NRSurError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DataFileNotFound(why) => write!(f, "{why}"),
            Self::Hdf5(err) => write!(f, "reading the surrogate data: {err}"),
            Self::InvalidData(why) => write!(f, "invalid surrogate data: {why}"),
            Self::InvalidMass => write!(f, "masses must be positive and finite"),
            Self::InvalidTimeStep(dt) => write!(f, "invalid sampling interval {dt}"),
            Self::InvalidMassRatio(q) => write!(f, "invalid mass ratio q = {q:.4} < 1"),
            Self::MassRatioTooLarge { q, q_max } => write!(
                f,
                "too much extrapolation: mass ratio q = {q:.4} > {q_max:.4}, the maximum \
                 allowed value"
            ),
            Self::SpinTooLarge {
                label,
                magnitude,
                chi_max,
            } => write!(
                f,
                "too much extrapolation: spin magnitude |chi{label}| = {magnitude:.4} > \
                 {chi_max:.4}, the maximum allowed value"
            ),
            Self::ReferenceBelowMinimum { f_ref, f_min } => {
                write!(f, "fRef ({f_ref}) cannot be lesser than fMin ({f_min})")
            }
            Self::ReferenceFrequencyTooHigh { omega } => {
                write!(
                    f,
                    "reference frequency omega_ref={omega:.4} > 0.2, too large"
                )
            }
            Self::FrequencyBelowSurrogate { omega, omega_min } => write!(
                f,
                "got omega_ref or omega_low={omega:.4} smaller than the minimum \
                 omega={omega_min:.4} for this configuration"
            ),
            Self::StartFrequency { f_min, f_start } => write!(
                f,
                "fMin should be 0 or >= {f_start:.8} for this configuration, got {f_min:.8}"
            ),
            Self::InvalidModeArray => {
                write!(f, "a requested mode is not available (ell must be 2..=4)")
            }
            Self::Failure(why) => write!(f, "NRSur7dq4 failed: {why}"),
        }
    }
}

impl std::error::Error for NRSurError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Hdf5(err) => Some(err),
            _ => None,
        }
    }
}

impl From<crate::hdf5::Hdf5Error> for NRSurError {
    fn from(err: crate::hdf5::Hdf5Error) -> Self {
        Self::Hdf5(err)
    }
}
