// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
// the unit tests that share tests/common refer to the crate by name
#[cfg(test)]
extern crate self as waveformrs;

/// The helpers of the integration tests, for the unit tests that compare
/// with the same reference data.
#[cfg(test)]
#[path = "../tests/common/mod.rs"]
pub(crate) mod test_common;

pub mod batch;
pub(crate) mod complex;
pub mod conditioning;
pub(crate) mod constants;
pub mod conversions;
pub mod error;
pub(crate) mod fastmath;
#[cfg(feature = "fft")]
pub mod fft;
pub mod grid;
pub mod hdf5;
pub mod imrphenomd;
pub mod imrphenomx;
pub mod interpolate;
pub mod model;
pub mod nrsur;
pub mod nrtidal;
pub mod pn;
#[cfg(test)]
mod reference;
pub mod series;
pub mod spintaylor;
pub mod taylor;
pub mod waveform;

pub use error::{Error, Result};
pub use model::{
    choose_fd_model, choose_fd_modes, choose_td_model, Approximant, Binary, FrequencyDomainModel,
    FrequencyMode, FrequencyModes, TimeDomainModel, WaveformOptions,
};
