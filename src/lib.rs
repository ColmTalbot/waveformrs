pub mod batch;
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
pub mod pn;
pub mod series;
pub mod spintaylor;
pub mod taylor;
pub mod waveform;

pub use error::{Error, Result};
pub use model::{
    choose_fd_model, choose_fd_modes, choose_td_model, Approximant, Binary, FrequencyDomainModel,
    FrequencyMode, FrequencyModes, TimeDomainModel, WaveformOptions,
};
