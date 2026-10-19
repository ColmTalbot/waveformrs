//! A common interface to the models of the crate, in the conventions of
//! LALSimulation's `SimInspiralChooseFDWaveform` and
//! `SimInspiralChooseTDWaveform`.
//!
//! * [`Binary`] holds the parameters of a binary as LAL's `ChooseFDWaveform`
//!   takes them (component masses, spin vectors at the reference frequency,
//!   distance, inclination, orbital reference phase and reference
//!   frequency).
//! * [`FrequencyDomainModel`] is implemented by [`TaylorF2`],
//!   [`IMRPhenomD`], [`IMRPhenomXAS`], [`IMRPhenomXHM`],
//!   [`IMRPhenomXP`](crate::imrphenomx::IMRPhenomXP),
//!   [`IMRPhenomXPHM`](crate::imrphenomx::IMRPhenomXPHM) and
//!   [`IMRPhenomXPNR`](crate::imrphenomx::IMRPhenomXPNR);
//!   [`TimeDomainModel`] by [`NRSur7dq4`]. All
//!   parameters are fixed when the model is built with `from_binary`, so
//!   the models are interchangeable.
//! * [`Approximant`] names the models as LAL does, and [`choose_fd_model`] /
//!   [`choose_td_model`] build a boxed model from an approximant, a binary
//!   and [`WaveformOptions`].
//! * Errors are [`crate::Error`].
//!
//! ```
//! use waveformrs::{choose_fd_model, Approximant, Binary, FrequencyDomainModel, WaveformOptions};
//!
//! let binary = Binary {
//!     mass_1: 36.0,
//!     mass_2: 29.0,
//!     spin_1: [0.4, 0.2, 0.3],
//!     spin_2: [-0.2, 0.3, -0.1],
//!     luminosity_distance: 500.0,
//!     inclination: 0.7,
//!     phi_ref: 1.2,
//!     f_ref: 20.0,
//! };
//! let frequencies: Vec<f64> = (0..100).map(|i| 20.0 + 2.0 * i as f64).collect();
//! let model = choose_fd_model(Approximant::IMRPhenomXPHM, &binary, &WaveformOptions::default())?;
//! let (hplus, hcross) = model.polarizations(&frequencies)?;
//! # assert_eq!(hplus.len(), 100);
//! # Ok::<(), waveformrs::Error>(())
//! ```

mod frequency_domain;
mod modes;
mod time_domain;

use std::fmt;
use std::str::FromStr;

use num_complex::Complex;

use crate::error::{Error, Result};
use crate::grid::FrequencyGrid;
use crate::imrphenomd::IMRPhenomD;
use crate::imrphenomx::{IMRPhenomXAS, IMRPhenomXHM, Polarizations, PrecessionVersion};
use crate::nrsur::{NRSur7dq4, NRSurWorkspace};
use crate::taylor::TaylorF2;

pub use crate::nrsur::{Mode as TimeSeriesMode, TimeSeriesModes, TimeSeriesPolarizations};
pub use modes::{FrequencyMode, FrequencyModes};

/// The parameters of a compact binary in the conventions of
/// `SimInspiralChooseFDWaveform` (also used by `SimInspiralChooseTDWaveform`
/// for NRSur7dq4). [`crate::imrphenomx::PrecessingBinary`] is the same
/// type.
///
/// The spins are the dimensionless spin vectors in the source frame at the
/// reference frequency, where the Newtonian orbital angular momentum `L` is
/// along `z` and the separation vector from the second body to the first
/// along `x`. The line of sight has polar angle `inclination` and azimuth
/// `pi / 2 - phi_ref` in this frame.
///
/// The aligned-spin models use only the `z` components, and (like LAL)
/// refuse non-zero `x` or `y` components.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Binary {
    /// Mass of the first body in solar masses.
    pub mass_1: f64,
    /// Mass of the second body in solar masses.
    pub mass_2: f64,
    /// Dimensionless spin of the first body.
    pub spin_1: [f64; 3],
    /// Dimensionless spin of the second body.
    pub spin_2: [f64; 3],
    /// Luminosity distance in Mpc.
    pub luminosity_distance: f64,
    /// Angle between the orbital angular momentum and the line of sight at
    /// the reference frequency (rad).
    pub inclination: f64,
    /// Orbital reference phase (rad), LAL's `phiRef`.
    pub phi_ref: f64,
    /// Reference gravitational-wave frequency in Hz, LAL's `f_ref`. Zero
    /// selects the first (minimum) frequency of the grid for the IMRPhenom
    /// models (and `f_min` for NRSur7dq4), and coalescence for TaylorF2, as
    /// in LALSimulation.
    pub f_ref: f64,
}

impl Binary {
    /// The total mass in solar masses.
    pub fn total_mass(&self) -> f64 {
        self.mass_1 + self.mass_2
    }

    /// Whether all in-plane (`x`, `y`) spin components are exactly zero, as
    /// LAL requires for aligned-spin approximants.
    pub fn has_aligned_spins(&self) -> bool {
        self.spin_1[0] == 0.0
            && self.spin_1[1] == 0.0
            && self.spin_2[0] == 0.0
            && self.spin_2[1] == 0.0
    }
}

/// The inclination and orbital reference phase the aligned-spin models use
/// in [`FrequencyDomainModel`] (zero for models built with their own
/// constructors).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Extrinsic {
    pub(crate) inclination: f64,
    pub(crate) phi_ref: f64,
}

impl Extrinsic {
    pub(crate) fn of(binary: &Binary) -> Self {
        Self {
            inclination: binary.inclination,
            phi_ref: binary.phi_ref,
        }
    }
}

/// The waveform models, named as LALSimulation's `Approximant`s.
///
/// [`Display`](fmt::Display) gives LAL's name and [`FromStr`] parses it
/// (case-sensitive, as `SimInspiralGetApproximantFromString`).
///
/// ```
/// use waveformrs::Approximant;
///
/// let approximant: Approximant = "IMRPhenomXPHM".parse()?;
/// assert_eq!(approximant, Approximant::IMRPhenomXPHM);
/// assert_eq!(approximant.to_string(), "IMRPhenomXPHM");
/// # Ok::<(), waveformrs::Error>(())
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
#[allow(clippy::upper_case_acronyms)]
pub enum Approximant {
    /// [`TaylorF2`] (frequency domain, aligned spins, (2,2) mode).
    TaylorF2,
    /// [`IMRPhenomD`] (frequency domain, aligned spins, (2,2) mode).
    IMRPhenomD,
    /// [`IMRPhenomXAS`] (frequency domain, aligned spins, (2,2) mode).
    IMRPhenomXAS,
    /// [`IMRPhenomXHM`] (frequency domain, aligned spins, higher modes).
    IMRPhenomXHM,
    /// [`IMRPhenomXP`](crate::imrphenomx::IMRPhenomXP) (frequency domain,
    /// precessing, (2,2) mode).
    IMRPhenomXP,
    /// [`IMRPhenomXPHM`](crate::imrphenomx::IMRPhenomXPHM) (frequency
    /// domain, precessing, higher modes).
    IMRPhenomXPHM,
    /// [`IMRPhenomXPNR`](crate::imrphenomx::IMRPhenomXPNR) (frequency
    /// domain, precessing, higher modes, NR-tuned angles and (2,2) mode).
    IMRPhenomXPNR,
    /// [`NRSur7dq4`] (time domain, precessing, `l <= 4`).
    NRSur7dq4,
}

impl Approximant {
    /// Every approximant of the crate.
    pub const ALL: [Approximant; 8] = [
        Self::TaylorF2,
        Self::IMRPhenomD,
        Self::IMRPhenomXAS,
        Self::IMRPhenomXHM,
        Self::IMRPhenomXP,
        Self::IMRPhenomXPHM,
        Self::IMRPhenomXPNR,
        Self::NRSur7dq4,
    ];

    /// LAL's name of the approximant.
    pub fn name(&self) -> &'static str {
        match self {
            Self::TaylorF2 => "TaylorF2",
            Self::IMRPhenomD => "IMRPhenomD",
            Self::IMRPhenomXAS => "IMRPhenomXAS",
            Self::IMRPhenomXHM => "IMRPhenomXHM",
            Self::IMRPhenomXP => "IMRPhenomXP",
            Self::IMRPhenomXPHM => "IMRPhenomXPHM",
            Self::IMRPhenomXPNR => "IMRPhenomXPNR",
            Self::NRSur7dq4 => "NRSur7dq4",
        }
    }

    /// Whether the model is a [`FrequencyDomainModel`] (otherwise it is a
    /// [`TimeDomainModel`]).
    pub fn is_frequency_domain(&self) -> bool {
        !matches!(self, Self::NRSur7dq4)
    }

    /// Whether the model accepts in-plane spins.
    pub fn is_precessing(&self) -> bool {
        matches!(
            self,
            Self::IMRPhenomXP | Self::IMRPhenomXPHM | Self::IMRPhenomXPNR | Self::NRSur7dq4
        )
    }
}

impl fmt::Display for Approximant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for Approximant {
    type Err = Error;

    fn from_str(name: &str) -> Result<Self> {
        Self::ALL
            .into_iter()
            .find(|approximant| approximant.name() == name)
            .ok_or_else(|| Error::UnknownApproximant(name.to_string()))
    }
}

/// Options of [`choose_fd_model`], [`choose_td_model`] and
/// [`choose_fd_modes`], the settings LAL takes from its `LALDict` (and the
/// switch for the automatic changes of the precession prescription, which
/// LAL does not have).
///
/// ```
/// use waveformrs::{imrphenomx::PrecessionVersion, WaveformOptions};
///
/// let options = WaveformOptions::default()
///     .with_precession_version(PrecessionVersion::Msa)
///     .with_mode_array(&[(2, 2), (2, -2), (3, 3), (3, -3)]);
/// # assert!(options.mode_array.is_some());
/// ```
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct WaveformOptions {
    /// The Euler angles of IMRPhenomXP and IMRPhenomXPHM (LAL's
    /// `PhenomXPrecVersion`). The default is the crate's default,
    /// [`PrecessionVersion::SpinTaylor`] (320); note that LALSimulation 6.2.0
    /// uses [`PrecessionVersion::Msa`] (223) when no version is given (see
    /// [`WaveformOptions::lalsimulation_defaults`]). Ignored by the other
    /// models (IMRPhenomXPNR always starts from the SpinTaylor angles
    /// connected to the PNR angles, LAL's `PhenomXPrecVersion` 330).
    pub precession_version: PrecessionVersion,
    /// LAL's `ModeArray`: the modes `(l, m)` of IMRPhenomXHM (with the sign
    /// of `m`), the co-precessing modes of IMRPhenomXPHM (`m` of either sign
    /// selects `|m|`, as in LAL) or the coorbital modes of NRSur7dq4. `None`
    /// selects all modes of the model. The (2,2)-mode models accept no mode
    /// array. IMRPhenomXPNR takes it as IMRPhenomXPHM.
    pub mode_array: Option<Vec<(u32, i32)>>,
    /// Whether IMRPhenomXP and IMRPhenomXPHM may change the Euler-angle
    /// prescription automatically as LALSimulation does (`true`, the
    /// default): SpinTaylor to NNLO for in-plane spins below `1e-7`,
    /// SpinTaylor to MSA when the SpinTaylor angles cannot be computed, MSA
    /// to NNLO when the MSA system cannot be initialised. With `false` the
    /// requested prescription is used or the evaluation fails with
    /// [`crate::imrphenomx::IMRPhenomXError::PrecessionVersionUnavailable`].
    /// For IMRPhenomXPNR, `false` requires the PNR angles with the
    /// SpinTaylor angles (see
    /// [`IMRPhenomXPNR::with_version_fallback`](crate::imrphenomx::IMRPhenomXPNR::with_version_fallback)).
    pub version_fallback: bool,
    /// IMRPhenomXPNR: whether to apply the fix of lalsuite merge request
    /// !2593 to the phase of the antisymmetric (2,2) mode (`true`, the
    /// default, as LALSimulation master after 6.2.1), or to reproduce
    /// LALSimulation 6.2.0 and 6.2.1 (`false`; see
    /// [`IMRPhenomXPNR::with_parity_fix`](crate::imrphenomx::IMRPhenomXPNR::with_parity_fix)).
    /// Ignored by the other models.
    pub pnr_parity_fix: bool,
}

impl Default for WaveformOptions {
    fn default() -> Self {
        Self {
            precession_version: PrecessionVersion::default(),
            mode_array: None,
            version_fallback: true,
            pnr_parity_fix: true,
        }
    }
}

impl WaveformOptions {
    /// The defaults of LALSimulation 6.2.0 when the `LALDict` is empty:
    /// the MSA Euler angles ([`PrecessionVersion::Msa`], LAL's default
    /// `PhenomXPrecVersion` 300 → 223), all modes and IMRPhenomXPNR without
    /// the fix of lalsuite !2593 ([`WaveformOptions::pnr_parity_fix`]).
    /// (LAL's multibanding, on by default, is not implemented: the models
    /// reproduce LAL with `PhenomXHMThresholdMband = PhenomXPHMThresholdMband
    /// = 0`.)
    pub fn lalsimulation_defaults() -> Self {
        Self {
            precession_version: PrecessionVersion::Msa,
            pnr_parity_fix: false,
            ..Self::default()
        }
    }

    /// Set [`WaveformOptions::precession_version`].
    pub fn with_precession_version(mut self, version: PrecessionVersion) -> Self {
        self.precession_version = version;
        self
    }

    /// Set [`WaveformOptions::mode_array`].
    pub fn with_mode_array(mut self, modes: &[(u32, i32)]) -> Self {
        self.mode_array = Some(modes.to_vec());
        self
    }

    /// Set [`WaveformOptions::version_fallback`].
    pub fn with_version_fallback(mut self, allow: bool) -> Self {
        self.version_fallback = allow;
        self
    }

    /// Set [`WaveformOptions::pnr_parity_fix`].
    pub fn with_pnr_parity_fix(mut self, fix: bool) -> Self {
        self.pnr_parity_fix = fix;
        self
    }
}

/// A frequency-domain waveform model with all parameters fixed.
///
/// The models reproduce `SimInspiralChooseFDWaveformSequence(phi_ref, m1,
/// m2, *spin_1, *spin_2, f_ref, distance, inclination, params, approximant,
/// frequencies)` of LALSimulation 6.2.0, with multibanding switched off for
/// IMRPhenomXHM and IMRPhenomXPHM (`PhenomXHMThresholdMband =
/// PhenomXPHMThresholdMband = 0`). In particular:
///
/// * `phi_ref` is the orbital reference phase for every model;
/// * `f_ref = 0` references the phase to the first frequency of each grid
///   (IMRPhenomD, IMRPhenomX*) or to coalescence (TaylorF2);
/// * the aligned-spin models reject non-zero in-plane spins
///   ([`Error::InPlaneSpins`]), as `ChooseFDWaveform` does.
///
/// The polarizations are `(h_+, h_x)` in strain per hertz at the given
/// frequencies (Hz, positive and increasing).
///
/// Several models have inherent methods with the same names and other
/// arguments (for example `IMRPhenomXHM::polarizations(f, phi_ref,
/// inclination)` and `IMRPhenomXPHM::modes()`, the list of co-precessing
/// modes). On a value of the concrete type the inherent method takes
/// precedence: call the trait methods through a generic parameter, a
/// `dyn FrequencyDomainModel`, or as `FrequencyDomainModel::polarizations(&model,
/// f)`.
///
/// The models are `Send + Sync`, so one model can be evaluated from several
/// threads (see [`crate::batch`]).
pub trait FrequencyDomainModel: Send + Sync + fmt::Debug {
    /// The model for `binary`, with LAL's default settings for the
    /// approximant (and the crate's default Euler angles for the precessing
    /// models, see [`WaveformOptions::precession_version`]). Returns an
    /// error where `SimInspiralChooseFDWaveform` would.
    fn from_binary(binary: &Binary) -> Result<Self>
    where
        Self: Sized;

    /// The approximant.
    fn approximant(&self) -> Approximant;

    /// `(h_+, h_x)` at `frequencies` (Hz).
    fn polarizations(&self, frequencies: &[f64]) -> Result<Polarizations> {
        let zero = Complex::new(0.0, 0.0);
        let mut hplus = vec![zero; frequencies.len()];
        let mut hcross = vec![zero; frequencies.len()];
        self.polarizations_into(frequencies, &mut hplus, &mut hcross)?;
        Ok((hplus, hcross))
    }

    /// [`FrequencyDomainModel::polarizations`] written into `hplus` and
    /// `hcross`.
    ///
    /// # Panics
    ///
    /// If `hplus` or `hcross` do not have the length of `frequencies`.
    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()>;

    /// [`FrequencyDomainModel::polarizations`] at the frequencies of a
    /// prepared [`FrequencyGrid`], written into `hplus` and `hcross`: the
    /// fast paths for many waveforms on the same frequencies, which agree
    /// with [`FrequencyDomainModel::polarizations`] to rounding (~1e-13
    /// relative) rather than bit for bit.
    ///
    /// # Panics
    ///
    /// If `hplus` or `hcross` do not have the length of the grid.
    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()>;

    /// The spherical-harmonic modes `h_lm` (without the `Y_lm`) at
    /// `frequencies` (Hz, positive and increasing), at `+f` and `-f`, in the
    /// conventions of `SimInspiralChooseFDModes` (see [`FrequencyModes`]):
    ///
    /// * IMRPhenomXHM: the modes of its mode array (by default `(l, +-m)`
    ///   for (2,2), (2,1), (3,3), (3,2), (4,4)), each as
    ///   `XLALSimIMRPhenomXHMFrequencySequenceOneMode`, including the
    ///   reference phase; `m < 0` at positive and `m > 0` at negative
    ///   frequencies.
    /// * IMRPhenomXPHM: the modes `(l, m)`, `|m| <= l`, in the inertial J
    ///   frame for every `l` of its co-precessing modes, each as
    ///   `XLALSimIMRPhenomXPHMFrequencySequenceOneMode` (with LAL's frequency
    ///   range of the SpinTaylor angles for each mode). Unlike LAL, a
    ///   co-precessing mode array without some `l` leaves out the J-frame
    ///   modes of that `l` instead of failing.
    /// * IMRPhenomXAS, IMRPhenomD, TaylorF2 (not supported by LAL's
    ///   `ChooseFDModes`): the modes (2,-2) and (2,2) in the conventions of
    ///   IMRPhenomXHM, with `h_{2,-2}(f) = h(f) / (2 sqrt(5 / 64 pi) e^{-i
    ///   pi})` for the face-on `h_+ = h`; for IMRPhenomXAS this is the (2,2)
    ///   mode of IMRPhenomXHM.
    /// * IMRPhenomXPNR: as IMRPhenomXPHM with the flags of IMRPhenomXPNR
    ///   (the PNR angles, the antisymmetric part of the (2,2) mode), each
    ///   J-frame mode with its own set-up as in LAL.
    /// * IMRPhenomXP (not supported by LAL): the J-frame modes `(2, m)` of
    ///   its (2,2) co-precessing mode, twisted up as for its polarizations.
    ///
    /// [`FrequencyModes::polarizations`] recovers the polarizations from the
    /// modes, to rounding (for IMRPhenomXPHM with the SpinTaylor angles, to
    /// the accuracy of the angles, whose integration range differs between
    /// the modes and the polarizations as in LAL).
    fn modes(&self, frequencies: &[f64]) -> Result<FrequencyModes>;
}

/// A time-domain waveform model with all parameters fixed (NRSur7dq4),
/// reproducing `SimInspiralChooseTDWaveform` and `SimInspiralChooseTDModes`.
pub trait TimeDomainModel: Send + Sync + fmt::Debug {
    /// Reusable buffers for repeated evaluation (one per thread).
    type Workspace: Default + Send;

    /// The model for `binary`. Returns an error where LAL would.
    fn from_binary(binary: &Binary) -> Result<Self>
    where
        Self: Sized;

    /// The approximant.
    fn approximant(&self) -> Approximant;

    /// `h_+` and `h_x` sampled every `delta_t` seconds from the (2,2)
    /// gravitational-wave frequency `f_min` (Hz; zero for the start of the
    /// model), as `SimInspiralChooseTDWaveform`.
    fn polarizations(&self, delta_t: f64, f_min: f64) -> Result<TimeSeriesPolarizations> {
        let mut out = TimeSeriesPolarizations::default();
        self.polarizations_into(&mut Self::Workspace::default(), delta_t, f_min, &mut out)?;
        Ok(out)
    }

    /// [`TimeDomainModel::polarizations`] reusing the buffers of
    /// `workspace` and `out`.
    fn polarizations_into(
        &self,
        workspace: &mut Self::Workspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesPolarizations,
    ) -> Result<()>;

    /// The inertial-frame modes `h_lm` (without the spherical harmonics),
    /// as `SimInspiralChooseTDModes`.
    fn modes(&self, delta_t: f64, f_min: f64) -> Result<TimeSeriesModes> {
        let mut out = TimeSeriesModes::default();
        self.modes_into(&mut Self::Workspace::default(), delta_t, f_min, &mut out)?;
        Ok(out)
    }

    /// [`TimeDomainModel::modes`] reusing the buffers of `workspace` and
    /// `out`.
    fn modes_into(
        &self,
        workspace: &mut Self::Workspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesModes,
    ) -> Result<()>;
}

impl<T: FrequencyDomainModel + ?Sized> FrequencyDomainModel for Box<T> {
    /// A boxed model cannot be built from a binary alone: use
    /// [`choose_fd_model`] (this returns [`Error::InvalidRequest`]).
    fn from_binary(_binary: &Binary) -> Result<Self> {
        Err(Error::InvalidRequest(
            "a boxed model cannot be built with from_binary: use choose_fd_model",
        ))
    }

    fn approximant(&self) -> Approximant {
        (**self).approximant()
    }

    fn polarizations(&self, frequencies: &[f64]) -> Result<Polarizations> {
        (**self).polarizations(frequencies)
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        (**self).polarizations_into(frequencies, hplus, hcross)
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        (**self).polarizations_on_grid(grid, hplus, hcross)
    }

    fn modes(&self, frequencies: &[f64]) -> Result<FrequencyModes> {
        (**self).modes(frequencies)
    }
}

/// The frequency-domain model `approximant` for `binary`, the analogue of
/// `SimInspiralChooseFDWaveform` (see [`FrequencyDomainModel`] for the
/// conventions): the settings of the `LALDict` are in `options`.
///
/// Returns [`Error::Unsupported`] for time-domain approximants and for
/// options the approximant does not take (a mode array for the (2,2)-mode
/// models), and the errors of the model's constructor.
pub fn choose_fd_model(
    approximant: Approximant,
    binary: &Binary,
    options: &WaveformOptions,
) -> Result<Box<dyn FrequencyDomainModel>> {
    let no_mode_array = || match options.mode_array {
        Some(_) => Err(Error::Unsupported {
            approximant,
            what: "a mode array",
        }),
        None => Ok(()),
    };
    Ok(match approximant {
        Approximant::TaylorF2 => {
            no_mode_array()?;
            Box::new(TaylorF2::from_binary(binary)?)
        }
        Approximant::IMRPhenomD => {
            no_mode_array()?;
            Box::new(IMRPhenomD::from_binary(binary)?)
        }
        Approximant::IMRPhenomXAS => {
            no_mode_array()?;
            Box::new(IMRPhenomXAS::from_binary(binary)?)
        }
        Approximant::IMRPhenomXHM => {
            let mut model = IMRPhenomXHM::from_binary(binary)?;
            if let Some(modes) = &options.mode_array {
                model = frequency_domain::xhm_with_mode_array(model, modes)?;
            }
            Box::new(model)
        }
        Approximant::IMRPhenomXP => {
            no_mode_array()?;
            Box::new(frequency_domain::xp(binary, options)?)
        }
        Approximant::IMRPhenomXPHM => Box::new(frequency_domain::xphm(binary, options)?),
        Approximant::IMRPhenomXPNR => Box::new(frequency_domain::xpnr(binary, options)?),
        Approximant::NRSur7dq4 => {
            return Err(Error::Unsupported {
                approximant,
                what: "frequency-domain waveforms (it is a time-domain model)",
            })
        }
    })
}

/// The time-domain model `approximant` for `binary`, the analogue of
/// `SimInspiralChooseTDWaveform`: NRSur7dq4 with the process-wide surrogate
/// data ([`crate::nrsur::NRSur7dq4Data::global`]) and the coorbital mode
/// array of `options`.
///
/// Returns [`Error::Unsupported`] for frequency-domain approximants.
pub fn choose_td_model(
    approximant: Approximant,
    binary: &Binary,
    options: &WaveformOptions,
) -> Result<Box<dyn TimeDomainModel<Workspace = NRSurWorkspace>>> {
    match approximant {
        Approximant::NRSur7dq4 => {
            let mut model = NRSur7dq4::from_binary(binary)?;
            if let Some(modes) = &options.mode_array {
                model = model.with_mode_array(modes)?;
            }
            Ok(Box::new(model))
        }
        _ => Err(Error::Unsupported {
            approximant,
            what: "time-domain waveforms (it is a frequency-domain model)",
        }),
    }
}

/// `SimInspiralChooseFDModes(m1, m2, *spin_1, *spin_2, delta_f, f_min,
/// f_max, f_ref, phi_ref, distance, inclination, params, approximant)` of
/// LALSimulation 6.2.0 for the approximants it supports, IMRPhenomXHM,
/// IMRPhenomXPHM and IMRPhenomXPNR (with `PhenomXHMThresholdMband =
/// PhenomXPHMThresholdMband = 0`).
///
/// The modes are on the uniform grid `f_k = k delta_f`, `k = 0 .. n - 1`,
/// with `n - 1` the next power of two of `max(f_max, f_cut) / delta_f`
/// (`f_max = 0` selects the cutoff `f_cut` at `Mf = 0.3`), and vanish below
/// `f_min` and above `min(f_max, f_cut)`. LAL's two-sided series of mode
/// `(l, m)` at `(i - n + 1) delta_f`, `i = 0 .. 2n - 2`, is
/// `negative[n - 1 - i]` for `i < n - 1` and `positive[i - n + 1]` after.
/// `f_ref = 0` selects `f_min`. The mode array of `options` is IMRPhenomXHM's
/// (with signs) or IMRPhenomXPHM's (IMRPhenomXPNR's) co-precessing one;
/// IMRPhenomXPHM and IMRPhenomXPNR return all J-frame modes of the `l` they
/// have (LAL's default `ModeArrayJframe`).
/// See [`FrequencyDomainModel::modes`] for the conventions.
pub fn choose_fd_modes(
    approximant: Approximant,
    binary: &Binary,
    delta_f: f64,
    f_min: f64,
    f_max: f64,
    options: &WaveformOptions,
) -> Result<FrequencyModes> {
    match approximant {
        Approximant::IMRPhenomXHM => {
            let mut model = IMRPhenomXHM::from_binary(binary)?;
            if let Some(modes) = &options.mode_array {
                model = frequency_domain::xhm_with_mode_array(model, modes)?;
            }
            Ok(model.frequency_modes_uniform(
                delta_f,
                f_min,
                f_max,
                binary.phi_ref,
                binary.inclination,
            )?)
        }
        Approximant::IMRPhenomXPHM => {
            Ok(frequency_domain::xphm(binary, options)?
                .jframe_modes_uniform(delta_f, f_min, f_max)?)
        }
        Approximant::IMRPhenomXPNR => {
            Ok(frequency_domain::xpnr(binary, options)?
                .jframe_modes_uniform(delta_f, f_min, f_max)?)
        }
        _ => Err(Error::Unsupported {
            approximant,
            what: "SimInspiralChooseFDModes (LAL supports IMRPhenomXHM, IMRPhenomXPHM and \
                   IMRPhenomXPNR)",
        }),
    }
}
