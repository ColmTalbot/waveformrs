//! Orbit-averaged SpinTaylorT4 evolution of precessing binaries.
//!
//! A port of `XLALSimInspiralSpinTaylorPNEvolveOrbit` from LALSimulation
//! 6.2.0 for the `SpinTaylorT4` approximant, including the adaptive
//! Runge-Kutta-Fehlberg integrator of LAL (`XLALAdaptiveRungeKutta4Hermite`)
//! and GSL 2.8 (`gsl_odeiv_*`) that it relies on.  This is the PN evolution
//! used to compute the "SpinTaylor" Euler angles of IMRPhenomXPHM
//! (`PrecVersion` 310/311/320/321/330).
//!
//! The implementation follows the C code operation by operation, including
//! the fused multiply-adds that clang emits for the arm64 builds of LAL (see
//! [`fp`]), so that the output agrees with LALSimulation to rounding
//! precision and has exactly the same length.
//!
//! # Units and conventions
//!
//! - Masses are in kg (as in LAL); [`SpinTaylorT4::from_solar_masses`]
//!   multiplies by `MSUN_SI` for convenience.
//! - Frequencies are gravitational-wave frequencies in Hz and the time step is
//!   in seconds.
//! - Spins are dimensionless (`chi_i`), all vectors are given in the frame of
//!   the caller, as in LAL.
//! - PN orders use LAL's conventions: twice the PN order, `-1` for all
//!   available terms.  The defaults are those IMRPhenomXPHM enforces
//!   (`phaseO = 7`, `spinO = 6`, `tideO = 12`, `lscorr = 0`); note that for
//!   SpinTaylorT4 `-1` gives exactly the same evolution.
//!
//! # Example
//!
//! ```
//! use waveformrs::spintaylor::SpinTaylorT4;
//!
//! let binary = SpinTaylorT4::from_solar_masses(36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1]);
//! let orbit = binary.evolve_orbit(1.0 / 512.0, 20.0, 400.0).unwrap();
//! assert_eq!(orbit.v.len(), orbit.lnhat_z.len());
//! ```

pub mod coefficients;
pub mod dynamics;
pub mod fp;
pub mod integrator;

use std::f64::consts::PI;
use std::fmt;

use dynamics::{mass_in_seconds, SpinTaylorT4System};
pub use dynamics::{SpinTaylorT4Coefficients, NUM_VARIABLES};
use integrator::{adaptive_rk4_hermite, IntegratorSettings};

use crate::constants::MSUN_SI;

/// Absolute tolerance of the adaptive integrator, `LAL_ST4_ABSOLUTE_TOLERANCE`.
pub const ABSOLUTE_TOLERANCE: f64 = 1.0e-12;
/// Relative tolerance of the adaptive integrator, `LAL_ST4_RELATIVE_TOLERANCE`.
pub const RELATIVE_TOLERANCE: f64 = 1.0e-12;

/// Errors raised by the SpinTaylor evolution, mirroring LAL's `XLAL_EINVAL`
/// and `XLAL_EFUNC` failures.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpinTaylorError {
    InvalidPhaseOrder(i32),
    InvalidSpinOrder(i32),
    InvalidTidalOrder(i32),
    /// `f_start` must be positive.
    InvalidStartFrequency(f64),
    /// `f_end` must be non-negative.
    InvalidEndFrequency(f64),
    /// The integration interval and the time step have different signs.
    InvalidInterval,
}

impl fmt::Display for SpinTaylorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPhaseOrder(o) => write!(f, "invalid PN phase order {o}"),
            Self::InvalidSpinOrder(o) => write!(
                f,
                "invalid PN spin order {o} (orbit-averaged dynamics support -1..=6)"
            ),
            Self::InvalidTidalOrder(o) => write!(f, "invalid PN tidal order {o}"),
            Self::InvalidStartFrequency(x) => write!(f, "fStart = {x} must be > 0"),
            Self::InvalidEndFrequency(x) => write!(f, "fEnd = {x} must be >= 0"),
            Self::InvalidInterval => {
                write!(f, "(tend_in - tinit) and deltat must have the same sign")
            }
        }
    }
}

impl std::error::Error for SpinTaylorError {}

/// Why the integration stopped (`integrator->returncode` in LAL).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Termination {
    /// The requested end frequency was crossed
    /// (`LALSIMINSPIRAL_ST_TEST_FREQBOUND`).
    FrequencyBound,
    /// The PN energy stopped decreasing (`LALSIMINSPIRAL_ST_TEST_ENERGY`).
    Energy,
    /// `d^2 omega / dt^2` changed sign
    /// (`LALSIMINSPIRAL_ST_TEST_OMEGADOUBLEDOT`).
    OmegaDoubleDot,
    /// The orbital frequency became NaN (`LALSIMINSPIRAL_ST_TEST_OMEGANAN`).
    OmegaNaN,
    /// `v >= 1` (`LALSIMINSPIRAL_ST_TEST_LARGEV`).
    LargeV,
    /// The derivative failed with a non-positive orbital frequency after
    /// all retries (`LALSIMINSPIRAL_ST_DERIVATIVE_OMEGANONPOS`).
    OmegaNonPositive,
    /// Any other code.
    Other(i32),
}

impl Termination {
    pub fn from_code(code: i32) -> Self {
        match code {
            dynamics::ST_TEST_FREQBOUND => Self::FrequencyBound,
            dynamics::ST_TEST_ENERGY => Self::Energy,
            dynamics::ST_TEST_OMEGADOUBLEDOT => Self::OmegaDoubleDot,
            dynamics::ST_TEST_OMEGANAN => Self::OmegaNaN,
            dynamics::ST_TEST_LARGEV => Self::LargeV,
            dynamics::ST_DERIVATIVE_OMEGANONPOS => Self::OmegaNonPositive,
            other => Self::Other(other),
        }
    }

    /// LAL considers `0`, the energy test and the frequency bound normal
    /// terminations and warns for anything else.
    pub fn is_normal(&self) -> bool {
        matches!(self, Self::Energy | Self::FrequencyBound | Self::Other(0))
    }
}

/// The evolved orbit: the same 14 time series that
/// `XLALSimInspiralSpinTaylorPNEvolveOrbit` returns, sampled every `delta_t`
/// seconds and ordered forward in time (also when integrating backwards).
#[derive(Clone, Debug, Default)]
pub struct SpinTaylorOrbit {
    /// Post-Newtonian parameter `v = (pi M f_GW)^(1/3)`.
    pub v: Vec<f64>,
    /// Orbital phase, zero at `f_start`.
    pub phi: Vec<f64>,
    /// Dimensionless spin of the first body.
    pub s1_x: Vec<f64>,
    pub s1_y: Vec<f64>,
    pub s1_z: Vec<f64>,
    /// Dimensionless spin of the second body.
    pub s2_x: Vec<f64>,
    pub s2_y: Vec<f64>,
    pub s2_z: Vec<f64>,
    /// Unit Newtonian orbital angular momentum.
    pub lnhat_x: Vec<f64>,
    pub lnhat_y: Vec<f64>,
    pub lnhat_z: Vec<f64>,
    /// Orbital-plane basis vector.
    pub e1_x: Vec<f64>,
    pub e1_y: Vec<f64>,
    pub e1_z: Vec<f64>,
    /// Sampling interval in seconds.
    pub delta_t: f64,
    /// Time of the first sample in seconds; the last sample is at `t = 0`.
    /// (LAL stores this as a GPS epoch, rounded to nanoseconds.)
    pub epoch: f64,
    /// LAL's integrator return code.
    pub return_code: i32,
}

impl SpinTaylorOrbit {
    pub fn len(&self) -> usize {
        self.v.len()
    }

    pub fn is_empty(&self) -> bool {
        self.v.is_empty()
    }

    pub fn termination(&self) -> Termination {
        Termination::from_code(self.return_code)
    }

    fn with_capacity(n: usize) -> Self {
        let v = || Vec::with_capacity(n);
        Self {
            v: v(),
            phi: v(),
            s1_x: v(),
            s1_y: v(),
            s1_z: v(),
            s2_x: v(),
            s2_y: v(),
            s2_z: v(),
            lnhat_x: v(),
            lnhat_y: v(),
            lnhat_z: v(),
            e1_x: v(),
            e1_y: v(),
            e1_z: v(),
            ..Default::default()
        }
    }
}

/// A binary to be evolved with the orbit-averaged SpinTaylorT4 equations.
#[derive(Clone, Copy, Debug)]
pub struct SpinTaylorT4 {
    /// Mass of the first body in kg.
    pub m1_si: f64,
    /// Mass of the second body in kg.
    pub m2_si: f64,
    /// Initial dimensionless spin of the first body.
    pub spin_1: [f64; 3],
    /// Initial dimensionless spin of the second body.
    pub spin_2: [f64; 3],
    /// Initial direction of the Newtonian orbital angular momentum.
    pub lnhat: [f64; 3],
    /// Initial orbital-plane basis vector.
    pub e1: [f64; 3],
    /// Tidal deformabilities `Lambda_i = lambda_i / m_i^5` (dimensionless).
    pub lambda_1: f64,
    pub lambda_2: f64,
    /// Quadrupole-monopole parameters (1 for black holes).
    pub quadparam_1: f64,
    pub quadparam_2: f64,
    /// Twice the PN order of the non-spinning phasing (`phaseO`).
    pub phase_order: i32,
    /// Twice the PN order of spin effects (`spinO`).
    pub spin_order: i32,
    /// Twice the PN order of tidal effects (`tideO`).
    pub tidal_order: i32,
    /// Include spin corrections to the orbital angular momentum (`lscorr`).
    pub lscorr: bool,
}

impl SpinTaylorT4 {
    /// A binary with masses in kg and dimensionless spins, with
    /// `LNhat = (0, 0, 1)`, `E1 = (1, 0, 0)`, black-hole matter parameters
    /// and the PN orders used by IMRPhenomXPHM.
    pub fn new(m1_si: f64, m2_si: f64, spin_1: [f64; 3], spin_2: [f64; 3]) -> Self {
        Self {
            m1_si,
            m2_si,
            spin_1,
            spin_2,
            lnhat: [0.0, 0.0, 1.0],
            e1: [1.0, 0.0, 0.0],
            lambda_1: 0.0,
            lambda_2: 0.0,
            quadparam_1: 1.0,
            quadparam_2: 1.0,
            phase_order: 7,
            spin_order: 6,
            tidal_order: 12,
            lscorr: false,
        }
    }

    /// As [`SpinTaylorT4::new`] with masses in solar masses.
    pub fn from_solar_masses(m1: f64, m2: f64, spin_1: [f64; 3], spin_2: [f64; 3]) -> Self {
        Self::new(m1 * MSUN_SI, m2 * MSUN_SI, spin_1, spin_2)
    }

    pub fn with_orientation(mut self, lnhat: [f64; 3], e1: [f64; 3]) -> Self {
        self.lnhat = lnhat;
        self.e1 = e1;
        self
    }

    pub fn with_matter(mut self, lambda: [f64; 2], quadparam: [f64; 2]) -> Self {
        [self.lambda_1, self.lambda_2] = lambda;
        [self.quadparam_1, self.quadparam_2] = quadparam;
        self
    }

    pub fn with_pn_orders(mut self, phase_order: i32, spin_order: i32, tidal_order: i32) -> Self {
        self.phase_order = phase_order;
        self.spin_order = spin_order;
        self.tidal_order = tidal_order;
        self
    }

    pub fn with_lscorr(mut self, lscorr: bool) -> Self {
        self.lscorr = lscorr;
        self
    }

    /// Set up the PN coefficients for an evolution from `f_start` to `f_end`.
    pub fn coefficients(
        &self,
        f_start: f64,
        f_end: f64,
    ) -> Result<SpinTaylorT4Coefficients, SpinTaylorError> {
        SpinTaylorT4Coefficients::new(
            self.m1_si,
            self.m2_si,
            f_start,
            f_end,
            self.lambda_1,
            self.lambda_2,
            self.quadparam_1,
            self.quadparam_2,
            self.spin_order,
            self.tidal_order,
            self.phase_order,
            self.lscorr,
        )
    }

    /// Port of `XLALSimInspiralSpinTaylorPNEvolveOrbit` for SpinTaylorT4.
    ///
    /// Evolves the binary from the GW frequency `f_start` (Hz) to `f_end`
    /// sampling every `delta_t` seconds.  If `0 < f_end < f_start` the
    /// integration runs backwards in time (the output is still ordered
    /// forward in time, ending at `f_start`); `f_end = 0` integrates forwards
    /// as far as possible.  As in LAL, a forward integration keeps the first
    /// sample beyond `f_end`, and a backward one the first sample below it.
    pub fn evolve_orbit(
        &self,
        delta_t: f64,
        f_start: f64,
        f_end: f64,
    ) -> Result<SpinTaylorOrbit, SpinTaylorError> {
        let coefficients = self.coefficients(f_start, f_end)?;
        self.evolve(&coefficients, delta_t)
    }

    fn evolve(
        &self,
        coefficients: &SpinTaylorT4Coefficients,
        delta_t: f64,
    ) -> Result<SpinTaylorOrbit, SpinTaylorError> {
        let f_start = coefficients.f_start;
        let f_end = coefficients.f_end;
        if f_start <= 0. || f_start.is_nan() {
            return Err(SpinTaylorError::InvalidStartFrequency(f_start));
        }
        if f_end < 0. {
            return Err(SpinTaylorError::InvalidEndFrequency(f_end));
        }
        let backwards = f_end < f_start && f_end != 0.;
        let sgn = if backwards { -1.0 } else { 1.0 };

        let m1sec = mass_in_seconds(self.m1_si);
        let m2sec = mass_in_seconds(self.m2_si);
        let msec = m1sec + m2sec;
        let mcsec = msec * (m1sec * m2sec / msec / msec).powf(0.6);

        // Estimate length of waveform using Newtonian t(f) formula
        let dt_start =
            (5.0 / 256.0) * PI.powf(-8.0 / 3.0) * (mcsec * f_start).powf(-5.0 / 3.0) / f_start;
        let dt_end = if f_end == 0. {
            0.
        } else {
            (5.0 / 256.0) * PI.powf(-8.0 / 3.0) * (mcsec * f_end).powf(-5.0 / 3.0) / f_end
        };
        let lengths = dt_start - dt_end;

        let norm1 = m1sec * m1sec / msec / msec;
        let norm2 = m2sec * m2sec / msec / msec;
        let (s1, s2) = (self.spin_1, self.spin_2);
        let mut yinit = [
            0.,
            PI * msec * f_start,
            self.lnhat[0],
            self.lnhat[1],
            self.lnhat[2],
            norm1 * s1[0],
            norm1 * s1[1],
            norm1 * s1[2],
            norm2 * s2[0],
            norm2 * s2[1],
            norm2 * s2[2],
            self.e1[0],
            self.e1[1],
            self.e1[2],
        ];

        let settings = IntegratorSettings {
            eps_abs: ABSOLUTE_TOLERANCE,
            eps_rel: RELATIVE_TOLERANCE,
            retries: 6,
            stop_on_test_only: true,
        };
        let mut system = SpinTaylorT4System {
            coefficients,
            prev_domega: 0.0,
        };
        // time is measured in units of the total mass
        let output = adaptive_rk4_hermite(
            &mut system,
            &settings,
            &mut yinit,
            0.0,
            lengths / msec,
            sgn * delta_t / msec,
        )
        .map_err(|_| SpinTaylorError::InvalidInterval)?;
        let samples = &output.y;
        let len = samples.len();

        // The integrator takes one adaptive stride past fEnd; keep the first
        // interpolated sample to cross it.
        let mut cutlen = len;
        if f_end != 0. && f_end < f_start {
            let w_end = PI * msec * f_end;
            while cutlen > 0 && samples[cutlen - 1][1] < w_end {
                cutlen -= 1;
            }
            if cutlen < len {
                cutlen += 1;
            }
        } else if f_end > f_start {
            let w_end = PI * msec * f_end;
            while cutlen > 0 && samples[cutlen - 1][1] > w_end {
                cutlen -= 1;
            }
            if cutlen < len {
                cutlen += 1;
            }
        }

        let mut orbit = SpinTaylorOrbit::with_capacity(cutlen);
        orbit.delta_t = delta_t;
        orbit.epoch = -1.0 * (cutlen as f64 - 1.0) * delta_t;
        orbit.return_code = output.return_code;
        let mut push = |y: &[f64; NUM_VARIABLES]| {
            orbit.phi.push(y[0]);
            orbit.v.push(y[1].cbrt());
            orbit.lnhat_x.push(y[2]);
            orbit.lnhat_y.push(y[3]);
            orbit.lnhat_z.push(y[4]);
            // spins are returned in the standard convention
            orbit.s1_x.push(y[5] / norm1);
            orbit.s1_y.push(y[6] / norm1);
            orbit.s1_z.push(y[7] / norm1);
            orbit.s2_x.push(y[8] / norm2);
            orbit.s2_y.push(y[9] / norm2);
            orbit.s2_z.push(y[10] / norm2);
            orbit.e1_x.push(y[11]);
            orbit.e1_y.push(y[12]);
            orbit.e1_z.push(y[13]);
        };
        if backwards {
            samples[..cutlen].iter().rev().for_each(&mut push);
        } else {
            samples[..cutlen].iter().for_each(&mut push);
        }
        Ok(orbit)
    }
}
