// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Orbit-averaged SpinTaylorT4 evolution of precessing binaries.
//!
//! A port of `XLALSimInspiralSpinTaylorPNEvolveOrbit` from LALSimulation
//! 6.2.0 for the `SpinTaylorT4` approximant, including the adaptive
//! Runge-Kutta-Fehlberg integrator of LAL (`XLALAdaptiveRungeKutta4Hermite`)
//! and GSL 2.8 (`gsl_odeiv_*`) that it relies on.  This is the PN evolution
//! used to compute the "SpinTaylor" Euler angles of IMRPhenomXPHM
//! (`PrecVersion` 310/311/320/321/330).
//!
//! The implementation follows the C code operation by operation (with fused
//! multiply-adds where the hardware has them, see [`fp`]), so that the
//! output agrees with LALSimulation to rounding precision and, in all but
//! rare cases, has the same length.
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
use integrator::{adaptive_rk4_hermite_into, hermite_interpolate, HermiteSink, IntegratorSettings};

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

/// Storage for the samples of an evolution, in the order of the integration
/// (see [`SpinTaylorT4::evolve_into`]).
pub(crate) trait OrbitSamples: HermiteSink<NUM_VARIABLES> {
    /// Number of samples recorded.
    fn len(&self) -> usize;
    /// The (dimensionless) orbital frequency `M omega` of sample `i`.
    fn omega(&self, i: usize) -> f64;
}

/// Every sample with its full state.
#[derive(Clone, Debug, Default)]
pub(crate) struct FullSamples(pub(crate) Vec<[f64; NUM_VARIABLES]>);

impl HermiteSink<NUM_VARIABLES> for FullSamples {
    fn reserve(&mut self, n: usize) {
        self.0.reserve(n);
    }

    #[inline]
    fn push(&mut self, _t: f64, y: &[f64; NUM_VARIABLES]) {
        self.0.push(*y);
    }
}

impl OrbitSamples for FullSamples {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn omega(&self, i: usize) -> f64 {
        self.0[i][1]
    }
}

/// One accepted step of an integration, from `t0` to `t1` (in units of the
/// total mass), with what the interpolant of the step needs.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DenseStep {
    t0: f64,
    t1: f64,
    y0: [f64; NUM_VARIABLES],
    y1: [f64; NUM_VARIABLES],
    k1: [f64; NUM_VARIABLES],
    k6: [f64; NUM_VARIABLES],
}

impl DenseStep {
    fn state(&self, t: f64) -> [f64; NUM_VARIABLES] {
        let h = self.t1 - self.t0;
        hermite_interpolate((t - self.t0) / h, h, &self.y0, &self.y1, &self.k1, &self.k6)
    }

    fn contains(&self, t: f64) -> bool {
        (self.t0 <= t && t <= self.t1) || (self.t1 <= t && t <= self.t0)
    }
}

/// Every accepted step of an integration, and `M omega` of the evenly
/// spaced samples (from which the length of LAL's orbit follows).
#[derive(Debug, Default)]
struct DenseSamples {
    omega: Vec<f64>,
    steps: Vec<DenseStep>,
}

impl HermiteSink<NUM_VARIABLES> for DenseSamples {
    fn reserve(&mut self, n: usize) {
        self.omega.reserve(n);
    }

    #[inline]
    fn push(&mut self, _t: f64, y: &[f64; NUM_VARIABLES]) {
        self.omega.push(y[1]);
    }

    #[inline]
    fn step(
        &mut self,
        t0: f64,
        t1: f64,
        y0: &[f64; NUM_VARIABLES],
        y1: &[f64; NUM_VARIABLES],
        k1: &[f64; NUM_VARIABLES],
        k6: &[f64; NUM_VARIABLES],
    ) {
        self.steps.push(DenseStep {
            t0,
            t1,
            y0: *y0,
            y1: *y1,
            k1: *k1,
            k6: *k6,
        });
    }
}

impl OrbitSamples for DenseSamples {
    fn len(&self) -> usize {
        self.omega.len()
    }

    fn omega(&self, i: usize) -> f64 {
        self.omega[i]
    }
}

/// An integration with its solution at any time: the interpolant of each
/// accepted step (the one LAL uses for the evenly spaced samples). Time is
/// in units of the total mass from the start of the integration, negative
/// for an integration backwards in time.
#[derive(Debug)]
pub(crate) struct DenseOrbit {
    coefficients: SpinTaylorT4Coefficients,
    steps: Vec<DenseStep>,
    evolution: Evolution,
}

impl DenseOrbit {
    /// The times covered, earliest first.
    pub(crate) fn span(&self) -> (f64, f64) {
        let (first, last) = (self.steps[0].t0, self.steps[self.steps.len() - 1].t1);
        (first.min(last), first.max(last))
    }

    /// The state at `t`, within [`Self::span`].
    pub(crate) fn state(&self, t: f64) -> [f64; NUM_VARIABLES] {
        // the steps are ordered by time in the direction of the integration
        let forward = !self.evolution.backwards;
        let i = self
            .steps
            .partition_point(|step| if forward { step.t1 < t } else { step.t1 > t });
        let step = &self.steps[i.min(self.steps.len() - 1)];
        debug_assert!(step.contains(t) || i == self.steps.len());
        step.state(t)
    }

    /// The dimensionless spins of a state.
    pub(crate) fn spins(&self, y: &[f64; NUM_VARIABLES]) -> ([f64; 3], [f64; 3]) {
        (self.evolution.s1(y), self.evolution.s2(y))
    }

    /// `d^2 (M omega) / dt^2` along the solution at `t`, from the derivatives
    /// on either side of the state along its direction of motion.
    fn omega_double_dot(&self, t: f64) -> f64 {
        const EPSILON: f64 = 1e-3;
        let y = self.state(t);
        let mut dydt = [0.0; NUM_VARIABLES];
        self.coefficients.derivatives(&y, &mut dydt);
        let shifted = |sign: f64| {
            let z: [f64; NUM_VARIABLES] = std::array::from_fn(|k| y[k] + sign * EPSILON * dydt[k]);
            let mut dz = [0.0; NUM_VARIABLES];
            self.coefficients.derivatives(&z, &mut dz);
            dz[1]
        };
        (shifted(1.0) - shifted(-1.0)) / (2.0 * EPSILON)
    }

    /// The time at which a forward integration stops, found exactly rather
    /// than at the end of the step after which the stopping test fails:
    /// where the frequency reaches its bound, where the PN energy stops
    /// decreasing, where `v = 1`, or where `d omega / dt` is largest,
    /// whichever stopped the integration. The end of the last step for any
    /// other reason.
    pub(crate) fn end_time(&self) -> f64 {
        let n = self.steps.len();
        let last = &self.steps[n - 1];
        let (mut low, mut high) = (last.t0, last.t1);
        let condition: Box<dyn Fn(f64) -> f64 + '_> = match self.evolution.return_code {
            dynamics::ST_TEST_FREQBOUND => {
                let omega_end = self.coefficients.omega_end();
                Box::new(move |t| self.state(t)[1] - omega_end)
            }
            dynamics::ST_TEST_ENERGY => {
                Box::new(|t| -self.coefficients.energy_test(&self.state(t)))
            }
            dynamics::ST_TEST_LARGEV => Box::new(|t| self.state(t)[1] - 1.0),
            dynamics::ST_TEST_OMEGADOUBLEDOT => {
                // The test compares d omega / dt at the ends of the last two
                // steps.
                low = self.steps[n.saturating_sub(2)].t0;
                Box::new(|t| -self.omega_double_dot(t))
            }
            _ => return high,
        };
        // the condition holds at `high` and not at `low`: bisect
        if !(condition(low) < 0.0 && condition(high) >= 0.0) {
            return high;
        }
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if middle == low || middle == high {
                break;
            }
            if condition(middle) >= 0.0 {
                high = middle;
            } else {
                low = middle;
            }
        }
        high
    }
}

/// An evolution recorded in some [`OrbitSamples`]: which samples form the
/// orbit and how a state converts to the series of [`SpinTaylorOrbit`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct Evolution {
    /// Number of samples of the orbit: the first `len` recorded.
    pub(crate) len: usize,
    /// Whether the integration ran backwards in time, in which case the
    /// orbit (ordered forward in time) is the recorded samples reversed.
    pub(crate) backwards: bool,
    norm1: f64,
    norm2: f64,
    /// LAL's integrator return code.
    pub(crate) return_code: i32,
}

impl Evolution {
    /// Index among the recorded samples of sample `i` of the orbit.
    #[inline]
    pub(crate) fn recorded(&self, i: usize) -> usize {
        if self.backwards {
            self.len - 1 - i
        } else {
            i
        }
    }

    /// `v = (M omega)^(1/3)`.
    #[inline]
    pub(crate) fn v(omega: f64) -> f64 {
        omega.cbrt()
    }

    /// The unit Newtonian orbital angular momentum of a state.
    #[inline]
    pub(crate) fn lnhat(y: &[f64; NUM_VARIABLES]) -> [f64; 3] {
        [y[2], y[3], y[4]]
    }

    /// The dimensionless spin of the first body of a state.
    #[inline]
    pub(crate) fn s1(&self, y: &[f64; NUM_VARIABLES]) -> [f64; 3] {
        [y[5] / self.norm1, y[6] / self.norm1, y[7] / self.norm1]
    }

    /// The dimensionless spin of the second body of a state.
    #[inline]
    pub(crate) fn s2(&self, y: &[f64; NUM_VARIABLES]) -> [f64; 3] {
        [y[8] / self.norm2, y[9] / self.norm2, y[10] / self.norm2]
    }

    /// The orbital-plane basis vector of a state.
    #[inline]
    pub(crate) fn e1(y: &[f64; NUM_VARIABLES]) -> [f64; 3] {
        [y[11], y[12], y[13]]
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
        let mut samples = FullSamples::default();
        let evolution = self.evolve_into(&coefficients, delta_t, &mut samples)?;
        let samples = &samples.0;
        let cutlen = evolution.len;
        let mut orbit = SpinTaylorOrbit::with_capacity(cutlen);
        orbit.delta_t = delta_t;
        orbit.epoch = -1.0 * (cutlen as f64 - 1.0) * delta_t;
        orbit.return_code = evolution.return_code;
        let mut push = |y: &[f64; NUM_VARIABLES]| {
            orbit.phi.push(y[0]);
            orbit.v.push(Evolution::v(y[1]));
            let [x, y_, z] = Evolution::lnhat(y);
            orbit.lnhat_x.push(x);
            orbit.lnhat_y.push(y_);
            orbit.lnhat_z.push(z);
            // spins are returned in the standard convention
            let [x, y_, z] = evolution.s1(y);
            orbit.s1_x.push(x);
            orbit.s1_y.push(y_);
            orbit.s1_z.push(z);
            let [x, y_, z] = evolution.s2(y);
            orbit.s2_x.push(x);
            orbit.s2_y.push(y_);
            orbit.s2_z.push(z);
            let [x, y_, z] = Evolution::e1(y);
            orbit.e1_x.push(x);
            orbit.e1_y.push(y_);
            orbit.e1_z.push(z);
        };
        if evolution.backwards {
            samples[..cutlen].iter().rev().for_each(&mut push);
        } else {
            samples[..cutlen].iter().for_each(&mut push);
        }
        Ok(orbit)
    }

    /// [`Self::evolve_orbit`] keeping the solution at any time rather than
    /// evenly spaced samples (the integration is the same).
    pub(crate) fn evolve_dense(
        &self,
        delta_t: f64,
        f_start: f64,
        f_end: f64,
    ) -> Result<DenseOrbit, SpinTaylorError> {
        let coefficients = self.coefficients(f_start, f_end)?;
        let mut samples = DenseSamples::default();
        let evolution = self.evolve_into(&coefficients, delta_t, &mut samples)?;
        if samples.steps.is_empty() {
            return Err(SpinTaylorError::InvalidInterval);
        }
        Ok(DenseOrbit {
            coefficients,
            steps: samples.steps,
            evolution,
        })
    }

    /// [`Self::evolve_orbit`] recording the samples in `samples` rather
    /// than in a [`SpinTaylorOrbit`], for callers that need only some of the
    /// series: the orbit is given by the returned [`Evolution`].
    pub(crate) fn evolve_orbit_into<K: OrbitSamples>(
        &self,
        delta_t: f64,
        f_start: f64,
        f_end: f64,
        samples: &mut K,
    ) -> Result<Evolution, SpinTaylorError> {
        let coefficients = self.coefficients(f_start, f_end)?;
        self.evolve_into(&coefficients, delta_t, samples)
    }

    fn evolve_into<K: OrbitSamples>(
        &self,
        coefficients: &SpinTaylorT4Coefficients,
        delta_t: f64,
        samples: &mut K,
    ) -> Result<Evolution, SpinTaylorError> {
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
        let return_code = adaptive_rk4_hermite_into(
            &mut system,
            &settings,
            &mut yinit,
            0.0,
            lengths / msec,
            sgn * delta_t / msec,
            samples,
        )
        .map_err(|_| SpinTaylorError::InvalidInterval)?;
        let len = samples.len();

        // The integrator takes one adaptive stride past fEnd; keep the first
        // interpolated sample to cross it.
        let mut cutlen = len;
        if f_end != 0. && f_end < f_start {
            let w_end = PI * msec * f_end;
            while cutlen > 0 && samples.omega(cutlen - 1) < w_end {
                cutlen -= 1;
            }
            if cutlen < len {
                cutlen += 1;
            }
        } else if f_end > f_start {
            let w_end = PI * msec * f_end;
            while cutlen > 0 && samples.omega(cutlen - 1) > w_end {
                cutlen -= 1;
            }
            if cutlen < len {
                cutlen += 1;
            }
        }
        Ok(Evolution {
            len: cutlen,
            backwards,
            norm1,
            norm2,
            return_code,
        })
    }
}
