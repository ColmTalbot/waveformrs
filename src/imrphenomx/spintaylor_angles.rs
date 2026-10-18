//! Numerical "SpinTaylor" Euler angles of IMRPhenomXP/XPHM
//! (`PhenomXPrecVersion` 320, arXiv:2412.16721), following
//! `LALSimIMRPhenomX_precession.c`:
//!
//! - [`inspiral_angles_spintaylor`] ports `IMRPhenomX_InspiralAngles_SpinTaylor`,
//!   the integration of the orbit-averaged SpinTaylorT4 equations (coarse
//!   backward/forward runs around `f_ref`, then a fine run at high
//!   frequency). This is the only call into [`crate::spintaylor`].
//! - [`SpinTaylorAngles::new`] ports `IMRPhenomX_InterpolateAlphaBeta_SpinTaylor`:
//!   splines of `alpha` and `cos(beta)` against `Mf` in the J frame, the
//!   precessing final spin from the spins at the end of the integration
//!   (`PhenomXPFinalSpinMod` 4) and the coefficients of the merger-ringdown
//!   continuation (`alphaMRD_coeff`, `betaMRD_coeff`).
//! - [`SpinTaylorAngles::interpolate_gamma`] ports
//!   `IMRPhenomX_InterpolateGamma_SpinTaylor`: `gamma` from the minimal
//!   rotation condition integrated with Boole's rule
//!   (`gamma_from_alpha_cosbeta`) on a uniform grid.
//! - [`SpinTaylorAngles::angles_imr`] ports `IMRPhenomXPSpinTaylorAnglesIMR`,
//!   the angles on the frequency grid of the waveform.
//!
//! Like the SpinTaylor integrator, the arithmetic mirrors the fused
//! multiply-adds of the arm64 builds of LALSimulation and GSL (see
//! [`crate::spintaylor::fp`]): on that platform all intermediate quantities
//! (PN solution, splines, final spin, merger-ringdown coefficients, `gamma`)
//! agree with LALSimulation bit for bit.

use std::f64::consts::PI;

use crate::constants::{C_SI, G_SI, REAL4_EPS};
use crate::fastmath;
use crate::imrphenomx::internals::{IMRPhenomXError, WaveformStruct};
use crate::imrphenomx::precession::{
    delta_mf_from_wfstruct, hz_to_mf, mf_to_hz, precessing_final_spin_2017, rotate_y_cs,
    rotate_z_cs, PrecessionStruct,
};
use crate::interpolate::CubicSpline;
use crate::spintaylor::fp::{fma, fms, fnms};
use crate::spintaylor::integrator::HermiteSink;
use crate::spintaylor::{Evolution, FullSamples, OrbitSamples, SpinTaylorT4, NUM_VARIABLES};

/// Knots of the `gamma` grid integrated beyond the last frequency at which
/// the spline is evaluated, see [`SpinTaylorAngles::interpolate_gamma`].
const GAMMA_SPLINE_MARGIN: usize = 64;

/// Steps of `gamma` integrated together by the vectorised merger-ringdown
/// integration of the grid paths.
const MRD_BLOCK: usize = 256;

/// Boole's rule for `-int cos(beta) dalpha` over four intervals of `step`
/// with the integrand `a` at the five points (`7 a0 + 32 a1 + 12 a2 + 32 a3
/// + 7 a4` as contracted by clang).
#[inline]
fn boole(a: &[f64; 5], step: f64) -> f64 {
    let sum = fma(
        7.,
        a[4],
        fma(32., a[3], fma(12., a[2], fma(7., a[0], 32. * a[1]))),
    );
    -2. * step / 45. * sum
}

/// The solution of the PN precession equations used for the angles
/// (`PhenomXPInspiralArrays`): `v = (pi M f)^(1/3)` and the Newtonian
/// orbital angular momentum direction in the source frame at `f_ref` at
/// every sample, and the dimensionless spins at the last sample (the only
/// ones the angles use, for the final spin).
#[derive(Clone, Debug, Default)]
pub(crate) struct PnArrays {
    pub(crate) v: Vec<f64>,
    pub(crate) lnhat: [Vec<f64>; 3],
    pub(crate) s1_end: [f64; 3],
    pub(crate) s2_end: [f64; 3],
}

impl PnArrays {
    pub(crate) fn len(&self) -> usize {
        self.v.len()
    }
}

/// Full states kept at each end of a run recorded in [`CompactSamples`].
const KEPT_STATES: usize = 32;

/// The samples of a PN integration with what the angles need at every
/// sample, `M omega` and `LNhat`, and the full states of the first and the
/// last [`KEPT_STATES`] samples. A long integration (millions of samples
/// for binaries of a few solar masses) then takes 32 bytes per sample
/// rather than 112, and the series of [`crate::spintaylor::SpinTaylorOrbit`]
/// that are not needed are never formed.
#[derive(Debug, Default)]
struct CompactSamples {
    omega: Vec<f64>,
    lnhat: [Vec<f64>; 3],
    head: Vec<[f64; NUM_VARIABLES]>,
    /// The last [`KEPT_STATES`] states, sample `i` at `i % KEPT_STATES`.
    tail: Vec<[f64; NUM_VARIABLES]>,
}

impl HermiteSink<NUM_VARIABLES> for CompactSamples {
    fn reserve(&mut self, n: usize) {
        self.omega.reserve(n);
        self.lnhat.iter_mut().for_each(|series| series.reserve(n));
    }

    #[inline]
    fn push(&mut self, _t: f64, y: &[f64; NUM_VARIABLES]) {
        let i = self.omega.len();
        self.omega.push(y[1]);
        let [x, y_, z] = Evolution::lnhat(y);
        self.lnhat[0].push(x);
        self.lnhat[1].push(y_);
        self.lnhat[2].push(z);
        if i < KEPT_STATES {
            self.head.push(*y);
            self.tail.push(*y);
        } else {
            self.tail[i % KEPT_STATES] = *y;
        }
    }
}

impl OrbitSamples for CompactSamples {
    fn len(&self) -> usize {
        self.omega.len()
    }

    fn omega(&self, i: usize) -> f64 {
        self.omega[i]
    }
}

impl CompactSamples {
    /// The full state of recorded sample `i`, if it was kept.
    fn state(&self, i: usize) -> Option<[f64; NUM_VARIABLES]> {
        let n = self.omega.len();
        if i < self.head.len() {
            Some(self.head[i])
        } else if i < n && i + KEPT_STATES >= n {
            Some(self.tail[i % KEPT_STATES])
        } else {
            None
        }
    }
}

/// The series of `SpinTaylorOrbit` that IMRPhenomX uses, at one sample.
#[derive(Clone, Copy, Debug, Default)]
struct Point {
    v: f64,
    s1: [f64; 3],
    s2: [f64; 3],
    lnhat: [f64; 3],
    e1: [f64; 3],
}

/// One run of `XLALSimInspiralSpinTaylorPNEvolveOrbit`.
struct Run {
    binary: SpinTaylorT4,
    delta_t: f64,
    f_start: f64,
    f_end: f64,
    samples: CompactSamples,
    evolution: Evolution,
}

impl Run {
    fn evolve(
        binary: &SpinTaylorT4,
        delta_t: f64,
        f_start: f64,
        f_end: f64,
    ) -> Result<Self, IMRPhenomXError> {
        let mut samples = CompactSamples::default();
        let evolution = binary
            .evolve_orbit_into(delta_t, f_start, f_end, &mut samples)
            .map_err(|_| {
                IMRPhenomXError::SpinTaylorFallbackToMsa(
                    "XLALSimInspiralSpinTaylorPNEvolveOrbit failed",
                )
            })?;
        Ok(Self {
            binary: *binary,
            delta_t,
            f_start,
            f_end,
            samples,
            evolution,
        })
    }

    /// Number of samples of the orbit.
    fn len(&self) -> usize {
        self.evolution.len
    }

    /// Sample `j` of the orbit (ordered forward in time).
    fn point(&self, j: usize) -> Point {
        let i = self.evolution.recorded(j);
        let y = self.samples.state(i).unwrap_or_else(|| {
            // Only states near the ends of the runs are read (see
            // `inspiral_angles_spintaylor`), so this is not expected:
            // integrate again keeping every state.
            let mut full = FullSamples::default();
            self.binary
                .evolve_orbit_into(self.delta_t, self.f_start, self.f_end, &mut full)
                .expect("the integration succeeded before");
            full.0[i]
        });
        Point {
            v: Evolution::v(y[1]),
            s1: self.evolution.s1(&y),
            s2: self.evolution.s2(&y),
            lnhat: Evolution::lnhat(&y),
            e1: Evolution::e1(&y),
        }
    }

    /// Append `v` and `LNhat` of the first `n` samples of the orbit.
    fn extend_series(&self, v: &mut Vec<f64>, lnhat: &mut [Vec<f64>; 3], n: usize) {
        for j in 0..n {
            let i = self.evolution.recorded(j);
            v.push(Evolution::v(self.samples.omega[i]));
            for (series, recorded) in lnhat.iter_mut().zip(&self.samples.lnhat) {
                series.push(recorded[i]);
            }
        }
    }

    /// `v` and `LNhat` of the first `n` samples of the orbit, reusing the
    /// recorded series, with room for `capacity` samples.
    fn into_series(self, n: usize, capacity: usize) -> (Vec<f64>, [Vec<f64>; 3]) {
        let Self {
            samples, evolution, ..
        } = self;
        let CompactSamples {
            mut omega,
            mut lnhat,
            ..
        } = samples;
        for series in std::iter::once(&mut omega).chain(lnhat.iter_mut()) {
            series.truncate(evolution.len);
            if evolution.backwards {
                series.reverse();
            }
            series.truncate(n);
            series.reserve_exact(capacity.saturating_sub(n));
        }
        omega.iter_mut().for_each(|w| *w = Evolution::v(*w));
        (omega, lnhat)
    }
}

/// The first, coarse part of the PN integration: a forward run from `f_ref`,
/// after a backward run from `f_ref` unless `f_ref = f_min`. The two runs
/// are joined by `appendTS` of `LALSimIMRPhenomX_precession.c`: the forward
/// run is copied over the last two samples of the backward one and the
/// result has length `len1 + len2 - 1`, ending with a zero sample.
struct Coarse {
    backward: Option<Run>,
    forward: Run,
}

impl Coarse {
    /// Number of samples of the backward run before the forward one.
    fn split(&self) -> usize {
        self.backward.as_ref().map_or(0, |run| run.len() - 2)
    }

    fn len(&self) -> usize {
        match &self.backward {
            None => self.forward.len(),
            Some(_) => self.split() + self.forward.len() + 1,
        }
    }

    fn point(&self, i: usize) -> Point {
        let split = self.split();
        match &self.backward {
            Some(backward) if i < split => backward.point(i),
            _ if i < split + self.forward.len() => self.forward.point(i - split),
            _ => Point::default(),
        }
    }

    /// `v` and `LNhat` of the first `n < len()` samples, with room for
    /// `capacity` samples.
    fn into_series(self, n: usize, capacity: usize) -> (Vec<f64>, [Vec<f64>; 3]) {
        let split = self.split().min(n);
        match self.backward {
            None => self.forward.into_series(n, capacity),
            Some(backward) => {
                let (mut v, mut lnhat) = backward.into_series(split, capacity);
                self.forward.extend_series(&mut v, &mut lnhat, n - split);
                (v, lnhat)
            }
        }
    }
}

/// `IMRPhenomX_InspiralAngles_SpinTaylor` for black-hole binaries: evolve the
/// SpinTaylorT4 equations from `f_min` (Hz, already including the buffer
/// below the first frequency of the waveform) through `f_ref` up to
/// `(fRING + 8 fDAMP) / M` of the aligned-spin model.
///
/// Every failure is reported as
/// [`IMRPhenomXError::SpinTaylorFallbackToMsa`]: LALSimulation then falls
/// back to the MSA angles.
///
/// Only `v` and `LNhat` are kept at every sample. The full states read are
/// near the ends of the runs: the fine run starts 10 samples before the end
/// of the coarse solution (within the first 11 samples of the backward run
/// if the forward run is shorter than 10 samples), and the angles use the
/// spins at the end of the fine run.
pub(crate) fn inspiral_angles_spintaylor(
    wf: &WaveformStruct,
    prec: &PrecessionStruct,
    f_min: f64,
) -> Result<PnArrays, IMRPhenomXError> {
    let fallback = IMRPhenomXError::SpinTaylorFallbackToMsa;
    let f_ref = wf.f_ref;
    if f_min <= 0.0 {
        return Err(fallback("fmin must be positive"));
    }
    if f_ref < f_min {
        return Err(fallback("fRef must be >= fmin"));
    }
    let (m1_si, m2_si) = (wf.m1_si, wf.m2_si);
    let pi_gm = PI * (m1_si + m2_si) * (G_SI / C_SI) / (C_SI * C_SI);
    let binary = SpinTaylorT4::new(m1_si, m2_si, prec.chi1, prec.chi2);

    // For versions 32* to work, the integration must start in the inspiral.
    let f_meco_hz = mf_to_hz(wf.f_meco, wf.m_tot);
    let mut f_min = f_min;
    if f_min > f_meco_hz && (prec.version == 320 || prec.version == 321) {
        f_min = f_meco_hz;
    }
    let f_cut = mf_to_hz(prec.f_ring_aligned + 8. * prec.f_damp_aligned, wf.m_tot);
    // PhenomXPSpinTaylorCoarseFactor, read as a float.
    let coarse_fac = 10.0f32 as f64;
    let delta_t_coarse = 0.5 * coarse_fac / f_cut;

    let coarse = if f_ref < REAL4_EPS || (f_ref - f_min).abs() < REAL4_EPS {
        // f_ref is (re)set to f_min: a single forward integration.
        Coarse {
            backward: None,
            forward: Run::evolve(&binary, delta_t_coarse, f_min, f_cut)?,
        }
    } else if (f_ref - f_min) > REAL4_EPS {
        // Integrate backward to fStart, then forward to the end of the waveform.
        let backward = Run::evolve(&binary, delta_t_coarse, f_ref, f_min - 0.5)?;
        if backward.len() <= 1 {
            return Err(fallback("backward PN integration returned a single sample"));
        }
        let forward = Run::evolve(&binary, delta_t_coarse, f_ref, f_cut)?;
        Coarse {
            backward: Some(backward),
            forward,
        }
    } else {
        return Err(fallback("fRef must be >= fmin"));
    };

    // At high frequencies, perform a fine integration, starting some bins
    // before the last one.
    let len_low = coarse.len();
    let nbuffer = 9.min(len_low - 1);
    let copy_length = len_low - 1 - nbuffer;
    let at = coarse.point(copy_length);
    let f_trans = at.v.powf(3.) / pi_gm;
    if f_trans.is_nan() || f_trans <= 0. {
        return Err(fallback(
            "transition frequency in PN integration is not positive",
        ));
    }
    let fine_binary =
        SpinTaylorT4::new(m1_si, m2_si, at.s1, at.s2).with_orientation(at.lnhat, at.e1);
    let fine = Run::evolve(&fine_binary, 0.5 / f_cut, f_trans, f_cut)?;
    let len_pn = copy_length + fine.len();
    if len_pn < 4 {
        return Err(fallback("too few points for spline interpolation"));
    }
    let (mut v, mut lnhat) = coarse.into_series(copy_length, len_pn);
    fine.extend_series(&mut v, &mut lnhat, fine.len());
    let end = fine.point(fine.len() - 1);
    let arrays = PnArrays {
        v,
        lnhat,
        s1_end: end.s1,
        s2_end: end.s2,
    };

    // The first frequency node must lie below the requested fmin to avoid
    // interpolation errors.
    let f_min_pn = arrays.v[0].powf(3.) / pi_gm;
    if f_min_pn < 0. || f_min_pn > f_min {
        return Err(fallback("PN solution starts above fmin"));
    }
    Ok(arrays)
}

/// `evaluate_QNMfit_fdamp21` of `LALSimIMRPhenomXHM_qnm.c`, the damping
/// frequency of the (2,1) mode (accumulated with fused multiply-adds like
/// the arm64 builds of LALSimulation, see `qnm.rs`).
fn fdamp21(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let num = [
        (-3.1358071947583093, af),
        (0.14456081596393977, x2),
        (1.2194717985037946, x3),
        (-0.2947372598589144, x4),
        (0.002943057145913646, x5),
    ]
    .iter()
    .fold(2.0696914454467294, |acc: f64, &(c, x): &(f64, f64)| {
        c.mul_add(x, acc)
    });
    let den = [
        (-219.81790388304876, af),
        (17.7141194900164, x2),
        (75.90115083917898, x3),
        (-18.975287709794745, x4),
    ]
    .iter()
    .fold(146.1779212636481, |acc: f64, &(c, x): &(f64, f64)| {
        c.mul_add(x, acc)
    });
    num / den
}

/// Coefficients of the merger-ringdown continuation of `alpha`,
/// `PhenomXPalphaMRD`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct AlphaMrd {
    a: f64,
    b: f64,
    c: f64,
}

impl AlphaMrd {
    /// `alphaMRD_coeff`: match `-(a + b / Mf^4 + c / Mf^2)` to the spline
    /// (value and derivative at `0.97 fmaxPN`, value at `0.99 fmaxPN`).
    fn new(spline: &CubicSpline, fmax_pn: f64) -> Result<Self, IMRPhenomXError> {
        let ftrans = fmax_pn;
        let f1 = 0.97 * ftrans;
        let f2 = 0.99 * ftrans;
        let f1sq = f1 * f1;
        let f2sq = f2 * f2;
        let f1cube = f1sq * f1;
        if !(spline.contains(f1) && spline.contains(f2)) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "alpha could not be interpolated for the MRD continuation",
            ));
        }
        let alpha1 = -spline.eval(f1);
        let dalpha1 = -spline.derivative(f1);
        let alpha2 = -spline.eval(f2);
        // The expressions of LAL with clang's contraction into fused
        // multiply-adds (see `crate::spintaylor::fp`).
        let f1_4 = f1.powf(4.);
        let f2_4 = f2.powf(4.);
        let den = 2. * (f1sq - f2sq).powf(2.);
        let a = fma(
            2. * f2_4,
            alpha2,
            fma(
                f1cube * (f1 - f2) * (f1 + f2),
                dalpha1,
                2. * fnms(2. * f1sq, f2sq, f1_4) * alpha1,
            ),
        ) / den;
        let b = (f1_4
            * f2sq
            * fma(
                f1 * (f1 - f2) * (f1 + f2),
                dalpha1,
                2. * f2sq * (-alpha1 + alpha2),
            ))
            / den;
        let c = (f1sq * fma(f1 * (-f1_4 + f2_4), dalpha1, 4. * f2_4 * (alpha1 - alpha2))) / den;
        Ok(Self { a, b, c })
    }

    /// `alphaMRD`.
    #[inline]
    pub(crate) fn alpha(&self, mf: f64) -> f64 {
        let mf2 = mf * mf;
        let mf4 = mf2 * mf2;
        -(self.a + self.b / mf4 + self.c / mf2)
    }

    /// `dalphaMRD`.
    #[inline]
    pub(crate) fn dalpha(&self, mf: f64) -> f64 {
        let mf2 = mf * mf;
        let mf3 = mf2 * mf;
        let mf5 = mf2 * mf3;
        4. * self.b / mf5 + 2. * self.c / mf3
    }

    /// Multiply each of `out` by [`Self::dalpha`] at the corresponding `mf`,
    /// with one division (for the grid paths; agrees to rounding).
    fn mul_dalpha_fast(&self, mf: &[f64], out: &mut [f64]) {
        for (out, &mf) in out.iter_mut().zip(mf) {
            let inv = 1.0 / mf;
            let inv2 = inv * inv;
            *out *= inv2 * inv * fma(4. * self.b, inv2, 2. * self.c);
        }
    }
}

/// Coefficients of the merger-ringdown continuation of `beta`,
/// `PhenomXPbetaMRD`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct BetaMrd {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    cosbeta_sign: f64,
    flat: bool,
    dfdamp: f64,
}

impl BetaMrd {
    /// `betaMRD_coeff`: match
    /// `exp(-kappa Mf) / Mf (a / Mf + b / Mf^2 + c / Mf^3) + d` to the spline
    /// of `cos(beta)` at `0.97 fmaxPN` and `0.98 fmaxPN`, with
    /// `kappa = 2 pi (fdamp21 - fdamp22)` and `d = 0` or `pi`.
    fn new(spline: &CubicSpline, fmax_pn: f64, fmax_inspiral: f64, wf: &WaveformStruct) -> Self {
        let fdamp = wf.f_damp;
        let dfdamp = fdamp21(wf.a_final) / wf.m_final - fdamp;
        let kappa = 2. * PI * dfdamp;
        let f1 = 0.97 * fmax_pn;
        let f2 = 0.98 * fmax_pn;
        let f1sq = f1 * f1;
        let f2sq = f2 * f2;
        let ef1 = (kappa * f1).exp();
        let ef2 = (kappa * f2).exp();

        // As in LAL, only the status of the last evaluation is checked.
        let eval = |x: f64| {
            if spline.contains(x) {
                spline.eval(x)
            } else {
                f64::NAN
            }
        };
        let cosbeta1 = eval(f1);
        let dcosbeta2 = if spline.contains(f2) {
            spline.derivative(f2)
        } else {
            f64::NAN
        };
        let cosbeta2 = eval(f2);
        let success = spline.contains(fmax_inspiral);
        let cosbetamax = eval(fmax_inspiral);

        let mut beta = Self {
            dfdamp,
            cosbeta_sign: 1.0f64.copysign(cosbetamax),
            ..Default::default()
        };
        if cosbeta1.abs() > 1. || cosbeta2.abs() > 1. || cosbetamax.abs() > 1. || !success {
            beta.flat = true;
        } else {
            let beta1 = cosbeta1.acos();
            let beta2 = cosbeta2.acos();
            let sqrtarg = fnms(cosbeta2, cosbeta2, 1.);
            let dbeta2 = -dcosbeta2 / (if sqrtarg <= 0. { 1. } else { sqrtarg }).sqrt();
            let off = if cosbetamax < 0. { PI } else { 0. };
            let den = (f1 - f2).powf(2.);
            // The expressions of LAL with clang's contraction into fused
            // multiply-adds (see `crate::spintaylor::fp`).
            let f1_4 = f1.powf(4.);
            let f2_3 = f2.powf(3.);
            let f2_4 = f2.powf(4.);
            let k2 = fma(f2, kappa, 2.);
            let k3 = fma(f2, kappa, 3.);
            let k4 = fma(f2, kappa, 4.);
            let w_a = fma(-f2, k3, f1 * k4);
            let z_a = fma(f2 * (-f1 + f2), dbeta2, w_a * (off - beta2));
            beta.a = fma(-(ef1 * f1_4), off - beta1, ef2 * f2_3 * z_a) / den;
            let w_b = fma(-f2sq, k2, f1sq * k4);
            let z_b = fms((f1 - f2) * f2 * (f1 + f2), dbeta2, w_b * (off - beta2));
            beta.b = fma(2. * ef1 * f1_4 * f2, off - beta1, ef2 * f2_3 * z_b) / den;
            let w_c = fma(-f2, k2, f1 * k3);
            let z_c = fma(f2 * (-f1 + f2), dbeta2, w_c * (off - beta2));
            beta.c = fma(-(ef1 * f1_4 * f2sq), off - beta1, ef2 * f1 * f2_4 * z_c) / den;
            beta.d = off;
            beta.flat = false;
        }
        beta
    }

    /// `betaMRD`.
    #[inline]
    pub(crate) fn beta(&self, mf: f64) -> f64 {
        if self.flat {
            self.cosbeta_sign.acos()
        } else {
            let kappa = 2. * PI * self.dfdamp;
            let mf2 = mf * mf;
            let mf3 = mf2 * mf;
            fma(
                (-mf * kappa).exp() / mf,
                self.a / mf + self.b / mf2 + self.c / mf3,
                self.d,
            )
        }
    }

    /// `cos(betaMRD)` at each of `mf`, written to `out`, with the inlined
    /// `exp` and `cos` in vectorisable loops and one division (for the grid
    /// paths; agrees with `self.beta(mf).cos()` to rounding).
    ///
    /// With `uniform` the points are equally spaced up to rounding, and
    /// `exp(-kappa Mf)` is the exponential at the first of every
    /// `EXP_RESTART` points times powers of `exp(-kappa dMf)` (the
    /// rounding of the powers stays below `EXP_RESTART` ulp).
    fn cos_beta_fast(&self, mf: &[f64], out: &mut [f64], uniform: bool) {
        /// Points between exact exponentials with `uniform`.
        const EXP_RESTART: usize = 64;
        let out = &mut out[..mf.len()];
        if self.flat {
            out.fill(self.cosbeta_sign.acos().cos());
            return;
        }
        let kappa = 2. * PI * self.dfdamp;
        if uniform && mf.len() > EXP_RESTART {
            let spacing = (mf[mf.len() - 1] - mf[0]) / (mf.len() - 1) as f64;
            let ratio = (-kappa * spacing).exp();
            let mut powers = [1.0; EXP_RESTART];
            for k in 1..EXP_RESTART {
                powers[k] = powers[k - 1] * ratio;
            }
            for (out, mf) in out.chunks_mut(EXP_RESTART).zip(mf.chunks(EXP_RESTART)) {
                let start = (-kappa * mf[0]).exp();
                for (out, power) in out.iter_mut().zip(&powers) {
                    *out = start * power;
                }
            }
        } else {
            for (out, &mf) in out.iter_mut().zip(mf) {
                *out = -mf * kappa;
            }
            fastmath::exp_in_place(out);
        }
        for (out, &mf) in out.iter_mut().zip(mf) {
            // exp(-kappa Mf) / Mf (a / Mf + b / Mf^2 + c / Mf^3) + d
            let inv = 1.0 / mf;
            let sum = inv * fma(inv, fma(self.c, inv, self.b), self.a);
            *out = fma(*out * inv, sum, self.d);
        }
        fastmath::cos_in_place(out);
    }
}

/// The SpinTaylor Euler angles of one binary (the angle-related members of
/// `IMRPhenomXPrecessionStruct`): splines of `alpha` and `cos(beta)` over
/// the PN solution, the merger-ringdown continuation and, once
/// [`SpinTaylorAngles::interpolate_gamma`] is called, the spline of `gamma`.
#[derive(Clone, Debug)]
pub(crate) struct SpinTaylorAngles {
    /// `alpha` (unwrapped) and `cos(beta)` against `Mf`.
    pub(crate) alpha_spline: CubicSpline,
    pub(crate) cosbeta_spline: CubicSpline,
    /// `gamma` against `Mf` on a uniform grid.
    pub(crate) gamma_spline: Option<CubicSpline>,
    /// End of the inspiral angles and start of the MRD continuation.
    #[allow(dead_code)] // kept as in LAL
    pub(crate) fmax_inspiral: f64,
    pub(crate) ftrans_mrd: f64,
    pub(crate) alpha_mrd: AlphaMrd,
    pub(crate) beta_mrd: BetaMrd,
    /// `alpha` and `gamma` at the reference frequency.
    pub(crate) alpha_ref: f64,
    pub(crate) gamma_ref: f64,
    /// Precession version (320 or 321 continue the angles analytically).
    version: u32,
}

impl SpinTaylorAngles {
    /// `IMRPhenomX_InterpolateAlphaBeta_SpinTaylor` with
    /// `PhenomXPFinalSpinMod` 4: also sets the final spin of `wf` from the
    /// spins at the end of the PN integration and updates `fRING`, `fDAMP`.
    ///
    /// The PN arrays are overwritten with `Mf`, `alpha` and `cos(beta)`,
    /// which become the knots and values of the splines.
    pub(crate) fn new(
        wf: &mut WaveformStruct,
        prec: &PrecessionStruct,
        arrays: PnArrays,
    ) -> Result<Self, IMRPhenomXError> {
        let len_pn = arrays.len();
        let i_max = len_pn - 1;
        let PnArrays {
            v: mut fgw,
            lnhat: [mut alpha, lnhat_y, mut cosbeta],
            s1_end,
            s2_end,
        } = arrays;
        // the direction of L at the end of the integration, in the source frame
        let lnf = [alpha[i_max], lnhat_y[i_max], cosbeta[i_max]];
        let (cz1, sz1) = ((-prec.phi_j_sf).cos(), (-prec.phi_j_sf).sin());
        let (cy, sy) = ((-prec.theta_j_sf).cos(), (-prec.theta_j_sf).sin());
        let (cz2, sz2) = ((-prec.kappa).cos(), (-prec.kappa).sin());
        let mut mfmax_pn = 0.;
        for i in 0..len_pn {
            let lnhat = [alpha[i], lnhat_y[i], cosbeta[i]];
            let lnhat = rotate_z_cs(cz2, sz2, rotate_y_cs(cy, sy, rotate_z_cs(cz1, sz1, lnhat)));
            let fgw_hz = fgw[i].powf(3.) / prec.pi_gm;
            let fgw_mf = hz_to_mf(fgw_hz, wf.m_tot);
            if fgw_hz.is_nan() || fgw_hz <= 0. {
                // LAL leaves the remaining nodes uninitialised.
                return Err(IMRPhenomXError::SpinTaylorAngles(
                    "non-positive frequency in the PN solution",
                ));
            }
            // Euler angles in the J frame (alpha before unwrapping)
            alpha[i] = lnhat[1].atan2(lnhat[0]);
            cosbeta[i] = lnhat[2];
            fgw[i] = fgw_mf;
            mfmax_pn = fgw_mf;
        }
        drop(lnhat_y);
        if fgw
            .windows(2)
            .any(|w| w[0] >= w[1] || w[0].is_nan() || w[1].is_nan())
        {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "PN frequencies are not strictly increasing",
            ));
        }

        let mut fmax_inspiral = mfmax_pn - wf.delta_mf;
        if fmax_inspiral > wf.f_ring - wf.f_damp {
            fmax_inspiral = 1.020 * wf.f_meco;
        }
        let ftrans_mrd = 0.98 * fmax_inspiral;

        unwrap_in_place(&mut alpha);
        let (alpha_spline, cosbeta_spline) = CubicSpline::pair(fgw, alpha, cosbeta);
        if !cosbeta_spline.contains(fmax_inspiral) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "cos(beta) could not be evaluated at the end of the inspiral",
            ));
        }
        let cosbetamax = cosbeta_spline.eval(fmax_inspiral);

        // Estimate the final spin from the spins at the end of the PN
        // integration (PhenomXPFinalSpinMod 4).
        let m1 = wf.m1_si / wf.m_tot_si;
        let m2 = wf.m2_si / wf.m_tot_si;
        let dot = |a: [f64; 3], b: [f64; 3]| crate::spintaylor::fp::dot(&a, &b);
        let lnorm = dot(lnf, lnf).sqrt();
        let (s1f, s2f) = (s1_end, s2_end);
        let dot_s1l = dot(s1f, lnf) / lnorm;
        let dot_s2l = dot(s2f, lnf) / lnorm;
        let perp = |s: [f64; 3], sl: f64, mass: f64| {
            let scale = mass * mass;
            std::array::from_fn::<f64, 3, _>(|k| scale * (s[k] - sl * lnf[k]))
        };
        let s1_perp = perp(s1f, dot_s1l, m1);
        let s2_perp = perp(s2f, dot_s2l, m2);
        let stot_perp: [f64; 3] = std::array::from_fn(|k| s1_perp[k] + s2_perp[k]);
        let s_perp_norm = dot(stot_perp, stot_perp).sqrt();
        let chi_perp_norm = s_perp_norm * (m1 + m2).powf(2.) / m1.powf(2.);
        let a_final = 1.0f64.copysign(cosbetamax)
            * precessing_final_spin_2017(wf.eta, dot_s1l, dot_s2l, chi_perp_norm);
        if a_final.is_nan() || a_final.abs() > 1.0 {
            // LAL evaluates the QNM fits outside of their domain.
            return Err(IMRPhenomXError::UnphysicalFinalSpin(a_final));
        }
        wf.set_final_spin(a_final);

        let (alpha_mrd, beta_mrd) = if matches!(prec.version, 320 | 321 | 330) {
            (
                AlphaMrd::new(&alpha_spline, fmax_inspiral)?,
                BetaMrd::new(&cosbeta_spline, fmax_inspiral, fmax_inspiral, wf),
            )
        } else {
            (AlphaMrd::default(), BetaMrd::default())
        };

        Ok(Self {
            alpha_spline,
            cosbeta_spline,
            gamma_spline: None,
            fmax_inspiral,
            ftrans_mrd,
            alpha_mrd,
            beta_mrd,
            alpha_ref: 0.0,
            gamma_ref: 0.0,
            version: prec.version,
        })
    }

    fn continues_analytically(&self) -> bool {
        self.version == 320 || self.version == 321
    }

    /// `gamma_from_alpha_cosbeta`: the change of `gamma` from `Mf` to
    /// `Mf + delta_mf` from the minimal rotation condition
    /// `dgamma/df = -cos(beta) dalpha/df`, with Boole's rule. `hint` caches
    /// the spline interval between calls (`gsl_interp_accel`).
    #[inline]
    pub(crate) fn gamma_step(
        &self,
        mf: f64,
        delta_mf: f64,
        hint: &mut usize,
    ) -> Result<f64, IMRPhenomXError> {
        self.gamma_step_cached(mf, delta_mf, hint, &mut BooleCache::default())
    }

    /// [`Self::gamma_step`] reusing the integrand at the last point of the
    /// previous step when it is the first point of this one (on a uniform
    /// grid `Mf_i + 4 (dMf / 4)` often equals `Mf_{i+1}` exactly), which
    /// gives the same result with fewer evaluations.
    #[inline]
    fn gamma_step_cached(
        &self,
        mf: f64,
        delta_mf: f64,
        hint: &mut usize,
        cache: &mut BooleCache,
    ) -> Result<f64, IMRPhenomXError> {
        let mf_high = mf + delta_mf;
        let step = (mf_high - mf) * 0.25;
        let mut alphadotcosbeta = [0.0; 5];
        let inspiral = mf <= self.ftrans_mrd;
        // fma(0, step, mf) = mf
        let first = cache.lookup(mf, inspiral);
        if inspiral {
            for (j, value) in alphadotcosbeta.iter_mut().enumerate() {
                if j == 0 {
                    if let Some(cached) = first {
                        *value = cached;
                        continue;
                    }
                }
                let mf_aux = fma(j as f64, step, mf);
                if !self.alpha_spline.contains(mf_aux) {
                    return Err(IMRPhenomXError::SpinTaylorAngles(
                        "could not evaluate splines for alpha and/or gamma angles",
                    ));
                }
                // alpha and cos(beta) share their knots.
                *hint = self.alpha_spline.locate(mf_aux, *hint);
                let cosbeta = self.cosbeta_spline.eval_in(*hint, mf_aux);
                let alphadot = self.alpha_spline.derivative_in(*hint, mf_aux);
                *value = cosbeta * alphadot;
            }
        } else {
            for (j, value) in alphadotcosbeta.iter_mut().enumerate() {
                if j == 0 {
                    if let Some(cached) = first {
                        *value = cached;
                        continue;
                    }
                }
                let mf_aux = fma(j as f64, step, mf);
                let cosbeta = self.beta_mrd.beta(mf_aux).cos();
                let alphadot = self.alpha_mrd.dalpha(mf_aux);
                *value = cosbeta * alphadot;
            }
        }
        *cache = BooleCache {
            mf: fma(4., step, mf),
            inspiral,
            value: alphadotcosbeta[4],
            valid: true,
        };
        Ok(boole(&alphadotcosbeta, step))
    }

    /// For the grid paths: `cos(beta)` of the merger-ringdown continuation
    /// at each of `mf`, written to `out` (see [`BetaMrd::cos_beta_fast`]).
    pub(crate) fn mrd_cos_beta_fast(&self, mf: &[f64], out: &mut [f64]) {
        self.beta_mrd.cos_beta_fast(mf, out, false);
    }

    /// For the grid paths: the merger-ringdown integrand `cos(beta)
    /// dalpha/df` of [`Self::gamma_step`] at each of `mf` (equally spaced
    /// with `uniform`), written to `out`.
    fn mrd_integrand_fast(&self, mf: &[f64], out: &mut [f64], uniform: bool) {
        self.beta_mrd.cos_beta_fast(mf, out, uniform);
        self.alpha_mrd.mul_dalpha_fast(mf, out);
    }

    /// For the grid paths: the steps of `gamma` in the merger-ringdown from
    /// each of `mf` to `mf + delta_mf` ([`Self::gamma_step`] above
    /// `ftrans_MRD`) with the vectorised integrand, added in turn to
    /// `gamma` (starting from `gamma_start`) and pushed to `out`. With
    /// `shared_ends` the steps follow each other (`mf + delta_mf` is the
    /// next `mf` up to rounding) and the integrand at the end of each step
    /// is reused for the next.
    fn push_mrd_gamma_fast(
        &self,
        mf: &[f64],
        delta_mf: impl Fn(usize) -> f64,
        shared_ends: bool,
        gamma_start: f64,
        out: &mut Vec<f64>,
    ) {
        let mut points = Vec::with_capacity(5 * MRD_BLOCK);
        let mut values = vec![0.0; 5 * MRD_BLOCK];
        let mut gamma = gamma_start;
        let mut previous_end: Option<f64> = None;
        let points_per_step = if shared_ends { 4 } else { 5 };
        for (block_index, block) in mf.chunks(MRD_BLOCK).enumerate() {
            let offset = block_index * MRD_BLOCK;
            points.clear();
            let first_point = shared_ends && previous_end.is_none();
            if first_point {
                points.push(block[0]);
            }
            for (k, &mf) in block.iter().enumerate() {
                let step = ((mf + delta_mf(offset + k)) - mf) * 0.25;
                for j in 5 - points_per_step..5 {
                    points.push(fma(j as f64, step, mf));
                }
            }
            let values = &mut values[..points.len()];
            // with shared ends the points are equally spaced
            self.mrd_integrand_fast(&points, values, shared_ends);
            let mut next = usize::from(first_point);
            if first_point {
                previous_end = Some(values[0]);
            }
            for (k, &mf) in block.iter().enumerate() {
                let step = ((mf + delta_mf(offset + k)) - mf) * 0.25;
                let a = if shared_ends {
                    let a0 = previous_end.expect("integrand at the start");
                    [
                        a0,
                        values[next],
                        values[next + 1],
                        values[next + 2],
                        values[next + 3],
                    ]
                } else {
                    [
                        values[next],
                        values[next + 1],
                        values[next + 2],
                        values[next + 3],
                        values[next + 4],
                    ]
                };
                next += points_per_step;
                previous_end = Some(a[4]);
                gamma += boole(&a, step);
                out.push(gamma);
            }
        }
    }

    /// `IMRPhenomX_InterpolateGamma_SpinTaylor`: integrate `gamma` on the
    /// uniform grid `Mf_min + i dMf` (with the spacing of
    /// `get_deltaF_from_wfstruct`) covering `[f_min, f_max]` (Hz), build its
    /// spline and evaluate `gamma` at the reference frequency.
    ///
    /// If the spline will only be evaluated up to `mf_needed`, the grid is
    /// truncated [`GAMMA_SPLINE_MARGIN`] knots beyond it. The coefficients of
    /// a natural cubic spline depend on a knot `k` positions away through a
    /// factor `~(2 - sqrt(3))^k`, `~1e-37` for 64 knots, so the spline below
    /// `mf_needed` is unchanged to the last bit while the (expensive)
    /// merger-ringdown part of the integration is skipped.
    ///
    /// With `fast` (for the grid paths) the merger-ringdown part (versions
    /// 320 and 321) is integrated with the vectorised integrand of
    /// [`Self::push_mrd_gamma_fast`], which agrees to rounding.
    pub(crate) fn interpolate_gamma(
        &mut self,
        wf: &WaveformStruct,
        f_min: f64,
        f_max: f64,
        mf_needed: Option<f64>,
        fast: bool,
    ) -> Result<(), IMRPhenomXError> {
        let mf_min = hz_to_mf(f_min, wf.m_tot);
        if f_min <= 0.0 || f_max <= 0.0 || f_max <= f_min || wf.f_ref < f_min {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "SpinTaylor angles need 0 < f_min <= f_ref and f_min < f_max",
            ));
        }
        let delta_mf = delta_mf_from_wfstruct(wf);
        let delta_f = mf_to_hz(delta_mf, wf.m_tot);
        // length of the frequency series expected by the user
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max / delta_f) as usize + 1;
        let output_length = i_stop - i_start;
        if output_length < 4 {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "too few points for the spline interpolation of gamma",
            ));
        }
        let mf_last = fma((output_length - 1) as f64, delta_mf, mf_min);
        if wf.mf_ref.is_nan() || wf.mf_ref < mf_min || wf.mf_ref > mf_last {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "gamma could not be evaluated at the reference frequency",
            ));
        }
        let n_knots = match mf_needed {
            Some(mf_needed) => {
                let last_needed = ((mf_needed - mf_min) / delta_mf).ceil().max(0.0) as usize;
                (last_needed + 1 + GAMMA_SPLINE_MARGIN).min(output_length)
            }
            None => output_length,
        };
        let mut frequencies = Vec::with_capacity(n_knots);
        let mut gamma = Vec::with_capacity(n_knots);
        frequencies.push(mf_min);
        gamma.push(0.);
        let mut hint = 0;
        let mut cache = BooleCache::default();
        let fast = fast && self.continues_analytically();
        for i in 1..n_knots {
            let mf = fma(i as f64, delta_mf, mf_min);
            if fast && mf > self.ftrans_mrd {
                // the merger-ringdown steps at the remaining knots
                let knots: Vec<f64> = (i..n_knots)
                    .map(|i| fma(i as f64, delta_mf, mf_min))
                    .collect();
                let start = gamma[i - 1];
                self.push_mrd_gamma_fast(&knots, |_| delta_mf, true, start, &mut gamma);
                frequencies.extend_from_slice(&knots);
                break;
            }
            frequencies.push(mf);
            let previous = gamma[i - 1];
            let deltagamma = if mf < self.ftrans_mrd || self.continues_analytically() {
                self.gamma_step_cached(mf, delta_mf, &mut hint, &mut cache)?
            } else {
                0.0
            };
            gamma.push(previous + deltagamma);
        }
        let spline = CubicSpline::from_vecs(frequencies, gamma);
        self.gamma_ref = spline.eval_unchecked(wf.mf_ref);
        self.gamma_spline = Some(spline);
        Ok(())
    }

    /// `alpha` at the reference frequency: from the spline in the inspiral,
    /// from the MRD continuation for versions 32*, otherwise the value at
    /// the transition frequency.
    pub(crate) fn set_alpha_ref(&mut self, mf_ref: f64) -> Result<(), IMRPhenomXError> {
        let error = IMRPhenomXError::SpinTaylorAngles("could not evaluate alpha(f_ref)");
        self.alpha_ref = if mf_ref <= self.ftrans_mrd {
            if !self.alpha_spline.contains(mf_ref) {
                return Err(error);
            }
            self.alpha_spline.eval(mf_ref)
        } else if self.continues_analytically() {
            self.alpha_mrd.alpha(mf_ref)
        } else {
            if !self.alpha_spline.contains(self.ftrans_mrd) {
                return Err(error);
            }
            self.alpha_spline.eval(self.ftrans_mrd)
        };
        Ok(())
    }

    /// `IMRPhenomXPSpinTaylorAnglesIMR`: build the spline of `gamma` between
    /// the first and last frequencies and return `(alpha, cos(beta), gamma)`
    /// on the grid `frequencies` (Hz, increasing), with `alpha` offset to
    /// `alpha0` at the reference frequency and `gamma` not yet offset (the
    /// twisting up subtracts `gamma_ref`).
    ///
    /// In the merger-ringdown `gamma` is integrated from point to point,
    /// using (like LAL) the step to the previous frequency forward from the
    /// current one.
    ///
    /// Angles are only returned up to the last frequency with `Mf <= mf_cut`
    /// (the waveform vanishes above), which leaves those values unchanged.
    ///
    /// With `fast` (for the grid paths) `cos(beta)` and the steps of `gamma`
    /// in the merger-ringdown (versions 320 and 321) use the vectorised
    /// functions of [`Self::push_mrd_gamma_fast`], which agree to rounding.
    pub(crate) fn angles_imr(
        &mut self,
        wf: &WaveformStruct,
        prec: &PrecessionStruct,
        frequencies: &[f64],
        mf_cut: f64,
        fast: bool,
    ) -> Result<EulerAngles, IMRPhenomXError> {
        let n = frequencies.len();
        let f_min = frequencies[0];
        let f_max = frequencies[n - 1];
        let mf_min = hz_to_mf(f_min, wf.m_tot);
        if f_min <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies("fmin must be positive"));
        }
        if f_max <= f_min {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "fmax must be larger than fmin",
            ));
        }
        if wf.f_ref < f_min {
            return Err(IMRPhenomXError::InvalidFrequencies("fRef must be >= fmin"));
        }
        // The gamma spline is evaluated below ftrans_MRD and at f_ref.
        let mf_needed = Some(self.ftrans_mrd.max(wf.mf_ref));
        self.interpolate_gamma(wf, f_min, f_max, mf_needed, fast)?;
        if !self.alpha_spline.contains(mf_min) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "could not evaluate alpha(f_min)",
            ));
        }
        let alphamin = self.alpha_spline.eval(mf_min);
        let cosbetamin = self.cosbeta_spline.eval(mf_min);
        self.set_alpha_ref(wf.mf_ref)?;
        let gamma_spline = self.gamma_spline.as_ref().expect("gamma spline");

        let alpha_offset = -self.alpha_ref + prec.alpha0;
        let mut angles = EulerAngles::with_capacity(n);
        angles.alpha.push(alphamin - self.alpha_ref + prec.alpha0);
        angles.cos_beta.push(cosbetamin);
        angles.gamma.push(0.);

        let mut hint = 0;
        let mut gamma_hint = 0;
        let mut step_hint = 0;
        for i in 1..n {
            let mf = hz_to_mf(frequencies[i], wf.m_tot);
            if mf > mf_cut {
                break;
            }
            if fast && self.continues_analytically() && mf > self.ftrans_mrd {
                // the merger-ringdown up to the cutoff
                let end =
                    i + frequencies[i..].partition_point(|&f| hz_to_mf(f, wf.m_tot) <= mf_cut);
                let mfs: Vec<f64> = frequencies[i..end]
                    .iter()
                    .map(|&f| hz_to_mf(f, wf.m_tot))
                    .collect();
                angles.alpha.extend(
                    mfs.iter()
                        .map(|&mf| self.alpha_mrd.alpha(mf) + alpha_offset),
                );
                let start = angles.cos_beta.len();
                angles.cos_beta.resize(start + mfs.len(), 0.0);
                self.beta_mrd
                    .cos_beta_fast(&mfs, &mut angles.cos_beta[start..], false);
                let delta_mf =
                    |k: usize| hz_to_mf(frequencies[i + k] - frequencies[i + k - 1], wf.m_tot);
                let gamma_start = angles.gamma[i - 1];
                self.push_mrd_gamma_fast(&mfs, delta_mf, false, gamma_start, &mut angles.gamma);
                break;
            }
            if mf < self.ftrans_mrd {
                let gamma_ok = gamma_spline.contains(mf);
                if !(self.alpha_spline.contains(mf) && gamma_ok) {
                    return Err(IMRPhenomXError::SpinTaylorAngles(
                        "interpolation of the SpinTaylor angles out of range",
                    ));
                }
                hint = self.alpha_spline.locate(mf, hint);
                gamma_hint = gamma_spline.locate(mf, gamma_hint);
                angles
                    .alpha
                    .push(self.alpha_spline.eval_in(hint, mf) + alpha_offset);
                angles.cos_beta.push(self.cosbeta_spline.eval_in(hint, mf));
                angles.gamma.push(gamma_spline.eval_in(gamma_hint, mf));
            } else if self.continues_analytically() {
                let delta_mf = hz_to_mf(frequencies[i] - frequencies[i - 1], wf.m_tot);
                angles.alpha.push(self.alpha_mrd.alpha(mf) + alpha_offset);
                angles.cos_beta.push(self.beta_mrd.beta(mf).cos());
                let previous = angles.gamma[i - 1];
                // A failed step keeps the previous value.
                let deltagamma = self.gamma_step(mf, delta_mf, &mut step_hint).unwrap_or(0.0);
                angles.gamma.push(previous + deltagamma);
            } else {
                // continue the angles into the ringdown with constant values
                angles.alpha.push(angles.alpha[i - 1]);
                angles.cos_beta.push(angles.cos_beta[i - 1]);
                angles.gamma.push(angles.gamma[i - 1]);
            }
        }
        Ok(angles)
    }
}

/// The integrand of the minimal rotation condition at the last point of a
/// Boole step of `gamma_from_alpha_cosbeta`.
#[derive(Clone, Copy, Debug, Default)]
struct BooleCache {
    mf: f64,
    inspiral: bool,
    value: f64,
    valid: bool,
}

impl BooleCache {
    /// The cached integrand if it was evaluated at exactly `mf` in the same
    /// branch (spline or merger-ringdown continuation).
    #[inline]
    fn lookup(&self, mf: f64, inspiral: bool) -> Option<f64> {
        (self.valid && self.inspiral == inspiral && self.mf.to_bits() == mf.to_bits())
            .then_some(self.value)
    }
}

/// Euler angles on a frequency grid.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct EulerAngles {
    /// Precession angle `alpha` of L about J.
    pub alpha: Vec<f64>,
    /// Cosine of the opening angle `beta` between L and J.
    pub cos_beta: Vec<f64>,
    /// Third Euler angle `gamma` (minimal rotation condition).
    pub gamma: Vec<f64>,
}

impl EulerAngles {
    fn with_capacity(n: usize) -> Self {
        Self {
            alpha: Vec::with_capacity(n),
            cos_beta: Vec::with_capacity(n),
            gamma: Vec::with_capacity(n),
        }
    }
}

/// `XLALSimIMRPhenomXUnwrapArray`, overwriting its input.
pub(crate) fn unwrap_in_place(values: &mut [f64]) {
    let Some(&first) = values.first() else {
        return;
    };
    let mut previous_input = first;
    for i in 1..values.len() {
        let input = values[i];
        let mut d = input - previous_input;
        d = if d > PI {
            d - 2. * PI
        } else if d < -PI {
            d + 2. * PI
        } else {
            d
        };
        values[i] = values[i - 1] + d;
        previous_input = input;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::MSUN_SI;
    use crate::imrphenomx::precession::PrecessionVersion;
    use crate::spintaylor::SpinTaylorOrbit;

    /// `IMRPhenomX_InspiralAngles_SpinTaylor` as LAL writes it: every series
    /// of each run, joined with `appendTS` and stitched to the fine run.
    fn reference_pn_arrays(wf: &WaveformStruct, prec: &PrecessionStruct, f_min: f64) -> PnArrays {
        let series = |o: &SpinTaylorOrbit| -> [Vec<f64>; 13] {
            [
                &o.v, &o.s1_x, &o.s1_y, &o.s1_z, &o.s2_x, &o.s2_y, &o.s2_z, &o.lnhat_x, &o.lnhat_y,
                &o.lnhat_z, &o.e1_x, &o.e1_y, &o.e1_z,
            ]
            .map(|s| s.clone())
        };
        let append_ts = |start: &[f64], end: &[f64]| {
            let mut out = start[..start.len() - 2].to_vec();
            out.extend_from_slice(end);
            out.push(0.0);
            out
        };
        let pi_gm = PI * (wf.m1_si + wf.m2_si) * (G_SI / C_SI) / (C_SI * C_SI);
        let binary = SpinTaylorT4::new(wf.m1_si, wf.m2_si, prec.chi1, prec.chi2);
        let f_cut = mf_to_hz(prec.f_ring_aligned + 8. * prec.f_damp_aligned, wf.m_tot);
        let dt = 0.5 * (10.0f32 as f64) / f_cut;
        let coarse = if (wf.f_ref - f_min).abs() < REAL4_EPS {
            series(&binary.evolve_orbit(dt, f_min, f_cut).unwrap())
        } else {
            let backward = series(&binary.evolve_orbit(dt, wf.f_ref, f_min - 0.5).unwrap());
            let forward = series(&binary.evolve_orbit(dt, wf.f_ref, f_cut).unwrap());
            std::array::from_fn(|k| append_ts(&backward[k], &forward[k]))
        };
        let len_low = coarse[0].len();
        let copy_length = len_low - 1 - 9.min(len_low - 1);
        let at = |k: usize| coarse[k][copy_length];
        let fine = SpinTaylorT4::new(
            wf.m1_si,
            wf.m2_si,
            [at(1), at(2), at(3)],
            [at(4), at(5), at(6)],
        )
        .with_orientation([at(7), at(8), at(9)], [at(10), at(11), at(12)]);
        let fine = series(
            &fine
                .evolve_orbit(0.5 / f_cut, at(0).powf(3.) / pi_gm, f_cut)
                .unwrap(),
        );
        let stitch = |k: usize| [&coarse[k][..copy_length], &fine[k][..]].concat();
        let last = |k: usize| *fine[k].last().unwrap();
        PnArrays {
            v: stitch(0),
            lnhat: [stitch(7), stitch(8), stitch(9)],
            s1_end: [last(1), last(2), last(3)],
            s2_end: [last(4), last(5), last(6)],
        }
    }

    /// The PN arrays built from the compact runs are LAL's, for a backward
    /// and a forward coarse run and for a single forward one.
    #[test]
    fn pn_arrays_match_the_full_series() {
        let (m1, m2) = (1.8, 1.3);
        let (chi1, chi2) = ([0.3, 0.1, 0.05], [-0.1, 0.2, 0.02]);
        let mut wf = WaveformStruct::new(
            m1, m2, chi1[2], chi2[2], 0.125, 30.0, 0.4, 25.0, 1024.0, 1.0, 0.0,
        )
        .unwrap();
        let prec = PrecessionStruct::new(
            &mut wf,
            m1 * MSUN_SI,
            m2 * MSUN_SI,
            chi1,
            chi2,
            PrecessionVersion::SpinTaylor,
            2,
            2,
        )
        .unwrap();
        for f_min in [25.0, wf.f_ref] {
            let arrays = inspiral_angles_spintaylor(&wf, &prec, f_min).unwrap();
            let reference = reference_pn_arrays(&wf, &prec, f_min);
            assert!(arrays.len() > 1000);
            assert_eq!(arrays.v, reference.v);
            assert_eq!(arrays.lnhat, reference.lnhat);
            assert_eq!(arrays.s1_end, reference.s1_end);
            assert_eq!(arrays.s2_end, reference.s2_end);
        }
    }

    /// A backward and a forward run recorded in [`CompactSamples`] give the
    /// same `v`, `LNhat` and full states as `SpinTaylorT4::evolve_orbit`,
    /// including states that were not kept (which are integrated again).
    #[test]
    fn compact_runs_match_the_orbit() {
        let binary = SpinTaylorT4::new(
            1.8 * MSUN_SI,
            1.3 * MSUN_SI,
            [0.3, 0.1, 0.05],
            [-0.1, 0.2, 0.02],
        );
        for (f_start, f_end) in [(30.0, 25.0), (30.0, 400.0)] {
            let delta_t = 1.0 / 2048.0;
            let orbit = binary.evolve_orbit(delta_t, f_start, f_end).unwrap();
            let run = Run::evolve(&binary, delta_t, f_start, f_end).unwrap();
            let n = run.len();
            assert_eq!(n, orbit.len());
            assert!(n > 4 * KEPT_STATES);
            let mut v = Vec::new();
            let mut lnhat = [Vec::new(), Vec::new(), Vec::new()];
            run.extend_series(&mut v, &mut lnhat, n);
            assert_eq!(v, orbit.v);
            assert_eq!(lnhat[0], orbit.lnhat_x);
            assert_eq!(lnhat[1], orbit.lnhat_y);
            assert_eq!(lnhat[2], orbit.lnhat_z);
            for j in [0, 1, KEPT_STATES, n / 2, n - KEPT_STATES - 1, n - 2, n - 1] {
                let point = run.point(j);
                assert_eq!(point.v, orbit.v[j]);
                assert_eq!(point.s1, [orbit.s1_x[j], orbit.s1_y[j], orbit.s1_z[j]]);
                assert_eq!(point.s2, [orbit.s2_x[j], orbit.s2_y[j], orbit.s2_z[j]]);
                let lnhat_j = [orbit.lnhat_x[j], orbit.lnhat_y[j], orbit.lnhat_z[j]];
                assert_eq!(point.lnhat, lnhat_j);
                assert_eq!(point.e1, [orbit.e1_x[j], orbit.e1_y[j], orbit.e1_z[j]]);
            }
            let (series_v, series_lnhat) = run.into_series(n - 3, n);
            assert_eq!(series_v, orbit.v[..n - 3]);
            assert_eq!(series_lnhat[2], orbit.lnhat_z[..n - 3]);
        }
    }

    #[test]
    fn unwrap_removes_jumps_of_two_pi() {
        let mut values = [3.0, -3.0, -0.5, 3.1, -3.1];
        unwrap_in_place(&mut values);
        let expected = [3.0, -3.0 + 2. * PI, -0.5 + 2. * PI, 3.1, -3.1 + 2. * PI];
        for (a, b) in values.iter().zip(expected) {
            assert!((a - b).abs() < 1e-14, "{a} {b}");
        }
    }
}
