//! Evaluation of IMRPhenomXP and IMRPhenomXPHM on a prepared
//! [`crate::grid::FrequencyGrid`], in passes over the grid.
//!
//! The per-frequency drivers evaluate, at each frequency, the Euler angles
//! for every `m'`, the twisting-up factors and every co-precessing mode. On
//! a grid the same quantities are computed one at a time over all
//! frequencies:
//!
//! 1. the Euler angles of each `m'` into arrays ([`AngleArrays`]), reusing
//!    the angles of another `m'` where the argument of the angle functions
//!    is bit-for-bit the same ([`angles_with_reuse`]);
//! 2. `exp(i alpha)` and the half-angle cosine and sine of `beta`
//!    ([`TwistArrays`]), with the vectorised `sin_cos`;
//! 3. each co-precessing mode, with the powers of `Mf` from the grid and
//!    `exp(-i m' epsilon)` folded into its phase, converted to complex
//!    values in one vectorised pass, and twisted up into the polarisations
//!    ([`add_twisted_mode`]).
//!
//! The angles are those of the per-frequency path (the same functions of
//! the same arguments); the modes and the twisting up agree with it to
//! rounding.

use std::f64::consts::PI;

use num_complex::Complex;

use crate::fastmath;
use crate::grid::{FrequencyGrid, ScalePowers};
use crate::imrphenomx::internals::UsefulPowers;
use crate::imrphenomx::precession::PrecessionStruct;
use crate::imrphenomx::wigner::half_angles;
use crate::imrphenomx::xhm::ModeGenerator;
use crate::imrphenomx::xphm::wigner_d;

/// The Euler angles `alpha`, `cos(beta)` and `epsilon` (including the
/// reference offsets) of one co-precessing `m'` at the grid points below
/// the cutoff.
#[derive(Clone, Debug, Default)]
pub(crate) struct AngleArrays {
    pub(crate) alpha: Vec<f64>,
    pub(crate) cos_beta: Vec<f64>,
    pub(crate) epsilon: Vec<f64>,
}

impl AngleArrays {
    pub(crate) fn with_capacity(n: usize) -> Self {
        Self {
            alpha: Vec::with_capacity(n),
            cos_beta: Vec::with_capacity(n),
            epsilon: Vec::with_capacity(n),
        }
    }

    #[inline]
    pub(crate) fn push(&mut self, alpha: f64, cos_beta: f64, epsilon: f64) {
        self.alpha.push(alpha);
        self.cos_beta.push(cos_beta);
        self.epsilon.push(epsilon);
    }

    /// The NNLO angles of `m'` at `omega = pi M f 2 / m'` for the first `n`
    /// frequencies `f` of `grid` (`m_sec` is `M` in seconds).
    pub(crate) fn nnlo(
        prec: &PrecessionStruct,
        eta: f64,
        grid: &FrequencyGrid,
        m_sec: f64,
        mprime: u32,
        n: usize,
    ) -> Self {
        let mut angles = Self {
            alpha: vec![0.0; n],
            cos_beta: vec![0.0; n],
            epsilon: vec![0.0; n],
        };
        prec.nnlo_angles_on_grid(
            eta,
            grid,
            PI * m_sec * 2.0 / mprime as f64,
            &mut angles.alpha,
            &mut angles.epsilon,
            &mut angles.cos_beta,
        );
        angles
    }
}

/// The Euler angles of the `m'` with the angle-function arguments
/// `arguments(m', i)` at the first `n` grid points (each non-decreasing in
/// `i`), evaluated by `evaluate(arguments, angles)` (which pushes the
/// angles of each argument). The angles of the first `m'` in `mprimes` are
/// evaluated at every point; those of the other `m'` are copied from it
/// wherever the argument is the same floating-point number (on a uniform
/// grid `2 Mf / m'` often coincides with the `Mf` of another grid point),
/// and evaluated elsewhere.
pub(crate) fn angles_with_reuse(
    mprimes: &[u32],
    n: usize,
    arguments: impl Fn(u32, usize) -> f64,
    mut evaluate: impl FnMut(&[f64], &mut AngleArrays),
) -> [Option<AngleArrays>; 5] {
    let mut out: [Option<AngleArrays>; 5] = Default::default();
    let Some((&first, others)) = mprimes.split_first() else {
        return out;
    };
    let reference_arguments: Vec<f64> = (0..n).map(|i| arguments(first, i)).collect();
    let mut reference = AngleArrays::with_capacity(n);
    evaluate(&reference_arguments, &mut reference);

    let mut pending_arguments = Vec::with_capacity(n);
    let mut pending_index = Vec::with_capacity(n);
    let mut pending = AngleArrays::with_capacity(n);
    for &mprime in others {
        let mut angles = AngleArrays {
            alpha: vec![0.0; n],
            cos_beta: vec![0.0; n],
            epsilon: vec![0.0; n],
        };
        pending_arguments.clear();
        pending_index.clear();
        let mut j = 0;
        for i in 0..n {
            let x = arguments(mprime, i);
            while j < n && reference_arguments[j] < x {
                j += 1;
            }
            if j < n && reference_arguments[j] == x {
                angles.alpha[i] = reference.alpha[j];
                angles.cos_beta[i] = reference.cos_beta[j];
                angles.epsilon[i] = reference.epsilon[j];
            } else {
                pending_arguments.push(x);
                pending_index.push(i);
            }
        }
        pending.alpha.clear();
        pending.cos_beta.clear();
        pending.epsilon.clear();
        evaluate(&pending_arguments, &mut pending);
        for (k, &i) in pending_index.iter().enumerate() {
            angles.alpha[i] = pending.alpha[k];
            angles.cos_beta[i] = pending.cos_beta[k];
            angles.epsilon[i] = pending.epsilon[k];
        }
        out[mprime as usize] = Some(angles);
    }
    out[first as usize] = Some(reference);
    out
}

/// The twisting-up factors of one `m'` that do not depend on the mode:
/// `exp(i alpha)` and `(cos(beta / 2), sin(beta / 2))`.
#[derive(Clone, Debug, Default)]
pub(crate) struct TwistArrays {
    rotation: Vec<Complex<f64>>,
    half_angles: Vec<[f64; 2]>,
}

impl TwistArrays {
    pub(crate) fn new(angles: &AngleArrays) -> Self {
        // exp(i alpha) = 1 exp(-i (-alpha))
        let mut rotation: Vec<Complex<f64>> = angles
            .alpha
            .iter()
            .map(|&alpha| Complex::new(1.0, -alpha))
            .collect();
        fastmath::polar_to_complex_conj(&mut rotation);
        let half_angles = angles
            .cos_beta
            .iter()
            .map(|&cos_beta| {
                let (c, s) = half_angles(cos_beta);
                [c, s]
            })
            .collect();
        Self {
            rotation,
            half_angles,
        }
    }
}

/// Fill `scratch` for [`fastmath::polar_to_complex_conj`] with the
/// co-precessing mode `exp(-i m' epsilon) h / 2` at the grid points
/// `offset + i`, where `h = A exp(i phi)` with `(A, phi) =
/// amplitude_phase(offset + i)`, stored as `(A / 2, m' epsilon - phi)`.
fn fill_polar(
    scratch: &mut [Complex<f64>],
    offset: usize,
    mprime: u32,
    epsilon: &[f64],
    amplitude_phase: impl Fn(usize) -> (f64, f64),
) {
    let m = mprime as f64;
    let epsilon = &epsilon[offset..offset + scratch.len()];
    for (i, (value, &epsilon)) in scratch.iter_mut().zip(epsilon).enumerate() {
        let (amplitude, phase) = amplitude_phase(offset + i);
        *value = Complex::new(0.5 * amplitude, m * epsilon - phase);
    }
}

/// `exp(-i m' epsilon) h / 2` at the first `scratch.len()` grid points,
/// written to `scratch`, for the mode `h = A exp(i phi)` with `(A, phi) =
/// amplitude_phase(i)` at grid point `i`.
pub(crate) fn coprecessing_mode(
    scratch: &mut [Complex<f64>],
    mprime: u32,
    epsilon: &[f64],
    amplitude_phase: impl Fn(usize) -> (f64, f64),
) {
    fill_polar(scratch, 0, mprime, epsilon, amplitude_phase);
    fastmath::polar_to_complex_conj(scratch);
}

/// [`coprecessing_mode`] for the IMRPhenomXHM mode `generator` with the
/// powers of `Mf = scale f` from `grid`. In the ringdown of the (3,2) mode
/// with mode mixing the two terms of the spheroidal-to-spherical rotation
/// are converted separately (`extra` holds the second) and added, rather
/// than going through the amplitude and phase of their sum.
pub(crate) fn coprecessing_mode_generator(
    scratch: &mut [Complex<f64>],
    extra: &mut Vec<Complex<f64>>,
    mprime: u32,
    epsilon: &[f64],
    generator: &ModeGenerator,
    grid: &FrequencyGrid,
    scale: &ScalePowers,
) {
    let n = scratch.len();
    let amplitude_phase = |i: usize| {
        let powers = grid.powers(scale, i);
        generator.amplitude_phase(powers.itself, &UsefulPowers::from_grid(&powers))
    };
    let Some((mode, amp0, start)) = generator.mixing_ringdown() else {
        return coprecessing_mode(scratch, mprime, epsilon, amplitude_phase);
    };
    let start = grid.frequencies()[..n].partition_point(|&f| scale.value * f < start);
    let (inspiral, ringdown) = scratch.split_at_mut(start);
    fill_polar(inspiral, 0, mprime, epsilon, amplitude_phase);
    extra.clear();
    extra.resize(ringdown.len(), Complex::new(0.0, 0.0));
    let m = mprime as f64;
    for (k, (first, second)) in ringdown.iter_mut().zip(extra.iter_mut()).enumerate() {
        let i = start + k;
        let powers = grid.powers(scale, i);
        let (lin, [(a1, phi1), (a2, phi2)]) =
            mode.mixing_ringdown_terms(&UsefulPowers::from_grid(&powers));
        let shift = m * epsilon[i] - lin;
        *first = Complex::new(0.5 * amp0 * a1, shift - phi1);
        *second = Complex::new(0.5 * amp0 * a2, shift - phi2);
    }
    fastmath::polar_to_complex_conj(scratch);
    fastmath::polar_to_complex_conj(extra);
    let [c1, c2] = mode.mixing_factors();
    for (h, second) in scratch[start..].iter_mut().zip(extra.iter()) {
        *h = c1 * *h + c2 * second;
    }
}

/// Add the co-precessing mode `(l, -m')`, `h = exp(-i m' epsilon) h_{l,-m'}
/// / 2` from [`coprecessing_mode`], to the polarisations: twisted up with
/// the Wigner-d matrices and the spherical harmonics `ylm` as in
/// `IMRPhenomXPHMTwistUp` (`hplus`, `hcross` at least as long as `h`).
pub(crate) fn add_twisted_mode(
    ell: u32,
    mprime: u32,
    ylm: &[f64],
    twist: &TwistArrays,
    h: &[Complex<f64>],
    hplus: &mut [Complex<f64>],
    hcross: &mut [Complex<f64>],
) {
    let n = h.len();
    let rotation = &twist.rotation[..n];
    let half = &twist.half_angles[..n];
    let (hplus, hcross) = (&mut hplus[..n], &mut hcross[..n]);
    match (ell, mprime) {
        (2, 2) => twist_up::<2, 2>(ylm, h, rotation, half, hplus, hcross),
        (2, 1) => twist_up::<2, 1>(ylm, h, rotation, half, hplus, hcross),
        (3, 3) => twist_up::<3, 3>(ylm, h, rotation, half, hplus, hcross),
        (3, 2) => twist_up::<3, 2>(ylm, h, rotation, half, hplus, hcross),
        (4, 4) => twist_up::<4, 4>(ylm, h, rotation, half, hplus, hcross),
        _ => unreachable!("mode ({ell},{mprime}) is not twisted up"),
    }
}

/// `IMRPhenomXPHMTwistUp` for mode `(L, M)` at every grid point, with the
/// sums over `k` of `crate::imrphenomx::xphm` grouped by `|k|`.
#[inline(always)]
fn twist_up<const L: usize, const M: usize>(
    ylm: &[f64],
    h: &[Complex<f64>],
    rotation: &[Complex<f64>],
    half: &[[f64; 2]],
    hplus: &mut [Complex<f64>],
    hcross: &mut [Complex<f64>],
) {
    let sigma = if L % 2 == 0 { 1.0 } else { -1.0 };
    // a_k = (-1)^(k + m') d_{-k} Y_k, b_k = sigma d_k Y_k with k = index - L
    let mut ya = [0.0; 9];
    let mut yb = [0.0; 9];
    for index in 0..=2 * L {
        let parity = (index + M + L) % 2;
        ya[index] = if parity == 0 { ylm[index] } else { -ylm[index] };
        yb[index] = sigma * ylm[index];
    }
    let n = h.len();
    let (rotation, half) = (&rotation[..n], &half[..n]);
    let (hplus, hcross) = (&mut hplus[..n], &mut hcross[..n]);
    // Only constant indices below, so that the arrays live in registers.
    for i in 0..n {
        // exp(i k alpha) by multiplication, as in LAL
        let e1 = rotation[i];
        let e2 = e1 * e1;
        let e3 = e1 * e2;
        let e4 = e1 * e3;
        let [cb, sb] = half[i];
        let d = wigner_d(L as u32, M as u32, &powers(cb), &powers(sb));
        // (a + b, b - a) at index L + k
        let pq = |k: isize| {
            let index = (L as isize + k) as usize;
            let a = ya[index] * d[2 * L - index];
            let b = yb[index] * d[index];
            (a + b, b - a)
        };
        let (p0, q0) = pq(0);
        let mut hp_re = p0;
        let mut hp_im = 0.0;
        let mut hc_re = 0.0;
        let mut hc_im = -q0;
        let mut add = |k: isize, e: Complex<f64>| {
            let (p_plus, q_plus) = pq(k);
            let (p_minus, q_minus) = pq(-k);
            hp_re += (p_plus + p_minus) * e.re;
            hc_im -= (q_plus + q_minus) * e.re;
            hp_im += (q_plus - q_minus) * e.im;
            hc_re += (p_plus - p_minus) * e.im;
        };
        add(1, e1);
        add(2, e2);
        if L >= 3 {
            add(3, e3);
        }
        if L >= 4 {
            add(4, e4);
        }
        let hi = h[i];
        hplus[i] += hi * Complex::new(hp_re, hp_im);
        hcross[i] += hi * Complex::new(hc_re, hc_im);
    }
}

/// `1, x, ..., x^8` by successive multiplication (as `cBetah2 = cBetah *
/// cBetah`, ... in LAL).
#[inline(always)]
fn powers(x: f64) -> [f64; 9] {
    let x2 = x * x;
    let x3 = x * x2;
    let x4 = x * x3;
    let x5 = x * x4;
    let x6 = x * x5;
    let x7 = x * x6;
    let x8 = x * x7;
    [1.0, x, x2, x3, x4, x5, x6, x7, x8]
}

/// Rotate the polarisations by `2 zeta` (the end of
/// `IMRPhenomXPHM_hplushcross` and `IMRPhenomXPGenerateFD`).
pub(crate) fn rotate_polarizations(
    zeta: f64,
    hplus: &mut [Complex<f64>],
    hcross: &mut [Complex<f64>],
) {
    if zeta.abs() > 0.0 {
        let cos_pol = (2.0 * zeta).cos();
        let sin_pol = (2.0 * zeta).sin();
        for (p, c) in hplus.iter_mut().zip(hcross.iter_mut()) {
            let (pp, pc) = (*p, *c);
            *p = pp * cos_pol + pc * sin_pol;
            *c = pc * cos_pol - pp * sin_pol;
        }
    }
}
