// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! Rotation of the waveform modes from the coorbital to the inertial frame
//! (`TransformModesCoorbitalToInertial` and `WignerDMatrices_Compute` in
//! `LALSimNRSurrogateUtilities.c`).
//!
//! LAL computes every intermediate quantity as a separate pass of GSL
//! vector operations over all coorbital times. Each time sample is
//! independent, so this implementation processes the samples in small
//! blocks with all intermediate values kept in registers or L1 cache,
//! while performing the same floating-point operations (each GSL vector
//! operation rounds, so no operations are fused here) in the same order.

use crate::nrsur::data::{mode_index, LMAX, N_MODES};
use crate::nrsur::simd::F64x4;

/// `|ra|^2` or `|rb|^2` below which LAL treats the quaternion component as
/// zero.
const EPS_SQR: f64 = 1.0e-24;

fn factorial(n: i32) -> f64 {
    if n <= 0 {
        1.0
    } else {
        factorial(n - 1) * n as f64
    }
}

fn factorial_ratio(n: i32, k: i32) -> f64 {
    if n <= k {
        1.0
    } else {
        factorial_ratio(n - 1, k) * n as f64
    }
}

fn binomial(n: i32, k: i32) -> f64 {
    factorial_ratio(n, k) / factorial(n - k)
}

/// One element `D^ell_{m, mp}` of the Wigner-D matrices in LAL's
/// parametrisation.
#[derive(Clone, Debug)]
struct WignerElement {
    ell: i32,
    m: i32,
    mp: i32,
    /// `sqrt((ell+m)!(ell-m)! / ((ell+mp)!(ell-mp)!))` (`wigner_coef`).
    coef: f64,
    /// Index of `ra^(m+mp)` and `rb^(m-mp)` in the complex power tables.
    ra_power: usize,
    rb_power: usize,
    /// Index of `ra^(m+mp) rb^(m-mp)` in the products of `(m, mp)`.
    product: usize,
    /// Exponent `ell - m` of `|ra|^2`.
    abs_ra_power: usize,
    /// First `rho`, the number of terms and the coefficients `c(rho) =
    /// (-1)^rho binomial(ell+mp, rho) binomial(ell-mp, ell-rho-m)`.
    rho_min: usize,
    n_rho: usize,
    c: [f64; LMAX + 1],
}

/// The parameter-independent parts of the Wigner-D matrices for
/// `2 <= ell <= lmax`, in LAL's order (`WignerDMatrix_Index`).
#[derive(Clone, Debug)]
pub(crate) struct WignerTables {
    lmax: i32,
    elements: Vec<WignerElement>,
    /// The elements whose sum over `rho` has `r + 1` terms, for each `r`.
    by_n_rho: [Vec<usize>; LMAX + 1],
}

impl WignerTables {
    pub(crate) fn new(lmax: usize) -> Self {
        let lmax = lmax as i32;
        let mut elements = Vec::new();
        for ell in 2..=lmax {
            for m in -ell..=ell {
                for mp in -ell..=ell {
                    let coef = (factorial(ell + m) * factorial(ell - m)
                        / (factorial(ell + mp) * factorial(ell - mp)))
                    .sqrt();
                    let rho_min = if mp > m { mp - m } else { 0 };
                    let rho_max = if mp > -m { ell - m } else { ell + mp };
                    let mut c = [0.0; LMAX + 1];
                    for rho in rho_min..=rho_max {
                        let value = binomial(ell + mp, rho) * binomial(ell - mp, ell - rho - m);
                        c[(rho - rho_min) as usize] = if rho % 2 == 1 { -value } else { value };
                    }
                    elements.push(WignerElement {
                        ell,
                        m,
                        mp,
                        coef,
                        ra_power: (2 * lmax + m + mp) as usize,
                        rb_power: (2 * lmax + m - mp) as usize,
                        product: ((m + lmax) * (2 * lmax + 1) + mp + lmax) as usize,
                        abs_ra_power: (ell - m) as usize,
                        rho_min: rho_min as usize,
                        n_rho: (rho_max - rho_min + 1) as usize,
                        c,
                    });
                }
            }
        }
        let mut by_n_rho: [Vec<usize>; LMAX + 1] = Default::default();
        for (k, e) in elements.iter().enumerate() {
            by_n_rho[e.n_rho - 1].push(k);
        }
        Self {
            lmax,
            elements,
            by_n_rho,
        }
    }
}

const N_COMPLEX_POWERS: usize = 4 * LMAX + 1;
const N_REAL_POWERS: usize = 2 * LMAX + 1;

/// `ComplexPowers_Compute` for one sample: `z^n` for `-2 lmax <= n <= 2
/// lmax`, stored at `n + 2 lmax`.
#[inline]
fn complex_powers(x: f64, y: f64) -> ([f64; N_COMPLEX_POWERS], [f64; N_COMPLEX_POWERS]) {
    let zero = 2 * LMAX;
    let mut re = [0.0; N_COMPLEX_POWERS];
    let mut im = [0.0; N_COMPLEX_POWERS];
    re[zero] = 1.0;
    re[zero + 1] = x;
    im[zero + 1] = y;
    for i in zero + 1..N_COMPLEX_POWERS - 1 {
        re[i + 1] = re[i] * x - im[i] * y;
        im[i + 1] = re[i] * y + im[i] * x;
    }
    for power in 1..=2 * LMAX {
        let i = zero + power;
        let j = zero - power;
        let mag = re[i] * re[i] + im[i] * im[i];
        re[j] = re[i] / mag;
        im[j] = -im[i] / mag;
    }
    (re, im)
}

/// `RealPowers_Compute` for one sample.
#[inline]
fn real_powers(x: f64) -> [f64; N_REAL_POWERS] {
    let mut p = [0.0; N_REAL_POWERS];
    p[0] = 1.0;
    p[1] = x;
    for n in 2..N_REAL_POWERS {
        p[n] = p[n - 1] * x;
    }
    p
}

/// The Wigner-D matrices for the quaternion `quat` of one sample, in the
/// order of [`WignerTables`] (`WignerDMatrices_Compute`).
fn wigner_d(tables: &WignerTables, quat: [f64; 4], d_re: &mut [f64], d_im: &mut [f64]) {
    use crate::spintaylor::fp::fma;
    // ra = q0 - i q3, rb = -q2 - i q1
    let mut index_type = 0;
    let rb_abs_sqr = fma(quat[2], quat[2], quat[1] * quat[1]);
    let (rb_re, rb_im, rb_mag_sqr) = if rb_abs_sqr < EPS_SQR {
        index_type = 2;
        (0.25, 0.25, 0.5)
    } else {
        (-quat[2], -quat[1], rb_abs_sqr)
    };
    let ra_abs_sqr = fma(quat[0], quat[0], quat[3] * quat[3]);
    let (ra_re, ra_im, ra_mag_sqr) = if ra_abs_sqr < EPS_SQR {
        index_type = 1;
        (0.25, 0.25, 0.5)
    } else {
        (quat[0], -quat[3], ra_abs_sqr)
    };
    let ratio = rb_mag_sqr / ra_mag_sqr;
    let (ra_pow_re, ra_pow_im) = complex_powers(ra_re, ra_im);
    let (rb_pow_re, rb_pow_im) = complex_powers(rb_re, rb_im);
    let abs_ra_pow = real_powers(ra_mag_sqr);
    let ratio_pow = real_powers(ratio);

    for (k, e) in tables.elements.iter().enumerate() {
        let (a_re, a_im) = (ra_pow_re[e.ra_power], ra_pow_im[e.ra_power]);
        let (b_re, b_im) = (rb_pow_re[e.rb_power], rb_pow_im[e.rb_power]);
        // complex_vector_mult, then scale by coef and |ra|^(2(ell-m))
        let mut re = a_re * b_re - a_im * b_im;
        let mut im = a_im * b_re + a_re * b_im;
        re *= e.coef;
        im *= e.coef;
        re *= abs_ra_pow[e.abs_ra_power];
        im *= abs_ra_pow[e.abs_ra_power];
        let mut sum = 0.0;
        for (i, &c) in e.c[..e.n_rho].iter().enumerate() {
            sum += ratio_pow[e.rho_min + i] * c;
        }
        d_re[k] = re * sum;
        d_im[k] = im * sum;
    }

    if index_type != 0 {
        // Fix the samples where ra or rb vanishes.
        for (k, e) in tables.elements.iter().enumerate() {
            let (mut zx, mut zy) = (0.0, 0.0);
            let fixed = if index_type == 1 && e.mp == -e.m {
                // (-1)^(ell+m+1) rb^(2m)
                Some((-quat[2], -quat[1], (e.ell + e.m) % 2 == 0))
            } else if index_type == 2 && e.mp == e.m {
                // ra^(2m)
                Some((quat[0], -quat[3], false))
            } else {
                None
            };
            if let Some((tmp_re, tmp_im, negate)) = fixed {
                zx = 1.0;
                zy = 0.0;
                for _ in 0..(2 * e.m).max(0) {
                    let c = zx * tmp_im;
                    zx = crate::spintaylor::fp::fms(zx, tmp_re, zy * tmp_im);
                    zy = fma(zy, tmp_re, c);
                }
                if negate {
                    zx = -zx;
                    zy = -zy;
                }
            }
            d_re[k] = zx;
            d_im[k] = zy;
        }
    }
}

/// Number of Wigner-D elements for `2 <= ell <= LMAX`.
pub(crate) const N_WIGNER: usize = {
    let mut n = 0;
    let mut ell = 2;
    while ell <= LMAX {
        n += (2 * ell + 1) * (2 * ell + 1);
        ell += 1;
    }
    n
};

/// Rotate the coorbital modes `h_coorb` (`[mode][time]`, real and imaginary
/// parts) to the inertial frame, writing `h` (same layout), given the
/// coprecessing-frame quaternions and the orbital phase at each time.
///
/// `cos_sin_phase[t]` holds `cos(m phi)`, `sin(m phi)` for `m = 1..=4`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn coorbital_to_inertial(
    tables: &WignerTables,
    n_times: usize,
    quat: [&[f64]; 4],
    cos_mphi: &[[f64; LMAX]],
    sin_mphi: &[[f64; LMAX]],
    h_coorb_re: &[f64],
    h_coorb_im: &[f64],
    h_re: &mut [f64],
    h_im: &mut [f64],
) {
    debug_assert_eq!(tables.lmax as usize, LMAX);
    let mut d_re = [0.0; N_WIGNER];
    let mut d_im = [0.0; N_WIGNER];
    let mut d_block = [(V::splat(0.0), V::splat(0.0)); N_WIGNER];
    let mut t = 0;
    while t < n_times {
        if t + LANES <= n_times {
            let q: [[f64; LANES]; 4] =
                std::array::from_fn(|k| std::array::from_fn(|l| quat[k][t + l]));
            if !(0..LANES).any(|l| is_special([q[0][l], q[1][l], q[2][l], q[3][l]])) {
                let mut input = BlockInput {
                    tables,
                    d: &mut d_block,
                    n_times,
                    t,
                    quat: q,
                    cos_mphi: &cos_mphi[t..t + LANES],
                    sin_mphi: &sin_mphi[t..t + LANES],
                    h_coorb_re,
                    h_coorb_im,
                };
                transform_block(&mut input, h_re, h_im);
                t += LANES;
                continue;
            }
        }
        let end = (t + LANES).min(n_times);
        for t in t..end {
            transform_sample(
                tables,
                n_times,
                t,
                [quat[0][t], quat[1][t], quat[2][t], quat[3][t]],
                cos_mphi[t],
                sin_mphi[t],
                h_coorb_re,
                h_coorb_im,
                h_re,
                h_im,
                &mut d_re,
                &mut d_im,
            );
        }
        t = end;
    }
}

/// Whether LAL treats `ra` or `rb` of this quaternion as zero.
#[inline]
fn is_special(quat: [f64; 4]) -> bool {
    use crate::spintaylor::fp::fma;
    fma(quat[2], quat[2], quat[1] * quat[1]) < EPS_SQR
        || fma(quat[0], quat[0], quat[3] * quat[3]) < EPS_SQR
}

/// `m` of each mode index.
const MODE_M: [i32; N_MODES] = {
    let mut out = [0; N_MODES];
    let mut ell = 2;
    let mut i = 0;
    while ell <= LMAX {
        let mut m = -(ell as i32);
        while m <= ell as i32 {
            out[i] = m;
            m += 1;
            i += 1;
        }
        ell += 1;
    }
    out
};

/// Number of time samples processed together by [`transform_block`].
const LANES: usize = F64x4::LANES;
type V = F64x4;

struct BlockInput<'a> {
    tables: &'a WignerTables,
    /// Scratch for the Wigner-D elements.
    d: &'a mut [(V, V); N_WIGNER],
    n_times: usize,
    t: usize,
    quat: [[f64; LANES]; 4],
    cos_mphi: &'a [[f64; LMAX]],
    sin_mphi: &'a [[f64; LMAX]],
    h_coorb_re: &'a [f64],
    h_coorb_im: &'a [f64],
}

/// [`transform_sample`] for `LANES` consecutive samples none of which is
/// special, with the same operations lane by lane. The product
/// `ra^(m+mp) rb^(m-mp)`, which does not depend on `ell`, is computed once
/// per `(m, mp)`.
fn transform_block(input: &mut BlockInput, h_re: &mut [f64], h_im: &mut [f64]) {
    use crate::spintaylor::fp::fma;
    let (n, t) = (input.n_times, input.t);
    let q = &input.quat;
    let ra_re = V::from_array(q[0]);
    let ra_im = -V::from_array(q[3]);
    let rb_re = -V::from_array(q[2]);
    let rb_im = -V::from_array(q[1]);
    let ra_mag = V::from_array(std::array::from_fn(|l| {
        fma(q[0][l], q[0][l], q[3][l] * q[3][l])
    }));
    let rb_mag = V::from_array(std::array::from_fn(|l| {
        fma(q[2][l], q[2][l], q[1][l] * q[1][l])
    }));
    let ratio = rb_mag / ra_mag;

    let (ra_pow_re, ra_pow_im) = complex_powers_lanes(ra_re, ra_im);
    let (rb_pow_re, rb_pow_im) = complex_powers_lanes(rb_re, rb_im);
    let abs_ra_pow = real_powers_lanes(ra_mag);
    let ratio_pow = real_powers_lanes(ratio);

    const NM: usize = 2 * LMAX + 1;
    let zero = 2 * LMAX;
    let prod: [(V, V); NM * NM] = std::array::from_fn(|index| {
        let (m, mp) = (index / NM, index % NM);
        // m + mp and m - mp relative to the zeroth power
        let a = m + mp;
        let b = zero + m - mp;
        let (a_re, a_im) = (ra_pow_re[a], ra_pow_im[a]);
        let (b_re, b_im) = (rb_pow_re[b], rb_pow_im[b]);
        (a_re * b_re - a_im * b_im, a_im * b_re + a_re * b_im)
    });

    // coprecessing modes
    let copr: [(V, V); N_MODES] = std::array::from_fn(|i| {
        let m = MODE_M[i];
        let re = V::from_slice(&input.h_coorb_re[i * n + t..]);
        let im = V::from_slice(&input.h_coorb_im[i * n + t..]);
        if m == 0 {
            (re, im)
        } else {
            let k = m.unsigned_abs() as usize - 1;
            let c = V::from_array(std::array::from_fn(|l| input.cos_mphi[l][k]));
            let s = V::from_array(std::array::from_fn(|l| input.sin_mphi[l][k]));
            let s = if m > 0 { s } else { -s };
            (re * c + im * s, -(re * s) + im * c)
        }
    });

    // The Wigner-D matrices: the elements are independent, so they are
    // evaluated grouped by the length of their sum over rho.
    let tables = input.tables;
    let d = &mut *input.d;
    wigner_group::<1>(tables, &prod, &abs_ra_pow, &ratio_pow, d);
    wigner_group::<2>(tables, &prod, &abs_ra_pow, &ratio_pow, d);
    wigner_group::<3>(tables, &prod, &abs_ra_pow, &ratio_pow, d);
    wigner_group::<4>(tables, &prod, &abs_ra_pow, &ratio_pow, d);
    wigner_group::<5>(tables, &prod, &abs_ra_pow, &ratio_pow, d);

    // h^{l,m} = sum_mp D^l_{m,mp} h_copr^{l,mp}, accumulated over mp in
    // order for all m of an ell at once (independent sums).
    let mut offset = 0;
    for ell in 2..=LMAX {
        let width = 2 * ell + 1;
        let base = mode_index(ell, -(ell as i32));
        let mut re = [V::splat(0.0); 2 * LMAX + 1];
        let mut im = [V::splat(0.0); 2 * LMAX + 1];
        for jp in 0..width {
            let (cr, ci) = copr[base + jp];
            for jm in 0..width {
                let (d_re, d_im) = input.d[offset + jm * width + jp];
                re[jm] = re[jm] + cr * d_re;
                im[jm] = im[jm] + cr * d_im;
                im[jm] = im[jm] + ci * d_re;
                re[jm] = re[jm] - ci * d_im;
            }
        }
        for jm in 0..width {
            let i = base + jm;
            h_re[i * n + t..i * n + t + LANES].copy_from_slice(&re[jm].to_array());
            h_im[i * n + t..i * n + t + LANES].copy_from_slice(&im[jm].to_array());
        }
        offset += width * width;
    }
}

/// The Wigner-D elements whose sum over `rho` has `R` terms.
#[inline(always)]
fn wigner_group<const R: usize>(
    tables: &WignerTables,
    prod: &[(V, V)],
    abs_ra_pow: &[V; N_REAL_POWERS],
    ratio_pow: &[V; N_REAL_POWERS],
    d: &mut [(V, V); N_WIGNER],
) {
    for &k in &tables.by_n_rho[R - 1] {
        let e = &tables.elements[k];
        let coef = V::splat(e.coef);
        let abs = abs_ra_pow[e.abs_ra_power];
        let mut sum = V::splat(0.0);
        for r in 0..R {
            sum = sum + ratio_pow[e.rho_min + r] * V::splat(e.c[r]);
        }
        let (p_re, p_im) = prod[e.product];
        d[k] = (p_re * coef * abs * sum, p_im * coef * abs * sum);
    }
}

/// [`complex_powers`] lane by lane.
#[inline(always)]
fn complex_powers_lanes(x: V, y: V) -> ([V; N_COMPLEX_POWERS], [V; N_COMPLEX_POWERS]) {
    let zero = 2 * LMAX;
    let mut re = [V::splat(0.0); N_COMPLEX_POWERS];
    let mut im = [V::splat(0.0); N_COMPLEX_POWERS];
    re[zero] = V::splat(1.0);
    re[zero + 1] = x;
    im[zero + 1] = y;
    for i in zero + 1..N_COMPLEX_POWERS - 1 {
        re[i + 1] = re[i] * x - im[i] * y;
        im[i + 1] = re[i] * y + im[i] * x;
    }
    for power in 1..=2 * LMAX {
        let i = zero + power;
        let j = zero - power;
        let mag = re[i] * re[i] + im[i] * im[i];
        re[j] = re[i] / mag;
        im[j] = -im[i] / mag;
    }
    (re, im)
}

/// [`real_powers`] lane by lane.
#[inline(always)]
fn real_powers_lanes(x: V) -> [V; N_REAL_POWERS] {
    let mut p = [V::splat(0.0); N_REAL_POWERS];
    p[0] = V::splat(1.0);
    p[1] = x;
    for n in 2..N_REAL_POWERS {
        p[n] = p[n - 1] * x;
    }
    p
}

/// The transformation of one sample (`TransformModesCoorbitalToInertial`
/// and `WignerDMatrices_Compute` restricted to time `t`), including LAL's
/// treatment of quaternions with vanishing `ra` or `rb`.
#[allow(clippy::too_many_arguments)]
fn transform_sample(
    tables: &WignerTables,
    n_times: usize,
    t: usize,
    quat: [f64; 4],
    cos_mphi: [f64; LMAX],
    sin_mphi: [f64; LMAX],
    h_coorb_re: &[f64],
    h_coorb_im: &[f64],
    h_re: &mut [f64],
    h_im: &mut [f64],
    d_re: &mut [f64; N_WIGNER],
    d_im: &mut [f64; N_WIGNER],
) {
    {
        // coprecessing modes: h_copr^{l,m} = exp(-i m phi) h_coorb^{l,m}
        let mut copr_re = [0.0; N_MODES];
        let mut copr_im = [0.0; N_MODES];
        for ell in 2..=LMAX {
            for m in -(ell as i32)..=ell as i32 {
                let i = mode_index(ell, m);
                let (re, im) = (h_coorb_re[i * n_times + t], h_coorb_im[i * n_times + t]);
                if m == 0 {
                    copr_re[i] = re;
                    copr_im[i] = im;
                } else {
                    let k = m.unsigned_abs() as usize - 1;
                    let c = cos_mphi[k];
                    let s = if m > 0 { sin_mphi[k] } else { -sin_mphi[k] };
                    copr_re[i] = re * c + im * s;
                    copr_im[i] = -(re * s) + im * c;
                }
            }
        }

        wigner_d(tables, quat, d_re, d_im);

        let mut k = 0;
        for ell in 2..=LMAX {
            let l = ell as i32;
            let base = mode_index(ell, -l);
            for m in -l..=l {
                let mut re = 0.0;
                let mut im = 0.0;
                for j in base..base + 2 * ell + 1 {
                    re += copr_re[j] * d_re[k];
                    im += copr_re[j] * d_im[k];
                    im += copr_im[j] * d_re[k];
                    re -= copr_im[j] * d_im[k];
                    k += 1;
                }
                let i = mode_index(ell, m);
                h_re[i * n_times + t] = re;
                h_im[i * n_times + t] = im;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_count() {
        assert_eq!(N_WIGNER, 25 + 49 + 81);
        assert_eq!(WignerTables::new(LMAX).elements.len(), N_WIGNER);
    }

    /// The Wigner-D matrices are unitary.
    #[test]
    fn unitary() {
        let tables = WignerTables::new(LMAX);
        let norm = (0.3f64 * 0.3 + 0.5 * 0.5 + 0.1 * 0.1 + 0.8 * 0.8).sqrt();
        let quat = [0.3 / norm, -0.5 / norm, 0.1 / norm, 0.8 / norm];
        let mut d_re = [0.0; N_WIGNER];
        let mut d_im = [0.0; N_WIGNER];
        wigner_d(&tables, quat, &mut d_re, &mut d_im);
        let mut offset = 0;
        for ell in 2..=LMAX {
            let n = 2 * ell + 1;
            for a in 0..n {
                for b in 0..n {
                    let mut re = 0.0;
                    let mut im = 0.0;
                    for c in 0..n {
                        let (x_re, x_im) = (d_re[offset + a * n + c], d_im[offset + a * n + c]);
                        let (y_re, y_im) = (d_re[offset + b * n + c], -d_im[offset + b * n + c]);
                        re += x_re * y_re - x_im * y_im;
                        im += x_re * y_im + x_im * y_re;
                    }
                    let expected = if a == b { 1.0 } else { 0.0 };
                    assert!((re - expected).abs() < 1e-12 && im.abs() < 1e-12);
                }
            }
            offset += n * n;
        }
    }
}
