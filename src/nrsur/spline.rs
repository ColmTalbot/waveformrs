// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Interpolation as done by `gsl_interp_cspline` and `gsl_interp_polynomial`,
//! specialised to the fixed knots of the surrogate.
//!
//! The natural cubic spline of GSL solves a symmetric tridiagonal system
//! whose matrix depends only on the knots. [`FixedKnots`] factorises it once
//! (when the surrogate data are loaded); fitting a spline to new values then
//! only needs the right-hand side and the two substitution sweeps, with the
//! same floating-point operations as `gsl_spline_init`, so the spline
//! values agree with GSL to rounding (see [`crate::interpolate`] for the
//! general spline and [`crate::spintaylor::fp`] for the fused multiply-adds).

use crate::spintaylor::fp::{fma, fms, fnms};

/// Knots of a natural cubic spline with the LDL^T factors of GSL's
/// tridiagonal system.
#[derive(Clone, Debug)]
pub(crate) struct FixedKnots {
    pub(crate) x: Vec<f64>,
    /// `x[i + 1] - x[i]`.
    h: Vec<f64>,
    /// `1 / h[i]` (zero for repeated knots, as in GSL).
    h_inv: Vec<f64>,
    /// Diagonal of the LDL^T factorisation (`alpha` in
    /// `gsl_linalg_solve_symm_tridiag`).
    alpha: Vec<f64>,
    /// Sub-diagonal of L (`gamma`).
    gamma: Vec<f64>,
}

impl FixedKnots {
    /// Factorise the spline system for strictly increasing knots `x` (at
    /// least four of them).
    pub(crate) fn new(x: &[f64]) -> Self {
        let size = x.len();
        assert!(size >= 4, "need at least four knots");
        let h: Vec<f64> = x.windows(2).map(|w| w[1] - w[0]).collect();
        let h_inv: Vec<f64> = h
            .iter()
            .map(|&h| if h != 0.0 { 1.0 / h } else { 0.0 })
            .collect();
        let sys_size = size - 2;
        let mut alpha = vec![0.0; sys_size];
        let mut gamma = vec![0.0; sys_size];
        let diag = |i: usize| 2.0 * (h[i + 1] + h[i]);
        let offdiag = |i: usize| h[i + 1];
        alpha[0] = diag(0);
        gamma[0] = offdiag(0) / alpha[0];
        for i in 1..sys_size - 1 {
            alpha[i] = fnms(offdiag(i - 1), gamma[i - 1], diag(i));
            gamma[i] = offdiag(i) / alpha[i];
        }
        alpha[sys_size - 1] = fnms(
            offdiag(sys_size - 2),
            gamma[sys_size - 2],
            diag(sys_size - 1),
        );
        Self {
            x: x.to_vec(),
            h,
            h_inv,
            alpha,
            gamma,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.x.len()
    }

    /// GSL's `c` coefficients (half the second derivatives) of the natural
    /// spline through `y`, written to `c` (length `len()`). `z` is scratch
    /// space of length `len()`.
    pub(crate) fn second_derivatives(&self, y: &[f64], c: &mut [f64], z: &mut [f64]) {
        let n = self.alpha.len();
        debug_assert!(y.len() == n + 2 && c.len() == n + 2 && z.len() >= n);
        let z = &mut z[..n];
        // right-hand side g of cspline_init, forward substitution
        for i in 0..n {
            let ydiff_i = y[i + 1] - y[i];
            let ydiff_ip1 = y[i + 2] - y[i + 1];
            let g = 3.0 * fms(ydiff_ip1, self.h_inv[i + 1], ydiff_i * self.h_inv[i]);
            z[i] = if i == 0 {
                g
            } else {
                fnms(self.gamma[i - 1], z[i - 1], g)
            };
        }
        for (z, alpha) in z.iter_mut().zip(&self.alpha) {
            *z /= alpha;
        }
        c[0] = 0.0;
        c[n + 1] = 0.0;
        c[n] = z[n - 1];
        for i in (0..n - 1).rev() {
            c[i + 1] = fnms(self.gamma[i], c[i + 2], z[i]);
        }
    }

    /// The polynomial coefficients `(b, d)` of interval `i` (`a = y[i]`,
    /// `c = c[i]`), as `coeff_calc` in GSL's `cspline.c`.
    #[inline(always)]
    pub(crate) fn interval(&self, i: usize, y: &[f64], c: &[f64]) -> (f64, f64) {
        let dx = self.h[i];
        let dy = y[i + 1] - y[i];
        let c_i = c[i];
        let c_ip1 = c[i + 1];
        let b = dy / dx - dx * fma(2.0, c_i, c_ip1) / 3.0;
        let d = (c_ip1 - c_i) / (3.0 * dx);
        (b, d)
    }

    /// [`Self::interval`] for the intervals `range`, written to `b[i]`,
    /// `d[i]`.
    pub(crate) fn intervals(
        &self,
        range: std::ops::Range<usize>,
        y: &[f64],
        c: &[f64],
        b: &mut [f64],
        d: &mut [f64],
    ) {
        let (lo, hi) = (range.start, range.end);
        let h = &self.h[lo..hi];
        let (y0, y1) = (&y[lo..hi], &y[lo + 1..hi + 1]);
        let (c0, c1) = (&c[lo..hi], &c[lo + 1..hi + 1]);
        for (k, (b, d)) in b[lo..hi].iter_mut().zip(&mut d[lo..hi]).enumerate() {
            let dx = h[k];
            let dy = y1[k] - y0[k];
            *b = dy / dx - dx * fma(2.0, c0[k], c1[k]) / 3.0;
            *d = (c1[k] - c0[k]) / (3.0 * dx);
        }
    }

    /// Index `i` with `x[i] <= t < x[i + 1]` (`len() - 2` at the last knot),
    /// as `gsl_interp_accel_find` returns it, searching up from `hint`.
    #[inline]
    pub(crate) fn locate_from(&self, t: f64, hint: usize) -> usize {
        let last = self.x.len() - 2;
        let mut i = hint.min(last);
        if t < self.x[i] {
            i = 0;
        }
        while i < last && self.x[i + 1] <= t {
            i += 1;
        }
        i
    }

    /// Whether GSL evaluates the spline at `t` without a domain error.
    #[inline]
    pub(crate) fn contains(&self, t: f64) -> bool {
        t >= self.x[0] && t <= self.x[self.x.len() - 1]
    }
}

/// Evaluate the cubic `a + delx (b + delx (c + delx d))` as GSL does.
#[inline(always)]
pub(crate) fn cubic(a: f64, b: f64, c: f64, d: f64, delx: f64) -> f64 {
    fma(delx, fma(delx, fma(delx, d, c), b), a)
}

/// Evaluation points of a spline on [`FixedKnots`]: the interval of each
/// point and its offset from the left knot.
#[derive(Clone, Debug)]
pub(crate) struct FixedPoints {
    pub(crate) index: Vec<usize>,
    pub(crate) delx: Vec<f64>,
}

impl FixedPoints {
    pub(crate) fn new(knots: &FixedKnots, t: &[f64]) -> Self {
        let mut hint = 0;
        let mut index = Vec::with_capacity(t.len());
        let mut delx = Vec::with_capacity(t.len());
        for &t in t {
            assert!(knots.contains(t), "evaluation point outside the knots");
            hint = knots.locate_from(t, hint);
            index.push(hint);
            delx.push(t - knots.x[hint]);
        }
        Self { index, delx }
    }
}

/// `gsl_interp_polynomial` through four points evaluated at `t`: Newton's
/// divided differences (`gsl_poly_dd_init`) and `gsl_poly_dd_eval`.
#[inline]
pub(crate) fn polynomial4(x: &[f64; 4], y: &[f64; 4], t: f64) -> f64 {
    let mut dd = *y;
    for j in (1..4).rev() {
        dd[j] = (y[j] - y[j - 1]) / (x[j] - x[j - 1]);
    }
    for i in 2..4 {
        for j in (i..4).rev() {
            dd[j] = (dd[j] - dd[j - 1]) / (x[j] - x[j - i]);
        }
    }
    let mut value = dd[3];
    for i in (0..3).rev() {
        value = fma(t - x[i], value, dd[i]);
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interpolate::CubicSpline;

    #[test]
    fn matches_general_spline() {
        let x: Vec<f64> = (0..40).map(|i| (i as f64 * 0.37).powf(1.3) - 3.0).collect();
        let y: Vec<f64> = x.iter().map(|x| (1.7 * x).sin() * x.exp().sqrt()).collect();
        let reference = CubicSpline::new(&x, &y);
        let knots = FixedKnots::new(&x);
        let mut c = vec![0.0; x.len()];
        let mut z = vec![0.0; x.len()];
        knots.second_derivatives(&y, &mut c, &mut z);
        let t: Vec<f64> = (0..500)
            .map(|i| x[0] + (x[x.len() - 1] - x[0]) * i as f64 / 499.0)
            .collect();
        let points = FixedPoints::new(&knots, &t);
        for (k, &t) in t.iter().enumerate() {
            let i = points.index[k];
            let (b, d) = knots.interval(i, &y, &c);
            let value = cubic(y[i], b, c[i], d, points.delx[k]);
            assert_eq!(value.to_bits(), reference.eval(t).to_bits(), "t = {t}");
        }
    }

    #[test]
    fn polynomial_is_exact_for_cubics() {
        let x = [0.1, 0.4, 1.3, 2.0];
        let f = |t: f64| 1.0 - 2.0 * t + 0.5 * t * t - 0.25 * t * t * t;
        let y = x.map(f);
        for t in [0.1, 0.5, 1.0, 1.9, 2.0] {
            assert!((polynomial4(&x, &y, t) - f(t)).abs() < 1e-14);
        }
    }
}
