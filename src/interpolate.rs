// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Natural cubic spline interpolation.
//!
//! This reproduces `gsl_interp_cspline` (natural boundary conditions) as used
//! throughout LALSimulation, including the order of floating point operations
//! in the tridiagonal solve, so that interpolated values agree with LAL to
//! rounding precision (with fused multiply-adds where the hardware has them,
//! [`crate::spintaylor::fp`]).

use crate::spintaylor::fp::{fma, fms, fnms};

/// A natural cubic spline through `(x_i, y_i)`.
#[derive(Clone, Debug)]
pub struct CubicSpline {
    x: Vec<f64>,
    y: Vec<f64>,
    /// Half the second derivative at each node (GSL's `c` array).
    c: Vec<f64>,
}

impl CubicSpline {
    /// Construct the spline. `x` must be strictly increasing and contain at
    /// least three points.
    pub fn new(x: &[f64], y: &[f64]) -> Self {
        Self::from_vecs(x.to_vec(), y.to_vec())
    }

    /// [`CubicSpline::new`] taking ownership of the knots and values.
    pub fn from_vecs(x: Vec<f64>, y: Vec<f64>) -> Self {
        let [c] = natural_c(&x, [&y]);
        Self { x, y, c }
    }

    /// The splines through `(x, y1)` and `(x, y2)`, the same as
    /// [`CubicSpline::new`] for each but sharing the work that depends only
    /// on the knots.
    pub fn pair(x: Vec<f64>, y1: Vec<f64>, y2: Vec<f64>) -> (Self, Self) {
        let [c1, c2] = natural_c(&x, [&y1, &y2]);
        let first = Self {
            x: x.clone(),
            y: y1,
            c: c1,
        };
        (first, Self { x, y: y2, c: c2 })
    }

    /// Index `i` such that `x[i] <= value < x[i + 1]` (clamped to the table).
    #[inline]
    fn find(&self, value: f64) -> usize {
        let mut lo = 0;
        let mut hi = self.x.len() - 1;
        while hi > lo + 1 {
            let mid = (hi + lo) / 2;
            if self.x[mid] > value {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        lo
    }

    #[inline]
    fn coefficients(&self, index: usize) -> (f64, f64, f64, f64) {
        let dx = self.x[index + 1] - self.x[index];
        let dy = self.y[index + 1] - self.y[index];
        let c_i = self.c[index];
        let c_ip1 = self.c[index + 1];
        let b = dy / dx - dx * fma(2.0, c_i, c_ip1) / 3.0;
        let d = (c_ip1 - c_i) / (3.0 * dx);
        (self.y[index], b, c_i, d)
    }

    /// Evaluate the spline. Values outside the knots return NaN (as GSL
    /// signals a domain error).
    pub fn eval(&self, value: f64) -> f64 {
        if !(value >= self.x[0] && value <= self.x[self.x.len() - 1]) {
            return f64::NAN;
        }
        self.eval_unchecked(value)
    }

    /// Evaluate the spline, extrapolating with the boundary polynomials.
    #[inline]
    pub fn eval_unchecked(&self, value: f64) -> f64 {
        let index = self.find(value);
        let (a, b, c, d) = self.coefficients(index);
        let delx = value - self.x[index];
        fma(delx, fma(delx, fma(delx, d, c), b), a)
    }

    /// Evaluate the first derivative of the spline.
    pub fn derivative(&self, value: f64) -> f64 {
        let index = self.find(value);
        let (_, b, c, d) = self.coefficients(index);
        let delx = value - self.x[index];
        fma(delx, fma(2.0, c, 3.0 * d * delx), b)
    }

    /// Evaluate the spline at a sorted sequence of points, walking the knots
    /// rather than bisecting for every point.
    pub fn eval_sorted(&self, values: &[f64], output: &mut [f64]) {
        let mut index = 0;
        let n = self.x.len();
        for (value, out) in values.iter().zip(output.iter_mut()) {
            while index + 2 < n && self.x[index + 1] <= *value {
                index += 1;
            }
            let (a, b, c, d) = self.coefficients(index);
            let delx = value - self.x[index];
            *out = fma(delx, fma(delx, fma(delx, d, c), b), a);
        }
    }

    pub fn x(&self) -> &[f64] {
        &self.x
    }

    /// The polynomial coefficients `(a, b, c, d)` of every interval, as
    /// [`CubicSpline::eval_in`] computes them, for evaluating the spline
    /// many times (`a + delx (b + delx (c + delx d))` with `delx = value -
    /// x[i]` gives the values of `eval_in` bit for bit).
    pub(crate) fn coefficient_table(&self) -> Vec<[f64; 4]> {
        (0..self.x.len() - 1)
            .map(|index| {
                let (a, b, c, d) = self.coefficients(index);
                [a, b, c, d]
            })
            .collect()
    }

    /// Whether `value` lies within the knots, i.e. GSL would evaluate the
    /// spline without a domain error.
    #[inline]
    pub fn contains(&self, value: f64) -> bool {
        value >= self.x[0] && value <= self.x[self.x.len() - 1]
    }

    /// Index `i` of the interval with `x[i] <= value < x[i + 1]` (clamped to
    /// the table), searching from the interval `hint` as
    /// `gsl_interp_accel_find` does. The result does not depend on `hint`.
    #[inline]
    pub fn locate(&self, value: f64, hint: usize) -> usize {
        let last = self.x.len() - 2;
        let mut index = hint.min(last);
        while index < last && self.x[index + 1] <= value {
            index += 1;
        }
        while index > 0 && self.x[index] > value {
            index -= 1;
        }
        index
    }

    /// Evaluate the spline polynomial of interval `index` at `value`.
    #[inline]
    pub fn eval_in(&self, index: usize, value: f64) -> f64 {
        let (a, b, c, d) = self.coefficients(index);
        let delx = value - self.x[index];
        fma(delx, fma(delx, fma(delx, d, c), b), a)
    }

    /// Evaluate the derivative of the spline polynomial of interval `index`
    /// at `value`.
    #[inline]
    pub fn derivative_in(&self, index: usize, value: f64) -> f64 {
        let (_, b, c, d) = self.coefficients(index);
        let delx = value - self.x[index];
        fma(delx, fma(2.0, c, 3.0 * d * delx), b)
    }
}

/// The `c` arrays (half the second derivatives) of the natural cubic
/// splines through `(x, y_k)`: `cspline_init` of GSL, which sets up the
/// tridiagonal system for the interior knots and solves it with
/// `gsl_linalg_solve_symm_tridiag` (an LDL^T factorisation).
///
/// The operations on each value are GSL's, in GSL's order, but the set-up,
/// the factorisation and the forward substitution are done in one pass, the
/// reciprocal of each knot spacing is computed once rather than twice, and
/// the factorisation (which depends only on `x`) is shared by the splines.
/// Only one array of the factorisation is stored, which matters for splines
/// with millions of knots.
fn natural_c<const K: usize>(x: &[f64], ys: [&[f64]; K]) -> [Vec<f64>; K] {
    let size = x.len();
    for y in ys {
        assert_eq!(size, y.len(), "x and y must have the same length");
    }
    assert!(size >= 3, "cubic spline needs at least three points");
    // the size of the system, for the interior knots 1 ..= n
    let n = size - 2;
    let mut cs: [Vec<f64>; K] = std::array::from_fn(|_| vec![0.0; size]);
    // gamma of the LDL^T factorisation, `offdiag[i] / alpha[i]`
    let mut gamma = vec![0.0; n.saturating_sub(1)];
    let reciprocal = |h: f64| if h != 0.0 { 1.0 / h } else { 0.0 };
    let mut h_i = x[1] - x[0];
    let mut g_i = reciprocal(h_i);
    // the forward substitution before the division by alpha
    let mut z = [0.0; K];
    for i in 0..n {
        let h_ip1 = x[i + 2] - x[i + 1];
        let g_ip1 = reciprocal(h_ip1);
        // diag[i] = 2 (h_ip1 + h_i) and offdiag[i - 1] = h_i
        let diag = 2.0 * (h_ip1 + h_i);
        let alpha = if i == 0 {
            diag
        } else {
            fnms(h_i, gamma[i - 1], diag)
        };
        if i + 1 < n {
            gamma[i] = h_ip1 / alpha;
        }
        for k in 0..K {
            let y = ys[k];
            let ydiff_i = y[i + 1] - y[i];
            let ydiff_ip1 = y[i + 2] - y[i + 1];
            let b = 3.0 * fms(ydiff_ip1, g_ip1, ydiff_i * g_i);
            z[k] = if i == 0 {
                b
            } else {
                fnms(gamma[i - 1], z[k], b)
            };
            cs[k][i + 1] = z[k] / alpha;
        }
        h_i = h_ip1;
        g_i = g_ip1;
    }
    // back substitution
    for c in cs.iter_mut() {
        for i in (0..n - 1).rev() {
            c[i + 1] = fnms(gamma[i], c[i + 2], c[i + 1]);
        }
    }
    cs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproduces_cubic_interior() {
        // A natural spline through a straight line is exact.
        let x: Vec<f64> = (0..10).map(|i| i as f64 * 0.5).collect();
        let y: Vec<f64> = x.iter().map(|x| 2.0 * x - 1.0).collect();
        let spline = CubicSpline::new(&x, &y);
        for v in [0.0, 0.3, 1.7, 4.49] {
            assert!((spline.eval(v) - (2.0 * v - 1.0)).abs() < 1e-14);
            assert!((spline.derivative(v) - 2.0).abs() < 1e-13);
        }
        assert!(spline.eval(-0.1).is_nan());
    }

    /// `cspline_init` with `gsl_linalg_solve_symm_tridiag` as separate
    /// steps, as GSL does them.
    fn gsl_c(x: &[f64], y: &[f64]) -> Vec<f64> {
        let sys_size = x.len() - 2;
        let (mut g, mut diag, mut offdiag) = (
            vec![0.0; sys_size],
            vec![0.0; sys_size],
            vec![0.0; sys_size],
        );
        for i in 0..sys_size {
            let h_i = x[i + 1] - x[i];
            let h_ip1 = x[i + 2] - x[i + 1];
            let g_i = if h_i != 0.0 { 1.0 / h_i } else { 0.0 };
            let g_ip1 = if h_ip1 != 0.0 { 1.0 / h_ip1 } else { 0.0 };
            offdiag[i] = h_ip1;
            diag[i] = 2.0 * (h_ip1 + h_i);
            g[i] = 3.0 * fms(y[i + 2] - y[i + 1], g_ip1, (y[i + 1] - y[i]) * g_i);
        }
        let mut c = vec![0.0; x.len()];
        let n = sys_size;
        if n == 1 {
            c[1] = g[0] / diag[0];
            return c;
        }
        let (mut gamma, mut alpha, mut z) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        alpha[0] = diag[0];
        gamma[0] = offdiag[0] / alpha[0];
        for i in 1..n - 1 {
            alpha[i] = fnms(offdiag[i - 1], gamma[i - 1], diag[i]);
            gamma[i] = offdiag[i] / alpha[i];
        }
        alpha[n - 1] = fnms(offdiag[n - 2], gamma[n - 2], diag[n - 1]);
        z[0] = g[0];
        for i in 1..n {
            z[i] = fnms(gamma[i - 1], z[i - 1], g[i]);
        }
        for i in 0..n {
            z[i] /= alpha[i];
        }
        c[n] = z[n - 1];
        for i in (0..n - 1).rev() {
            c[i + 1] = fnms(gamma[i], c[i + 2], z[i]);
        }
        c
    }

    /// The fused construction gives GSL's coefficients to rounding, also
    /// for two splines sharing their knots.
    #[test]
    fn constructors_match_gsl() {
        for n in [3, 4, 5, 50] {
            let x: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).powf(1.3)).collect();
            let y1: Vec<f64> = x.iter().map(|x| x.sin()).collect();
            let y2: Vec<f64> = x.iter().map(|x| (0.3 * x).exp()).collect();
            let (a, b) = CubicSpline::pair(x.clone(), y1.clone(), y2.clone());
            assert_eq!(a.c, gsl_c(&x, &y1));
            assert_eq!(b.c, gsl_c(&x, &y2));
            assert_eq!(CubicSpline::new(&x, &y1).c, gsl_c(&x, &y1));
            assert_eq!(
                CubicSpline::from_vecs(x.clone(), y2.clone()).c,
                gsl_c(&x, &y2)
            );
        }
    }

    #[test]
    fn sorted_matches_scalar() {
        let x: Vec<f64> = (0..20).map(|i| (i as f64 * 0.3).powi(2)).collect();
        let y: Vec<f64> = x.iter().map(|x| x.sin()).collect();
        let spline = CubicSpline::new(&x, &y);
        let values: Vec<f64> = (0..100).map(|i| i as f64 * 0.32).collect();
        let mut out = vec![0.0; values.len()];
        spline.eval_sorted(&values, &mut out);
        for (v, o) in values.iter().zip(&out) {
            assert_eq!(spline.eval_unchecked(*v), *o);
        }
    }
}
