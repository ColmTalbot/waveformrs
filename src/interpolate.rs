//! Natural cubic spline interpolation.
//!
//! This reproduces `gsl_interp_cspline` (natural boundary conditions) as used
//! throughout LALSimulation, including the order of floating point operations
//! in the tridiagonal solve, so that interpolated values agree with LAL to
//! rounding precision.

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
        let size = x.len();
        assert_eq!(size, y.len(), "x and y must have the same length");
        assert!(size >= 3, "cubic spline needs at least three points");
        let max_index = size - 1;
        let sys_size = max_index - 1;
        let mut c = vec![0.0; size];
        let mut g = vec![0.0; sys_size];
        let mut diag = vec![0.0; sys_size];
        let mut offdiag = vec![0.0; sys_size];
        for i in 0..sys_size {
            let h_i = x[i + 1] - x[i];
            let h_ip1 = x[i + 2] - x[i + 1];
            let ydiff_i = y[i + 1] - y[i];
            let ydiff_ip1 = y[i + 2] - y[i + 1];
            let g_i = if h_i != 0.0 { 1.0 / h_i } else { 0.0 };
            let g_ip1 = if h_ip1 != 0.0 { 1.0 / h_ip1 } else { 0.0 };
            offdiag[i] = h_ip1;
            diag[i] = 2.0 * (h_ip1 + h_i);
            g[i] = 3.0 * (ydiff_ip1 * g_ip1 - ydiff_i * g_i);
        }
        if sys_size == 1 {
            c[1] = g[0] / diag[0];
        } else {
            solve_symmetric_tridiagonal(&diag, &offdiag, &g, &mut c[1..=sys_size]);
        }
        Self {
            x: x.to_vec(),
            y: y.to_vec(),
            c,
        }
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
        let b = dy / dx - dx * (c_ip1 + 2.0 * c_i) / 3.0;
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
        a + delx * (b + delx * (c + delx * d))
    }

    /// Evaluate the first derivative of the spline.
    pub fn derivative(&self, value: f64) -> f64 {
        let index = self.find(value);
        let (_, b, c, d) = self.coefficients(index);
        let delx = value - self.x[index];
        b + delx * (2.0 * c + 3.0 * d * delx)
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
            *out = a + delx * (b + delx * (c + delx * d));
        }
    }

    pub fn x(&self) -> &[f64] {
        &self.x
    }
}

/// Port of `gsl_linalg_solve_symm_tridiag` (LDL^T factorisation).
fn solve_symmetric_tridiagonal(diag: &[f64], offdiag: &[f64], b: &[f64], x: &mut [f64]) {
    let n = diag.len();
    let mut gamma = vec![0.0; n];
    let mut alpha = vec![0.0; n];
    let mut z = vec![0.0; n];
    alpha[0] = diag[0];
    gamma[0] = offdiag[0] / alpha[0];
    for i in 1..n - 1 {
        alpha[i] = diag[i] - offdiag[i - 1] * gamma[i - 1];
        gamma[i] = offdiag[i] / alpha[i];
    }
    if n > 1 {
        alpha[n - 1] = diag[n - 1] - offdiag[n - 2] * gamma[n - 2];
    }
    z[0] = b[0];
    for i in 1..n {
        z[i] = b[i] - gamma[i - 1] * z[i - 1];
    }
    for i in 0..n {
        z[i] /= alpha[i];
    }
    x[n - 1] = z[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = z[i] - gamma[i] * x[i + 1];
    }
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
