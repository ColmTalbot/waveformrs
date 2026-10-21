// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Values with their first and second derivatives with respect to the
//! frequency (second-order forward-mode differentiation).
//!
//! IMRPhenomXHM fixes some coefficients from derivatives of the ansatzes at
//! the region boundaries that LALSimulation takes by finite differences
//! (steps of 1e-9 for the amplitudes, 1e-7 for the second derivative of the
//! (3,2) ringdown phase). Those amplify the last-bit rounding of the
//! ansatzes by up to ~1e14. Evaluating the same ansatzes on [`Jet`]s gives
//! the derivatives the finite differences approximate, to rounding.

use std::ops::{Add, Div, Mul, Neg, Sub};

use num_complex::Complex;

/// `v + d1 e + d2 e^2 / 2`: a value with its first and second derivatives.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Jet {
    pub(crate) v: f64,
    pub(crate) d1: f64,
    pub(crate) d2: f64,
}

impl Jet {
    pub(crate) const fn new(v: f64, d1: f64, d2: f64) -> Self {
        Self { v, d1, d2 }
    }

    /// The independent variable at `x`.
    pub(crate) const fn variable(x: f64) -> Self {
        Self::new(x, 1., 0.)
    }

    pub(crate) const fn constant(c: f64) -> Self {
        Self::new(c, 0., 0.)
    }

    /// `x^p` of the independent variable `x`.
    pub(crate) fn power(x: f64, p: f64) -> Self {
        Self::new(x.powf(p), p * x.powf(p - 1.), p * (p - 1.) * x.powf(p - 2.))
    }

    /// `g(self)` for `g` with value `g0` and derivatives `g1`, `g2`.
    fn chain(self, g0: f64, g1: f64, g2: f64) -> Self {
        Self::new(g0, g1 * self.d1, g1 * self.d2 + g2 * self.d1 * self.d1)
    }

    pub(crate) fn exp(self) -> Self {
        let e = self.v.exp();
        self.chain(e, e, e)
    }

    pub(crate) fn ln(self) -> Self {
        let r = 1. / self.v;
        self.chain(self.v.ln(), r, -r * r)
    }

    pub(crate) fn recip(self) -> Self {
        let r = 1. / self.v;
        self.chain(r, -r * r, 2. * r * r * r)
    }

    pub(crate) fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        self.chain(s, 0.5 / s, -0.25 / (s * self.v))
    }

    /// `sqrt(self^2 + other^2)`.
    pub(crate) fn hypot(self, other: Self) -> Self {
        (self * self + other * other).sqrt()
    }
}

impl Add for Jet {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.v + o.v, self.d1 + o.d1, self.d2 + o.d2)
    }
}

impl Add<f64> for Jet {
    type Output = Self;
    fn add(self, c: f64) -> Self {
        Self::new(self.v + c, self.d1, self.d2)
    }
}

impl Sub for Jet {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.v - o.v, self.d1 - o.d1, self.d2 - o.d2)
    }
}

impl Sub<f64> for Jet {
    type Output = Self;
    fn sub(self, c: f64) -> Self {
        Self::new(self.v - c, self.d1, self.d2)
    }
}

impl Neg for Jet {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.v, -self.d1, -self.d2)
    }
}

impl Mul for Jet {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.v * o.v,
            self.d1 * o.v + self.v * o.d1,
            self.d2 * o.v + 2. * self.d1 * o.d1 + self.v * o.d2,
        )
    }
}

impl Mul<f64> for Jet {
    type Output = Self;
    fn mul(self, c: f64) -> Self {
        Self::new(self.v * c, self.d1 * c, self.d2 * c)
    }
}

impl Mul<Jet> for f64 {
    type Output = Jet;
    fn mul(self, j: Jet) -> Jet {
        j * self
    }
}

impl Div for Jet {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let q = self.v / o.v;
        let q1 = (self.d1 - q * o.d1) / o.v;
        let q2 = (self.d2 - 2. * q1 * o.d1 - q * o.d2) / o.v;
        Self::new(q, q1, q2)
    }
}

impl Div<f64> for Jet {
    type Output = Self;
    fn div(self, c: f64) -> Self {
        self * (1. / c)
    }
}

impl Div<Jet> for f64 {
    type Output = Jet;
    fn div(self, j: Jet) -> Jet {
        Jet::constant(self) / j
    }
}

/// `c A exp(i phi)` with its first and second derivatives, for a complex
/// constant `c` and real `A`, `phi`.
pub(crate) fn polar(c: Complex<f64>, amplitude: Jet, phase: Jet) -> [Complex<f64>; 3] {
    let i = Complex::i();
    let rotation = c * Complex::from_polar(1., phase.v);
    let (a, p) = (amplitude, phase);
    [
        rotation * a.v,
        rotation * (a.d1 + i * (a.v * p.d1)),
        rotation * (a.d2 - a.v * p.d1 * p.d1 + i * (2. * a.d1 * p.d1 + a.v * p.d2)),
    ]
}

/// The first and second derivatives of `arg h` from `h` and its
/// derivatives.
pub(crate) fn arg_derivatives(h: [Complex<f64>; 3]) -> (f64, f64) {
    let log1 = h[1] / h[0];
    (log1.im, (h[2] / h[0] - log1 * log1).im)
}

/// The first derivative of `|h|` from `h` and its derivatives.
pub(crate) fn abs_derivative(h: [Complex<f64>; 3]) -> f64 {
    (h[0].conj() * h[1]).re / h[0].norm()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Analytic derivatives against a fourth-order central difference.
    #[test]
    fn jet_rules() {
        let g = |x: Jet| ((x * x + 1.).ln() / x).exp() * x.sqrt() - 3. / x.hypot(Jet::constant(2.));
        let g0 = |x: f64| ((x * x + 1.).ln() / x).exp() * x.sqrt() - 3. / x.hypot(2.);
        let x = 0.7;
        let h = 1e-3;
        let d1 = (g0(x - 2. * h) - 8. * g0(x - h) + 8. * g0(x + h) - g0(x + 2. * h)) / (12. * h);
        let d2 = (-g0(x - 2. * h) + 16. * g0(x - h) - 30. * g0(x) + 16. * g0(x + h)
            - g0(x + 2. * h))
            / (12. * h * h);
        let j = g(Jet::variable(x));
        assert!((j.v - g0(x)).abs() < 1e-15);
        assert!((j.d1 - d1).abs() < 1e-10, "{} {}", j.d1, d1);
        assert!((j.d2 - d2).abs() < 1e-7, "{} {}", j.d2, d2);
        let p = Jet::power(x, -7. / 6.);
        let lnp = (-7. / 6. * Jet::variable(x).ln()).exp();
        assert!((p.d1 - lnp.d1).abs() < 1e-14 && (p.d2 - lnp.d2).abs() < 1e-13);
    }

    #[test]
    fn polar_and_arg() {
        let c = Complex::new(0.3, -1.1);
        let x = Jet::variable(0.4);
        let a = (x * 2.).exp();
        let p = x * x * 3.;
        let h = polar(c, a, p);
        let (dphi, d2phi) = arg_derivatives(h);
        assert!((dphi - 6. * 0.4).abs() < 1e-14);
        assert!((d2phi - 6.).abs() < 1e-13);
        let dabs = abs_derivative(h);
        assert!((dabs - c.norm() * 2. * (0.8f64).exp()).abs() < 1e-13);
    }
}
