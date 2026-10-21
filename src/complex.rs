// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Elementary functions of complex numbers.

use num_complex::Complex;

/// `exp(i phi) = cos(phi) + i sin(phi)`.
#[inline]
pub(crate) fn cexp_i(phi: f64) -> Complex<f64> {
    let (sin, cos) = phi.sin_cos();
    Complex::new(cos, sin)
}

/// The principal square root, with the real part non-negative and the
/// imaginary part of the sign of `z.im`. Computed as in C99's `csqrt`:
/// `t = sqrt((|x| + |z|) / 2)` adds two non-negative numbers, so it has no
/// cancellation, and the other part is `|y| / (2 t)`.
pub(crate) fn csqrt(z: Complex<f64>) -> Complex<f64> {
    let (x, y) = (z.re, z.im);
    if x == 0.0 && y == 0.0 {
        return Complex::new(0.0, y);
    }
    if y.is_infinite() {
        return Complex::new(f64::INFINITY, y);
    }
    if x.is_nan() || y.is_nan() {
        return Complex::new(f64::NAN, f64::NAN);
    }
    if x.is_infinite() {
        return if x > 0.0 {
            Complex::new(x, 0.0f64.copysign(y))
        } else {
            Complex::new(0.0, f64::INFINITY.copysign(y))
        };
    }
    // halving first keeps |x| + |z| finite near the largest doubles
    let t = (0.5 * x.abs() + 0.5 * x.hypot(y)).sqrt();
    if x >= 0.0 {
        Complex::new(t, y / (2.0 * t))
    } else {
        Complex::new(y.abs() / (2.0 * t), t.copysign(y))
    }
}

/// The principal cube root `|z|^(1/3) exp(i arg(z) / 3)`, the value of C's
/// `cpow(z, 1/3)`, with the modulus from `cbrt` rather than `exp(log|z| /
/// 3)`.
pub(crate) fn principal_cbrt(z: Complex<f64>) -> Complex<f64> {
    if z.im == 0.0 && z.re >= 0.0 {
        return Complex::new(z.re.cbrt(), z.im);
    }
    let (r, theta) = z.to_polar();
    r.cbrt() * cexp_i(theta / 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: Complex<f64>, b: Complex<f64>, ulps: f64) -> bool {
        (a - b).norm() <= ulps * f64::EPSILON * b.norm().max(f64::MIN_POSITIVE)
    }

    #[test]
    fn exponential() {
        for phi in [0.0, 0.3, -2.5, 1234.5678, -98765.4321] {
            let z = cexp_i(phi);
            assert_eq!((z.re, z.im), (phi.cos(), phi.sin()));
        }
    }

    #[test]
    fn square_root() {
        for (re, im) in [
            (4.0, 0.0),
            (-4.0, 0.0),
            (-4.0, -0.0),
            (3.0, 4.0),
            (-3.0, 4.0),
            (-3.0, -4.0),
            (1e-300, 1e-300),
            (1e300, -1e300),
            (-1e-8, 1.0),
            (1.0, 1e-20),
        ] {
            let z = Complex::new(re, im);
            let w = csqrt(z);
            assert!(w.re >= 0.0, "{z}");
            assert_eq!(w.im.is_sign_negative(), im.is_sign_negative(), "{z}");
            assert!(close(w * w, z, 4.0), "{z}: {w}");
        }
        assert_eq!(csqrt(Complex::new(-4.0, 0.0)), Complex::new(0.0, 2.0));
        assert_eq!(
            csqrt(Complex::new(0.0, -0.0)).im.to_bits(),
            (-0.0f64).to_bits()
        );
        assert_eq!(csqrt(Complex::new(-f64::INFINITY, 1.0)).im, f64::INFINITY);
    }

    #[test]
    fn cube_root() {
        for (re, im) in [
            (8.0, 0.0),
            (-8.0, 0.0),
            (3.0, 4.0),
            (-3.0, -4.0),
            (1e-200, -3e-200),
            (2.5e7, 1e-9),
        ] {
            let z = Complex::new(re, im);
            let w = principal_cbrt(z);
            assert!(close(w * w * w, z, 8.0), "{z}: {w}");
            // the principal branch: arg(w) = arg(z) / 3
            assert!(
                (w.arg() - z.arg() / 3.0).abs() <= 4.0 * f64::EPSILON,
                "{z}: {w}"
            );
        }
        assert_eq!(
            principal_cbrt(Complex::new(27.0, 0.0)),
            Complex::new(3.0, 0.0)
        );
    }
}
