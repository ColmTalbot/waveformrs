// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! The C library's complex exponential.
//!
//! LALSimulation builds the modes as `amp * cexp(I * phi)`. For a purely
//! imaginary argument the platform `cexp` does not always return the
//! correctly rounded `(cos(phi), sin(phi))` (on macOS it differs in the last
//! bit for about 0.2% of arguments). The (3,2) mode of IMRPhenomXHM
//! differentiates such phases numerically with a step of 1e-7, so `cexp` is
//! called through FFI to reproduce LAL exactly.

use num_complex::Complex;

/// C `double complex`, which has the same layout and calling convention as
/// a struct of two doubles on the supported platforms.
#[repr(C)]
#[derive(Clone, Copy)]
struct CComplex {
    re: f64,
    im: f64,
}

extern "C" {
    fn cexp(z: CComplex) -> CComplex;
}

/// `cexp(I * phi)` as evaluated by C: the argument is `(0 * phi, phi)`.
#[inline]
pub(crate) fn cexp_i(phi: f64) -> Complex<f64> {
    // SAFETY: `cexp` is a pure function of its argument.
    let z = unsafe {
        cexp(CComplex {
            re: 0.0 * phi,
            im: phi,
        })
    };
    Complex::new(z.re, z.im)
}

#[cfg(test)]
mod tests {
    use super::cexp_i;

    #[test]
    fn unit_modulus() {
        for &phi in &[0.0, 0.3, -2.5, 1234.5678, -98765.4321] {
            let z = cexp_i(phi);
            assert!((z.re - phi.cos()).abs() <= 2.0 * f64::EPSILON);
            assert!((z.im - phi.sin()).abs() <= 2.0 * f64::EPSILON);
        }
    }
}
