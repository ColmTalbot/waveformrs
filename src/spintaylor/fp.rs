// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Fused multiply-adds.
//!
//! `a * b + c` with a single rounding ([`f64::mul_add`]) is at least as
//! accurate as a separately rounded product and sum, and on targets with an
//! FMA instruction it is also as fast or faster. [`fma`], [`fms`] and
//! [`fnms`] therefore fuse wherever the target has the instruction
//! ([`FUSED_MULTIPLY_ADD`]: aarch64, and x86-64 built with the `fma` target
//! feature, e.g. `-C target-cpu=native` on Haswell or later). On other
//! targets `f64::mul_add` would be a call of the C library's software `fma`,
//! which is slow and stops loops from vectorising, so the product and the
//! sum are rounded separately there.
//!
//! The expressions of the crate were translated from C, and the helpers mark
//! where a product is fused with a sum; which product of `a * b + c * d` is
//! fused follows clang's choice, which is as good as any other.

/// Whether [`fma`] and friends evaluate `a * b + c` with a single rounding:
/// where `f64::mul_add` is a single instruction.
pub const FUSED_MULTIPLY_ADD: bool = cfg!(any(target_arch = "aarch64", target_feature = "fma"));

/// `a * b + c`, fused where the target has an FMA instruction.
#[inline(always)]
pub(crate) fn fma(a: f64, b: f64, c: f64) -> f64 {
    if FUSED_MULTIPLY_ADD {
        a.mul_add(b, c)
    } else {
        a * b + c
    }
}

/// `a * b - c`, fused where the target has an FMA instruction.
#[inline(always)]
pub(crate) fn fms(a: f64, b: f64, c: f64) -> f64 {
    if FUSED_MULTIPLY_ADD {
        a.mul_add(b, -c)
    } else {
        a * b - c
    }
}

/// `c - a * b`, fused where the target has an FMA instruction.
#[inline(always)]
pub(crate) fn fnms(a: f64, b: f64, c: f64) -> f64 {
    if FUSED_MULTIPLY_ADD {
        (-a).mul_add(b, c)
    } else {
        c - a * b
    }
}

/// Dot product `ax*bx + ay*by + az*bz`.
#[inline(always)]
pub(crate) fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    fma(a[2], b[2], fma(a[0], b[0], a[1] * b[1]))
}

/// Cross product (`XLALSimInspiralVectorCrossProduct`).
#[inline(always)]
pub(crate) fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        fms(a[1], b[2], a[2] * b[1]),
        fms(a[2], b[0], a[0] * b[2]),
        fms(a[0], b[1], a[1] * b[0]),
    ]
}
