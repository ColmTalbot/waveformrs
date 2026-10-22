// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Inlined, vectorisable elementary functions for the grid evaluation paths.
//!
//! The per-frequency evaluation of the models matches LALSimulation, which
//! calls the C library. When many waveforms are evaluated on a prepared
//! [`crate::grid::FrequencyGrid`], the remaining cost is dominated by
//! `sin_cos` of the phase; [`sin_cos`] here is branch-free so the compiler
//! can inline it and vectorise loops over frequencies.
//!
//! Accuracy: the argument is reduced modulo pi/2 with a three-part
//! Cody-Waite reduction whose first step is exact for `|x| < 2^30 pi/2`
//! (about 1.7e9); the kernels are the fdlibm `__kernel_sin`/`__kernel_cos`
//! polynomials (errors below one ulp on `[-pi/4, pi/4]`). The results agree
//! with the C library to a few ulp. Larger arguments fall back to the
//! standard library.
//!
//! Targets without a fused multiply-add instruction (baseline x86-64, i.e.
//! without `-C target-feature=+fma`) would turn `f64::mul_add` into a call
//! of the C library's software `fma`, which is slow and stops the loops
//! from vectorising. There the reduction splits the quadrant number instead
//! so that every product is exact (see [`reduce_pio2`]), and the polynomial
//! steps use a separately rounded multiply and add. Rounding to an integer
//! uses the `1.5 * 2^52` shift off aarch64, which (unlike `round_ties_even`
//! without SSE4.1, and the conversion to `i64` without AVX-512) vectorises.

// The constants are quoted with fdlibm's published digits.
#![allow(clippy::excessive_precision)]

use std::f64::consts::FRAC_2_PI;

// pi/2 split into three parts: PIO2_1 and PIO2_2 have 33 significant bits
// (fdlibm's pio2_1, pio2_2), PIO2_3 is the remainder (pio2_3 + pio2_3t).
const PIO2_1: f64 = 1.570_796_326_734_125_614_17; // 0x3FF921FB54400000
const PIO2_2: f64 = 6.077_100_506_303_965_976_6e-11; // 0x3DD0B4611A600000
const PIO2_3: f64 = 2.022_266_248_795_950_631_54e-21; // 0x3BA3198A2E037073 (pio2_2t)

// fdlibm __kernel_sin coefficients
const S1: f64 = -1.666_666_666_666_663_243_48e-01;
const S2: f64 = 8.333_333_333_322_489_461_09e-03;
const S3: f64 = -1.984_126_982_985_794_931_22e-04;
const S4: f64 = 2.755_731_370_707_006_767_89e-06;
const S5: f64 = -2.505_076_025_340_686_341_95e-08;
const S6: f64 = 1.589_690_995_211_550_102_21e-10;
// fdlibm __kernel_cos coefficients
const C1: f64 = 4.166_666_666_666_660_190_37e-02;
const C2: f64 = -1.388_888_888_887_411_009_60e-03;
const C3: f64 = 2.480_158_728_947_672_941_57e-05;
const C4: f64 = -2.755_731_435_139_066_330_35e-07;
const C5: f64 = 2.087_572_321_298_174_827_19e-09;
const C6: f64 = -1.135_964_755_778_819_482_94e-11;

/// Largest argument handled by the inlined reduction.
const REDUCTION_LIMIT: f64 = 1.0e9;

/// Whether `f64::mul_add` is a single instruction on this target. The
/// functions below take this as their `FUSED` parameter (the tests check
/// both variants on every target).
const HARDWARE_FMA: bool = crate::spintaylor::fp::FUSED_MULTIPLY_ADD;

/// `a * b + c`, fused or separately rounded (only used where both are
/// accurate enough).
#[inline(always)]
fn mul_add<const FUSED: bool>(a: f64, b: f64, c: f64) -> f64 {
    if FUSED {
        a.mul_add(b, c)
    } else {
        a * b + c
    }
}

/// `1.5 * 2^52`: for `|y| < 2^51`, `y + SHIFT` rounds `y` to the nearest
/// integer (ties to even) and holds that integer in its low mantissa bits.
const SHIFT: f64 = 6_755_399_441_055_744.0;

/// `y` rounded to the nearest integer, ties to even, as an `f64` and as the
/// low bits of an integer (two's complement; the bits above the 51st are not
/// meaningful), for `|y| < 2^51`.
#[inline(always)]
fn round_to_integer(y: f64) -> (f64, u64) {
    if cfg!(target_arch = "aarch64") {
        // frintn and fcvtzs, both vectorised
        let k = y.round_ties_even();
        (k, k as i64 as u64)
    } else {
        round_to_integer_shift(y)
    }
}

/// [`round_to_integer`] with the `1.5 * 2^52` shift.
#[inline(always)]
fn round_to_integer_shift(y: f64) -> (f64, u64) {
    let shifted = y + SHIFT;
    (shifted - SHIFT, shifted.to_bits())
}

/// `x - k pi/2` for the integer `k = round(x 2/pi)`, `|x| < 1e9`.
///
/// The first step, `x - k PIO2_1`, is exact. With fused multiply-adds each
/// step rounds once. Without them, `k` is split as `k_hi + k_lo` with `k_hi`
/// a multiple of 2^10 (at most 20 significant bits for `|k| < 2^30`) and
/// `|k_lo| <= 2^9`, so that the products with the 33-bit `PIO2_1` and
/// `PIO2_2` are exact. `x - k_hi PIO2_1` is a multiple of `ulp(x)` below
/// 2^10 in magnitude and `x - k PIO2_1` one below 1, so both differences are
/// exact as well, and only the second step rounds twice instead of once.
#[inline(always)]
fn reduce_pio2<const FUSED: bool>(x: f64, k: f64) -> f64 {
    if FUSED {
        let r = (-k).mul_add(PIO2_1, x);
        let r = (-k).mul_add(PIO2_2, r);
        (-k).mul_add(PIO2_3, r)
    } else {
        const TWO_10: f64 = 1024.0;
        let k_hi = round_to_integer_shift(k * (1.0 / TWO_10)).0 * TWO_10;
        let k_lo = k - k_hi;
        let r = (x - k_hi * PIO2_1) - k_lo * PIO2_1;
        let r = (r - k_hi * PIO2_2) - k_lo * PIO2_2;
        r - k * PIO2_3
    }
}

/// `(sin(x), cos(x))`, inlined; branch-free for `|x| < 1e9`.
#[inline(always)]
pub(crate) fn sin_cos(x: f64) -> (f64, f64) {
    if x.abs() >= REDUCTION_LIMIT || x.is_nan() {
        return x.sin_cos();
    }
    sin_cos_reduced(x)
}

/// [`sin_cos`] without the range check: requires `|x| < 1e9`.
#[inline(always)]
fn sin_cos_reduced(x: f64) -> (f64, f64) {
    sin_cos_reduced_with::<HARDWARE_FMA>(x)
}

#[inline(always)]
fn sin_cos_reduced_with<const FUSED: bool>(x: f64) -> (f64, f64) {
    let (k, q) = round_to_integer(x * FRAC_2_PI);
    let r = reduce_pio2::<FUSED>(x, k);
    let (s, c) = kernel_sin_cos(r);
    // quadrant k mod 4: (s, c) -> (s, c), (c, -s), (-s, -c), (-c, s),
    // written with selects so that loops vectorise
    let odd = q & 1 != 0;
    // flip the sign bit of both results in quadrants 2 and 3
    let sign = (q & 2) << 62;
    let sin = if odd { c } else { s };
    let cos = if odd { -s } else { c };
    (
        f64::from_bits(sin.to_bits() ^ sign),
        f64::from_bits(cos.to_bits() ^ sign),
    )
}

/// Replace each `(amplitude, phase)` stored as `Complex::new(amplitude,
/// phase)` by `amplitude * exp(-i phase)`, in a loop the compiler can
/// vectorise.
pub(crate) fn polar_to_complex_conj(values: &mut [num_complex::Complex<f64>]) {
    let limit_ok = values.iter().all(|v| v.im.abs() < REDUCTION_LIMIT);
    if limit_ok {
        for v in values.iter_mut() {
            let (s, c) = sin_cos_reduced(v.im);
            let a = v.re;
            *v = num_complex::Complex::new(a * c, -a * s);
        }
    } else {
        for v in values.iter_mut() {
            let (s, c) = sin_cos(v.im);
            let a = v.re;
            *v = num_complex::Complex::new(a * c, -a * s);
        }
    }
}

/// sin and cos on `[-pi/4, pi/4]` (fdlibm kernels without the tail term).
#[inline(always)]
fn kernel_sin_cos(r: f64) -> (f64, f64) {
    let z = r * r;
    let w = z * z;
    // sin: r + r^3 (S1 + z (S2 + z (S3 + z (S4 + z (S5 + z S6)))))
    let sp = S2 + z * (S3 + z * S4) + z * w * (S5 + z * S6);
    let sin = r + r * z * (S1 + z * sp);
    // cos: 1 - z/2 + z^2 (C1 + z C2 + ...), with fdlibm's correction for
    // the rounding of 1 - z/2
    let cp = z * (C1 + z * C2 + w * C3) + w * w * (C4 + z * C5 + w * C6);
    let hz = 0.5 * z;
    let one_minus_hz = 1.0 - hz;
    let cos = one_minus_hz + (((1.0 - one_minus_hz) - hz) + z * cp);
    (sin, cos)
}

/// Replace each value by its cosine, in a loop the compiler can vectorise.
pub(crate) fn cos_in_place(values: &mut [f64]) {
    if values.iter().all(|x| x.abs() < REDUCTION_LIMIT) {
        for x in values.iter_mut() {
            *x = sin_cos_reduced(*x).1;
        }
    } else {
        for x in values.iter_mut() {
            *x = sin_cos(*x).1;
        }
    }
}

/// Largest `|x|` handled by the inlined exponential (`2^n` stays normal).
const EXP_LIMIT: f64 = 700.0;

/// Replace each value by its exponential, in a loop the compiler can
/// vectorise: `x = n ln(2) + r` with `|r| <= ln(2) / 2`, `exp(r)` by its
/// Taylor polynomial of degree 13 (truncation below 1e-17) and the factor
/// `2^n` from the exponent bits. Within an ulp or two of the C library.
pub(crate) fn exp_in_place(values: &mut [f64]) {
    if values.iter().all(|x| x.abs() < EXP_LIMIT) {
        for x in values.iter_mut() {
            *x = exp_reduced(*x);
        }
    } else {
        for x in values.iter_mut() {
            *x = x.exp();
        }
    }
}

/// `exp(x)` for `|x| < EXP_LIMIT`.
#[inline(always)]
fn exp_reduced(x: f64) -> f64 {
    exp_reduced_with::<HARDWARE_FMA>(x)
}

#[inline(always)]
fn exp_reduced_with<const FUSED: bool>(x: f64) -> f64 {
    // ln(2) in two parts (fdlibm's ln2_hi, ln2_lo)
    const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-1;
    const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;
    // 1 / k! for k = 13 down to 2
    const INV_FACTORIAL: [f64; 12] = [
        1.605_904_383_682_161_3e-10,
        2.087_675_698_786_81e-9,
        2.505_210_838_544_172e-8,
        2.755_731_922_398_589e-7,
        2.755_731_922_398_589_3e-6,
        2.480_158_730_158_730_2e-5,
        1.984_126_984_126_984e-4,
        1.388_888_888_888_889e-3,
        8.333_333_333_333_333e-3,
        4.166_666_666_666_666_4e-2,
        1.666_666_666_666_666_6e-1,
        0.5,
    ];
    let (n, bits) = round_to_integer(x * std::f64::consts::LOG2_E);
    // n LN2_HI is exact (|n| < 2^11, LN2_HI has 32 significant bits) and so
    // is x - n LN2_HI: the reduction is accurate without fused operations
    let r = mul_add::<FUSED>(-n, LN2_LO, mul_add::<FUSED>(-n, LN2_HI, x));
    let mut p = INV_FACTORIAL[0];
    for &c in &INV_FACTORIAL[1..] {
        p = mul_add::<FUSED>(p, r, c);
    }
    // 1 + r + r^2 p
    let p = mul_add::<FUSED>(r * r, p, r) + 1.0;
    // 2^n, from the low bits of n
    let scale = f64::from_bits(bits.wrapping_add(1023) << 52);
    p * scale
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The variants with and without fused multiply-adds (only one of
    /// which is used on a given target) against the C library.
    #[test]
    fn fused_and_unfused_variants_agree_with_the_c_library() {
        let mut worst_sin_cos = [0.0f64; 2];
        let mut worst_exp = [0.0f64; 2];
        let mut differ = 0usize;
        let mut n = 0usize;
        // up to the reduction limit, with denser sampling at small |x|
        let mut x = 1e-3;
        while x < REDUCTION_LIMIT {
            for x in [x, -x, x * 1.000_000_7, x + 0.5] {
                let reference = x.sin_cos();
                let fused = sin_cos_reduced_with::<true>(x);
                let unfused = sin_cos_reduced_with::<false>(x);
                for (worst, (s, c)) in worst_sin_cos.iter_mut().zip([fused, unfused]) {
                    *worst = worst
                        .max((s - reference.0).abs())
                        .max((c - reference.1).abs());
                }
                differ += usize::from(fused != unfused);
                n += 1;
            }
            x *= 1.000_37;
        }
        let mut x = -EXP_LIMIT;
        while x < EXP_LIMIT {
            let reference = x.exp();
            for (worst, e) in worst_exp
                .iter_mut()
                .zip([exp_reduced_with::<true>(x), exp_reduced_with::<false>(x)])
            {
                *worst = worst.max(((e - reference) / reference).abs());
            }
            x += 0.000_731;
        }
        println!(
            "sin_cos: worst abs error {:.2e} (fused) {:.2e} (unfused), the variants differ in \
             {differ} of {n}; exp: worst rel error {:.2e} (fused) {:.2e} (unfused)",
            worst_sin_cos[0], worst_sin_cos[1], worst_exp[0], worst_exp[1]
        );
        for worst in worst_sin_cos {
            assert!(worst < 4.0 * f64::EPSILON, "{worst:e}");
        }
        for worst in worst_exp {
            assert!(worst < 3.0 * f64::EPSILON, "{worst:e}");
        }
        // rounding with the shift is the same as round_ties_even
        for y in [
            0.5,
            1.5,
            2.5,
            -0.5,
            -2.5,
            3.49,
            -7.51,
            6.4e8 + 0.5,
            -1e12 - 0.5,
        ] {
            let (k, bits) = round_to_integer_shift(y);
            assert_eq!(k, y.round_ties_even(), "{y}");
            assert_eq!(bits & 0xfff, (k as i64 as u64) & 0xfff, "{y}");
        }
    }

    #[test]
    fn exp_and_cos_agree_with_the_c_library() {
        let mut worst: f64 = 0.0;
        let x: Vec<f64> = (0..200_001).map(|i| -350.0 + 0.0035 * i as f64).collect();
        let mut values = x.clone();
        exp_in_place(&mut values);
        for (&x, &e) in x.iter().zip(&values) {
            let e0 = x.exp();
            worst = worst.max(((e - e0) / e0).abs());
        }
        assert!(worst < 3.0 * f64::EPSILON, "exp: {worst:e}");
        let mut values = x.clone();
        cos_in_place(&mut values);
        for (&x, &c) in x.iter().zip(&values) {
            assert!((c - x.cos()).abs() < 4.0 * f64::EPSILON);
        }
        // out of range and NaN go to the C library
        let mut values = vec![0.5, 800.0, f64::NAN];
        exp_in_place(&mut values);
        assert_eq!(values[..2], [0.5f64.exp(), 800.0f64.exp()]);
        assert!(values[2].is_nan());
    }

    fn ulps(a: f64, b: f64) -> f64 {
        ((a - b) / f64::EPSILON / b.abs().max(f64::MIN_POSITIVE)).abs()
    }

    #[test]
    fn agrees_with_the_c_library() {
        let mut worst_abs: f64 = 0.0;
        let mut worst_ulp: f64 = 0.0;
        let mut x = -3.0e6;
        while x < 3.0e6 {
            let (s, c) = sin_cos(x);
            let (s0, c0) = x.sin_cos();
            worst_abs = worst_abs.max((s - s0).abs()).max((c - c0).abs());
            if s0.abs() > 1e-3 {
                worst_ulp = worst_ulp.max(ulps(s, s0));
            }
            if c0.abs() > 1e-3 {
                worst_ulp = worst_ulp.max(ulps(c, c0));
            }
            x += 0.987_654_321;
        }
        for &x in &[
            0.0,
            1e-300,
            -1e-10,
            0.785,
            std::f64::consts::FRAC_PI_2,
            3.0,
            1e8 + 0.5,
        ] {
            let (s, c) = sin_cos(x);
            let (s0, c0) = x.sin_cos();
            worst_abs = worst_abs.max((s - s0).abs()).max((c - c0).abs());
        }
        assert!(worst_abs < 4.0 * f64::EPSILON, "{worst_abs:e}");
        assert!(worst_ulp < 3.0, "{worst_ulp}");
        let (s, c) = sin_cos(f64::NAN);
        assert!(s.is_nan() && c.is_nan());
        assert_eq!(sin_cos(2e9), 2e9_f64.sin_cos());
    }
}
