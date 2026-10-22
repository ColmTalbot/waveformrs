// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Quasi-normal mode fits for the (2,2) mode, `LALSimIMRPhenomX_qnm.c`
//! (`QNMfits == 1`, the default build of LALSimulation).
//!
//! Both fits are in units of the remnant mass; divide by the final mass in
//! units of the total mass to get frequencies in units of the total mass.
//!
//! The denominators of these rational fits suffer from strong cancellation
//! for large spins, which amplifies rounding differences by up to ~1e4. The
//! sums are therefore accumulated with `fp::fma`, fused where the target
//! LALSimulation build contracts them (arm64), which matches LAL to the last
//! bit.

/// Accumulate `terms[0] + sum_i coefficients[i] * powers[i]` left to right
/// with (where LAL fuses them) fused multiply-adds.
#[inline]
fn fma_sum<const N: usize>(constant: f64, coefficients: [f64; N], powers: [f64; N]) -> f64 {
    coefficients
        .iter()
        .zip(powers)
        .fold(constant, |acc, (&c, x)| {
            crate::spintaylor::fp::fma(c, x, acc)
        })
}

/// `evaluate_QNMfit_fring22`: ringdown frequency of the (2,2) mode for
/// dimensionless final spin `af`, |af| <= 1.
pub(crate) fn fring22(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    let x7 = x4 * x3;
    let numerator = fma_sum(
        0.05947169566573468,
        [
            -0.14989771215394762,
            0.09535606290986028,
            0.02260924869042963,
            -0.02501704155363241,
            -0.005852438240997211,
            0.0027489038393367993,
            0.0005821983163192694,
        ],
        [af, x2, x3, x4, x5, x6, x7],
    );
    let denominator = fma_sum(
        1.0,
        [
            -2.8570126619966296,
            2.373335413978394,
            -0.6036964688511505,
            0.0873798215084077,
        ],
        [af, x2, x4, x6],
    );
    numerator / denominator
}

/// `evaluate_QNMfit_fdamp22`: damping frequency of the (2,2) mode for
/// dimensionless final spin `af`, |af| <= 1.
pub(crate) fn fdamp22(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    let numerator = fma_sum(
        0.014158792290965177,
        [
            -0.036989395871554566,
            0.026822526296575368,
            0.0008490933750566702,
            -0.004843996907020524,
            -0.00014745235759327472,
            0.0001504546201236794,
        ],
        [af, x2, x3, x4, x5, x6],
    );
    let denominator = fma_sum(
        1.0,
        [
            -2.5900842798681376,
            1.8952576220623967,
            -0.31416610693042507,
            0.009002719412204133,
        ],
        [af, x2, x4, x6],
    );
    numerator / denominator
}
