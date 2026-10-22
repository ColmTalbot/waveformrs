// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Coefficients of the intermediate amplitude polynomial of the IMRPhenomXHM
//! 122019 release, `IMRPhenomXHM_Intermediate_Amp_delta0` ...
//! `IMRPhenomXHM_Intermediate_Amp_delta5` in
//! `LALSimIMRPhenomXHM_intermediate.c`, for the reconstructions the release
//! can reach (101, 102, 1032, 104, 1042, 105).
//!
//! The intermediate amplitude is the inverse of a polynomial in `f` fixed by
//! values `v_i` of the inverse amplitude at frequencies `f_i` and derivatives
//! `d1`, `d4` at the boundaries. Generated from the LAL sources with LAL's
//! fused multiply-adds (see [`fma`]); long sums are written as
//! accumulations, which does not change the operations.

use crate::imrphenomx::internals::fma;

/// The inputs of the reconstruction (`pAmp->v1` ... `pAmp->d4`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct PolyInputs {
    pub(crate) d1: f64,
    pub(crate) d4: f64,
    pub(crate) v1: f64,
    pub(crate) v2: f64,
    pub(crate) v3: f64,
    pub(crate) v4: f64,
    pub(crate) f1: f64,
    pub(crate) f2: f64,
    pub(crate) f3: f64,
    pub(crate) f4: f64,
}

/// `delta0` ... `delta5` for reconstruction `version`
/// (`Update_Intermediate_Amplitude_Coefficients`).
///
/// # Panics
///
/// If `version` is not a reconstruction of the 122019 release.
pub(crate) fn coefficients(p: &PolyInputs, version: i32) -> [f64; 6] {
    type Delta = fn(&PolyInputs) -> f64;
    let deltas: [Delta; 6] = match version {
        101 => [
            delta0_101, delta1_101, delta2_101, delta3_101, delta4_101, delta5_101,
        ],
        102 => [
            delta0_102, delta1_102, delta2_102, delta3_102, delta4_102, delta5_102,
        ],
        1032 => [
            delta0_1032,
            delta1_1032,
            delta2_1032,
            delta3_1032,
            delta4_1032,
            delta5_1032,
        ],
        104 => [
            delta0_104, delta1_104, delta2_104, delta3_104, delta4_104, delta5_104,
        ],
        1042 => [
            delta0_1042,
            delta1_1042,
            delta2_1042,
            delta3_1042,
            delta4_1042,
            delta5_1042,
        ],
        105 => [
            delta0_105, delta1_105, delta2_105, delta3_105, delta4_105, delta5_105,
        ],
        _ => panic!("IMRPhenomXHM: intermediate amplitude version {version} is not valid"),
    };
    deltas.map(|delta| delta(p))
}

/// `IMRPhenomXHM_Intermediate_Amp_delta0`, reconstruction 101.
pub(crate) fn delta0_101(p: &PolyInputs) -> f64 {
    let f1mf4 = p.f1 - p.f4;
    fma(-p.f4, p.v1, p.f1 * p.v4) / f1mf4
}

/// `IMRPhenomXHM_Intermediate_Amp_delta0`, reconstruction 102.
pub(crate) fn delta0_102(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f42 = p.f4 * p.f4;
    let f1mf4 = p.f1 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    fma(
        -(2.0 * p.f1 * p.f4),
        p.v4,
        fma(f12, p.v4, fma(-(p.d4 * p.f1 * f1mf4), p.f4, f42 * p.v1)),
    ) / f1mf42
}

/// `IMRPhenomXHM_Intermediate_Amp_delta0`, reconstruction 1032.
pub(crate) fn delta0_1032(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f1mf4 = p.f1 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf43 = f1mf42 * f1mf4;
    fma(
        -(3.0 * f12 * p.f4),
        p.v4,
        fma(
            f13,
            p.v4,
            fma(
                -f43,
                p.v1,
                fma(
                    3.0 * p.f1 * f42,
                    p.v1,
                    fma(
                        p.d4 * f12 * p.f4,
                        -p.f1 + p.f4,
                        p.d1 * p.f1 * (-p.f1 + p.f4) * f42,
                    ),
                ),
            ),
        ),
    ) / f1mf43
}

/// `IMRPhenomXHM_Intermediate_Amp_delta0`, reconstruction 104.
pub(crate) fn delta0_104(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f42 = p.f4 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf4 = p.f1 - p.f4;
    let f2mf4 = p.f2 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    fma(
        f12 * f1mf22
            * p.f2
            * fma(
                4.0,
                f42,
                fma(-(3.0 * p.f2), p.f4, fma(p.f1, p.f2, -(2.0 * p.f1 * p.f4))),
            ),
        p.v4,
        fma(
            f42,
            fma(
                p.f2 * f2mf42
                    * fma(
                        -p.f2,
                        p.f4,
                        fma(2.0 * p.f1, p.f4, fma(-4.0, f12, 3.0 * p.f1 * p.f2)),
                    ),
                p.v1,
                f12 * f1mf43 * p.v2,
            ),
            fma(
                -(p.d4 * f12 * f1mf22 * f1mf4 * p.f2 * f2mf4),
                p.f4,
                p.d1 * p.f1 * f1mf2 * f1mf4 * p.f2 * f2mf42 * f42,
            ),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta0`, reconstruction 1042.
pub(crate) fn delta0_1042(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f1mf3 = p.f1 - p.f3;
    let f3mf4 = p.f3 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf32 = f1mf3 * f1mf3;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf42 * f1mf4;
    let mut acc0 = p.d1 * p.f1 * f1mf3 * f1mf4 * p.f3 * f3mf42 * f42;
    acc0 = fma(-(p.d4 * f12 * f1mf32 * f1mf4 * p.f3 * f3mf4), p.f4, acc0);
    acc0 = fma(-(4.0 * f12 * f33 * f42), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f34 * f42, p.v1, acc0);
    acc0 = fma(8.0 * f12 * f32 * f43, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f33 * f43), p.v1, acc0);
    acc0 = fma(-(f34 * f43), p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * p.f3 * f44), p.v1, acc0);
    acc0 = fma(-(p.f1 * f32 * f44), p.v1, acc0);
    acc0 = fma(2.0 * f33 * f44, p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * p.f3 * f45, p.v1, acc0);
    acc0 = fma(-(f32 * f45), p.v1, acc0);
    acc0 = fma(f15 * f42, p.v3, acc0);
    acc0 = fma(-(3.0 * f14 * f43), p.v3, acc0);
    acc0 = fma(3.0 * f13 * f44, p.v3, acc0);
    acc0 = fma(-(f12 * f45), p.v3, acc0);
    acc0 = fma(f15 * f32, p.v4, acc0);
    acc0 = fma(-(2.0 * f14 * f33), p.v4, acc0);
    acc0 = fma(f13 * f34, p.v4, acc0);
    acc0 = fma(-(2.0 * f15 * p.f3 * p.f4), p.v4, acc0);
    acc0 = fma(f14 * f32 * p.f4, p.v4, acc0);
    acc0 = fma(4.0 * f13 * f33 * p.f4, p.v4, acc0);
    acc0 = fma(-(3.0 * f12 * f34 * p.f4), p.v4, acc0);
    acc0 = fma(4.0 * f14 * p.f3 * f42, p.v4, acc0);
    acc0 = fma(-(8.0 * f13 * f32 * f42), p.v4, acc0);
    acc0 = fma(4.0 * f12 * f33 * f42, p.v4, acc0);
    acc0 / (f1mf32 * f1mf43 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta0`, reconstruction 105.
pub(crate) fn delta0_105(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f16 = f15 * p.f1;
    let f17 = f16 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f25 = f24 * p.f2;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f35 = f34 * p.f3;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f46 = f45 * p.f4;
    let f47 = f46 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf3 = p.f1 - p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f2mf3 = p.f2 - p.f3;
    let f2mf4 = p.f2 - p.f4;
    let f3mf4 = p.f3 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f1mf32 = f1mf3 * f1mf3;
    let f1mf42 = f1mf4 * f1mf4;
    let f2mf42 = f2mf4 * f2mf4;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf42 * f1mf4;
    let mut acc0 =
        -(p.d1 * p.f1 * f1mf2 * f1mf3 * f1mf4 * p.f2 * f2mf3 * f2mf42 * p.f3 * f3mf42 * f42);
    acc0 = fma(
        -(p.d4 * f12 * f1mf22 * f1mf32 * f1mf4 * p.f2 * f2mf3 * f2mf4 * p.f3 * f3mf4),
        p.f4,
        acc0,
    );
    acc0 = fma(5.0 * f13 * f24 * f33 * f42, p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * f25 * f33 * f42), p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * f23 * f34 * f42), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f25 * f34 * f42, p.v1, acc0);
    acc0 = fma(4.0 * f12 * f23 * f35 * f42, p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f24 * f35 * f42), p.v1, acc0);
    acc0 = fma(-(10.0 * f13 * f24 * f32 * f43), p.v1, acc0);
    acc0 = fma(8.0 * f12 * f25 * f32 * f43, p.v1, acc0);
    acc0 = fma(5.0 * f12 * f24 * f33 * f43, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f25 * f33 * f43), p.v1, acc0);
    acc0 = fma(10.0 * f13 * f22 * f34 * f43, p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * f23 * f34 * f43), p.v1, acc0);
    acc0 = fma(-(f25 * f34 * f43), p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * f22 * f35 * f43), p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f23 * f35 * f43, p.v1, acc0);
    acc0 = fma(f24 * f35 * f43, p.v1, acc0);
    acc0 = fma(5.0 * f13 * f24 * p.f3 * f44, p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * f25 * p.f3 * f44), p.v1, acc0);
    acc0 = fma(15.0 * f13 * f23 * f32 * f44, p.v1, acc0);
    acc0 = fma(-(10.0 * f12 * f24 * f32 * f44), p.v1, acc0);
    acc0 = fma(-(p.f1 * f25 * f32 * f44), p.v1, acc0);
    acc0 = fma(-(15.0 * f13 * f22 * f33 * f44), p.v1, acc0);
    acc0 = fma(5.0 * p.f1 * f24 * f33 * f44, p.v1, acc0);
    acc0 = fma(2.0 * f25 * f33 * f44, p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * p.f2 * f34 * f44), p.v1, acc0);
    acc0 = fma(10.0 * f12 * f22 * f34 * f44, p.v1, acc0);
    acc0 = fma(-(5.0 * p.f1 * f23 * f34 * f44), p.v1, acc0);
    acc0 = fma(4.0 * f12 * p.f2 * f35 * f44, p.v1, acc0);
    acc0 = fma(p.f1 * f22 * f35 * f44, p.v1, acc0);
    acc0 = fma(-(2.0 * f23 * f35 * f44), p.v1, acc0);
    acc0 = fma(-(10.0 * f13 * f23 * p.f3 * f45), p.v1, acc0);
    acc0 = fma(5.0 * f12 * f24 * p.f3 * f45, p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * f25 * p.f3 * f45, p.v1, acc0);
    acc0 = fma(-(f12 * f23 * f32 * f45), p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * f24 * f32 * f45, p.v1, acc0);
    acc0 = fma(-(f25 * f32 * f45), p.v1, acc0);
    acc0 = fma(10.0 * f13 * p.f2 * f33 * f45, p.v1, acc0);
    acc0 = fma(f12 * f22 * f33 * f45, p.v1, acc0);
    acc0 = fma(-(3.0 * f24 * f33 * f45), p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * p.f2 * f34 * f45), p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * f22 * f34 * f45), p.v1, acc0);
    acc0 = fma(3.0 * f23 * f34 * f45, p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * p.f2 * f35 * f45), p.v1, acc0);
    acc0 = fma(f22 * f35 * f45, p.v1, acc0);
    acc0 = fma(5.0 * f13 * f22 * p.f3 * f46, p.v1, acc0);
    acc0 = fma(2.0 * f12 * f23 * p.f3 * f46, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f24 * p.f3 * f46), p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * p.f2 * f32 * f46), p.v1, acc0);
    acc0 = fma(-(p.f1 * f23 * f32 * f46), p.v1, acc0);
    acc0 = fma(2.0 * f24 * f32 * f46, p.v1, acc0);
    acc0 = fma(-(2.0 * f12 * p.f2 * f33 * f46), p.v1, acc0);
    acc0 = fma(p.f1 * f22 * f33 * f46, p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * p.f2 * f34 * f46, p.v1, acc0);
    acc0 = fma(-(2.0 * f22 * f34 * f46), p.v1, acc0);
    acc0 = fma(-(3.0 * f12 * f22 * p.f3 * f47), p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * f23 * p.f3 * f47, p.v1, acc0);
    acc0 = fma(3.0 * f12 * p.f2 * f32 * f47, p.v1, acc0);
    acc0 = fma(-(f23 * f32 * f47), p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * p.f2 * f33 * f47), p.v1, acc0);
    acc0 = fma(f22 * f33 * f47, p.v1, acc0);
    acc0 = fma(-(f17 * f33 * f42), p.v2, acc0);
    acc0 = fma(2.0 * f16 * f34 * f42, p.v2, acc0);
    acc0 = fma(-(f15 * f35 * f42), p.v2, acc0);
    acc0 = fma(2.0 * f17 * f32 * f43, p.v2, acc0);
    acc0 = fma(-(f16 * f33 * f43), p.v2, acc0);
    acc0 = fma(-(4.0 * f15 * f34 * f43), p.v2, acc0);
    acc0 = fma(3.0 * f14 * f35 * f43, p.v2, acc0);
    acc0 = fma(-(f17 * p.f3 * f44), p.v2, acc0);
    acc0 = fma(-(4.0 * f16 * f32 * f44), p.v2, acc0);
    acc0 = fma(8.0 * f15 * f33 * f44, p.v2, acc0);
    acc0 = fma(-(3.0 * f13 * f35 * f44), p.v2, acc0);
    acc0 = fma(3.0 * f16 * p.f3 * f45, p.v2, acc0);
    acc0 = fma(-(8.0 * f14 * f33 * f45), p.v2, acc0);
    acc0 = fma(4.0 * f13 * f34 * f45, p.v2, acc0);
    acc0 = fma(f12 * f35 * f45, p.v2, acc0);
    acc0 = fma(-(3.0 * f15 * p.f3 * f46), p.v2, acc0);
    acc0 = fma(4.0 * f14 * f32 * f46, p.v2, acc0);
    acc0 = fma(f13 * f33 * f46, p.v2, acc0);
    acc0 = fma(-(2.0 * f12 * f34 * f46), p.v2, acc0);
    acc0 = fma(f14 * p.f3 * f47, p.v2, acc0);
    acc0 = fma(-(2.0 * f13 * f32 * f47), p.v2, acc0);
    acc0 = fma(f12 * f33 * f47, p.v2, acc0);
    acc0 = fma(f17 * f23 * f42, p.v3, acc0);
    acc0 = fma(-(2.0 * f16 * f24 * f42), p.v3, acc0);
    acc0 = fma(f15 * f25 * f42, p.v3, acc0);
    acc0 = fma(-(2.0 * f17 * f22 * f43), p.v3, acc0);
    acc0 = fma(f16 * f23 * f43, p.v3, acc0);
    acc0 = fma(4.0 * f15 * f24 * f43, p.v3, acc0);
    acc0 = fma(-(3.0 * f14 * f25 * f43), p.v3, acc0);
    acc0 = fma(f17 * p.f2 * f44, p.v3, acc0);
    acc0 = fma(4.0 * f16 * f22 * f44, p.v3, acc0);
    acc0 = fma(-(8.0 * f15 * f23 * f44), p.v3, acc0);
    acc0 = fma(3.0 * f13 * f25 * f44, p.v3, acc0);
    acc0 = fma(-(3.0 * f16 * p.f2 * f45), p.v3, acc0);
    acc0 = fma(8.0 * f14 * f23 * f45, p.v3, acc0);
    acc0 = fma(-(4.0 * f13 * f24 * f45), p.v3, acc0);
    acc0 = fma(-(f12 * f25 * f45), p.v3, acc0);
    acc0 = fma(3.0 * f15 * p.f2 * f46, p.v3, acc0);
    acc0 = fma(-(4.0 * f14 * f22 * f46), p.v3, acc0);
    acc0 = fma(-(f13 * f23 * f46), p.v3, acc0);
    acc0 = fma(2.0 * f12 * f24 * f46, p.v3, acc0);
    acc0 = fma(-(f14 * p.f2 * f47), p.v3, acc0);
    acc0 = fma(2.0 * f13 * f22 * f47, p.v3, acc0);
    acc0 = fma(-(f12 * f23 * f47), p.v3, acc0);
    acc0 = fma(
        f12 * f1mf22
            * f1mf32
            * p.f2
            * f2mf3
            * p.f3
            * fma(
                p.f4,
                fma(
                    -5.0,
                    f42,
                    fma(-3.0 * p.f2, p.f3, 4.0 * (p.f2 + p.f3) * p.f4),
                ),
                p.f1 * fma(3.0, f42, fma(p.f2, p.f3, -(2.0 * (p.f2 + p.f3) * p.f4))),
            ),
        p.v4,
        acc0,
    );
    acc0 / (f1mf22 * f1mf32 * f1mf43 * f2mf3 * f2mf42 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta1`, reconstruction 101.
pub(crate) fn delta1_101(p: &PolyInputs) -> f64 {
    let f1mf4 = p.f1 - p.f4;
    (p.v1 - p.v4) / f1mf4
}

/// `IMRPhenomXHM_Intermediate_Amp_delta1`, reconstruction 102.
pub(crate) fn delta1_102(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f42 = p.f4 * p.f4;
    let f1mf42 = (p.f1 - p.f4) * (p.f1 - p.f4);
    fma(p.d4, f12 - f42, 2.0 * p.f4 * (-p.v1 + p.v4)) / f1mf42
}

/// `IMRPhenomXHM_Intermediate_Amp_delta1`, reconstruction 1032.
pub(crate) fn delta1_1032(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f42 = p.f4 * p.f4;
    let f1mf4 = p.f1 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf43 = f1mf42 * f1mf4;
    fma(
        p.d4 * p.f1 * f1mf4,
        fma(2.0, p.f4, p.f1),
        -(p.f4
            * fma(
                p.d1,
                fma(-2.0, f12, p.f1 * p.f4) + f42,
                6.0 * p.f1 * (p.v1 - p.v4),
            )),
    ) / f1mf43
}

/// `IMRPhenomXHM_Intermediate_Amp_delta1`, reconstruction 104.
pub(crate) fn delta1_104(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf4 = p.f1 - p.f4;
    let f2mf4 = p.f2 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    fma(
        p.d4 * p.f1 * f1mf22 * f1mf4 * f2mf4,
        fma(2.0 * p.f2, p.f4, p.f1 * (p.f2 + p.f4)),
        p.f4 * fma(
            -(p.d1 * f1mf2 * f1mf4 * f2mf42),
            fma(2.0 * p.f1, p.f2, (p.f1 + p.f2) * p.f4),
            -(2.0
                * p.f1
                * fma(
                    p.f1,
                    fma(
                        4.0 * f23,
                        -p.v1 + p.v4,
                        fma(2.0 * f43, -p.v1 + p.v2, 6.0 * f22 * p.f4 * (p.v1 - p.v4)),
                    ),
                    fma(
                        2.0 * f13 * p.f4,
                        -p.v2 + p.v4,
                        fma(
                            4.0 * f23 * p.f4,
                            -p.v1 + p.v4,
                            fma(
                                f14,
                                p.v2 - p.v4,
                                fma(f44, p.v1 - p.v2, 3.0 * f24 * (p.v1 - p.v4)),
                            ),
                        ),
                    ),
                )),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta1`, reconstruction 1042.
pub(crate) fn delta1_1042(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f1mf3 = p.f1 - p.f3;
    let f3mf4 = p.f3 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf32 = f1mf3 * f1mf3;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf42 * f1mf4;
    let mut acc0 = -(2.0 * p.d4 * f14 * f33);
    acc0 = fma(p.d4 * f15, f32, acc0);
    acc0 = fma(p.d4 * f13, f34, acc0);
    acc0 = fma(p.d4 * f14 * f32, p.f4, acc0);
    acc0 = fma(-(2.0 * p.d1 * f13 * f33), p.f4, acc0);
    acc0 = fma(-(2.0 * p.d4 * f13 * f33), p.f4, acc0);
    acc0 = fma(2.0 * p.d1 * f12 * f34, p.f4, acc0);
    acc0 = fma(p.d4 * f12 * f34, p.f4, acc0);
    acc0 = fma(-(p.d4 * f15), f42, acc0);
    acc0 = fma(3.0 * p.d1 * f13 * f32, f42, acc0);
    acc0 = fma(p.d4 * f13 * f32, f42, acc0);
    acc0 = fma(-(2.0 * p.d1 * f12 * f33), f42, acc0);
    acc0 = fma(2.0 * p.d4 * f12 * f33, f42, acc0);
    acc0 = fma(-(p.d1 * p.f1 * f34), f42, acc0);
    acc0 = fma(-(2.0 * p.d4 * p.f1 * f34), f42, acc0);
    acc0 = fma(p.d4 * f14, f43, acc0);
    acc0 = fma(-(p.d1 * f12 * f32), f43, acc0);
    acc0 = fma(-(3.0 * p.d4 * f12 * f32), f43, acc0);
    acc0 = fma(2.0 * p.d1 * p.f1 * f33, f43, acc0);
    acc0 = fma(2.0 * p.d4 * p.f1 * f33, f43, acc0);
    acc0 = fma(-(p.d1 * f34), f43, acc0);
    acc0 = fma(-(p.d1 * f13), f44, acc0);
    acc0 = fma(-(p.d1 * p.f1 * f32), f44, acc0);
    acc0 = fma(2.0 * p.d1 * f33, f44, acc0);
    acc0 = fma(p.d1 * f12, f45, acc0);
    acc0 = fma(-(p.d1 * f32), f45, acc0);
    acc0 = fma(8.0 * f12 * f33 * p.f4, p.v1, acc0);
    acc0 = fma(-(6.0 * p.f1 * f34 * p.f4), p.v1, acc0);
    acc0 = fma(-(12.0 * f12 * f32 * f42), p.v1, acc0);
    acc0 = fma(8.0 * p.f1 * f33 * f42, p.v1, acc0);
    acc0 = fma(4.0 * f12 * f44, p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * f45), p.v1, acc0);
    acc0 = fma(-(2.0 * f15 * p.f4), p.v3, acc0);
    acc0 = fma(4.0 * f14 * f42, p.v3, acc0);
    acc0 = fma(-(4.0 * f12 * f44), p.v3, acc0);
    acc0 = fma(2.0 * p.f1 * f45, p.v3, acc0);
    acc0 = fma(2.0 * f15 * p.f4, p.v4, acc0);
    acc0 = fma(-(8.0 * f12 * f33 * p.f4), p.v4, acc0);
    acc0 = fma(6.0 * p.f1 * f34 * p.f4, p.v4, acc0);
    acc0 = fma(-(4.0 * f14 * f42), p.v4, acc0);
    acc0 = fma(12.0 * f12 * f32 * f42, p.v4, acc0);
    acc0 = fma(-(8.0 * p.f1 * f33 * f42), p.v4, acc0);
    acc0 / (f1mf32 * f1mf43 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta1`, reconstruction 105.
pub(crate) fn delta1_105(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f16 = f15 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f25 = f24 * p.f2;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f35 = f34 * p.f3;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf3 = p.f1 - p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f2mf3 = p.f2 - p.f3;
    let f2mf4 = p.f2 - p.f4;
    let f3mf4 = p.f3 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f1mf32 = f1mf3 * f1mf3;
    let f2mf42 = f2mf4 * f2mf4;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 = f13
        * p.f4
        * fma(
            5.0 * f24 * p.f4,
            -p.v3 + p.v4,
            fma(
                4.0 * f35,
                -p.v2 + p.v4,
                fma(
                    4.0 * f25,
                    p.v3 - p.v4,
                    fma(f45, -p.v2 + p.v3, 5.0 * f34 * p.f4 * (p.v2 - p.v4)),
                ),
            ),
        );
    acc0 = fma(
        f16,
        fma(
            2.0 * f23,
            -p.v3 + p.v4,
            fma(
                3.0 * f32 * p.f4,
                -p.v2 + p.v4,
                fma(
                    3.0 * f22 * p.f4,
                    p.v3 - p.v4,
                    fma(f43, p.v2 - p.v3, 2.0 * f33 * (p.v2 - p.v4)),
                ),
            ),
        ),
        acc0,
    );
    let mut acc1 = 2.0 * f35 * (p.v2 - p.v4);
    acc1 = fma(3.0 * f45, p.v2 - p.v3, acc1);
    acc1 = fma(5.0 * f34 * p.f4, p.v2 - p.v4, acc1);
    acc1 = fma(10.0 * f23 * f42, p.v3 - p.v4, acc1);
    acc1 = fma(10.0 * f33 * f42, -p.v2 + p.v4, acc1);
    acc1 = fma(2.0 * f25, -p.v3 + p.v4, acc1);
    acc1 = fma(5.0 * f24 * p.f4, -p.v3 + p.v4, acc1);
    acc0 = fma(f14, acc1, acc0);
    let mut acc2 = 2.0 * f33 * p.f4 * (p.v2 - p.v4);
    acc2 = fma(3.0 * f44, -p.v2 + p.v3, acc2);
    acc2 = fma(5.0 * f32 * f42, p.v2 - p.v4, acc2);
    acc2 = fma(4.0 * f24, p.v3 - p.v4, acc2);
    acc2 = fma(4.0 * f34, -p.v2 + p.v4, acc2);
    acc2 = fma(2.0 * f23 * p.f4, -p.v3 + p.v4, acc2);
    acc2 = fma(5.0 * f22 * f42, -p.v3 + p.v4, acc2);
    acc0 = fma(f15, acc2, acc0);
    acc0 = fma(
        -(5.0 * f12),
        fma(
            f22 * p.f4,
            fma(
                4.0 * f33 * p.f4,
                -p.v1 + p.v4,
                fma(f44, p.v1 - p.v3, 3.0 * f34 * (p.v1 - p.v4)),
            ),
            fma(
                f24,
                fma(
                    3.0 * f32 * p.f4,
                    -p.v1 + p.v4,
                    fma(f43, p.v1 - p.v3, 2.0 * f33 * (p.v1 - p.v4)),
                ),
                fma(
                    -(f32 * f3mf42 * f43),
                    p.v1 - p.v2,
                    2.0 * f23
                        * fma(
                            f34,
                            -p.v1 + p.v4,
                            fma(f44, -p.v1 + p.v3, 2.0 * f32 * f42 * (p.v1 - p.v4)),
                        ),
                ),
            ),
        ),
        acc0,
    );
    acc0 = fma(
        p.f1,
        fma(
            3.0 * f22 * p.f4,
            fma(
                5.0 * f34 * p.f4,
                -p.v1 + p.v4,
                fma(f45, p.v1 - p.v3, 4.0 * f35 * (p.v1 - p.v4)),
            ),
            fma(
                -(5.0 * f24 * p.f4),
                fma(
                    3.0 * f32 * p.f4,
                    -p.v1 + p.v4,
                    fma(f43, p.v1 - p.v3, 2.0 * f33 * (p.v1 - p.v4)),
                ),
                fma(
                    4.0 * f25,
                    fma(
                        3.0 * f32 * p.f4,
                        -p.v1 + p.v4,
                        fma(f43, p.v1 - p.v3, 2.0 * f33 * (p.v1 - p.v4)),
                    ),
                    fma(
                        -(f32 * f3mf42 * fma(4.0, p.f3, 3.0 * p.f4) * f43),
                        p.v1 - p.v2,
                        2.0 * f23
                            * fma(
                                4.0 * f35,
                                -p.v1 + p.v4,
                                fma(f45, -p.v1 + p.v3, 5.0 * f34 * p.f4 * (p.v1 - p.v4)),
                            ),
                    ),
                ),
            ),
        ),
        acc0,
    );
    acc0 = fma(
        -2.0,
        fma(
            f23 * p.f4,
            fma(
                5.0 * f34 * p.f4,
                -p.v1 + p.v4,
                fma(f45, p.v1 - p.v3, 4.0 * f35 * (p.v1 - p.v4)),
            ),
            fma(
                f25,
                fma(
                    4.0 * f33 * p.f4,
                    -p.v1 + p.v4,
                    fma(f44, p.v1 - p.v3, 3.0 * f34 * (p.v1 - p.v4)),
                ),
                fma(
                    -(f33 * f3mf42 * f44),
                    p.v1 - p.v2,
                    f24 * fma(
                        3.0 * f35,
                        -p.v1 + p.v4,
                        fma(2.0 * f45, -p.v1 + p.v3, 5.0 * f33 * f42 * (p.v1 - p.v4)),
                    ),
                ),
            ),
        ),
        acc0,
    );
    fma(
        p.d4 * p.f1 * f1mf22 * f1mf32 * f1mf4 * f2mf3 * f2mf4 * f3mf4,
        fma(
            p.f1 * (p.f2 + p.f3),
            p.f4,
            fma(p.f1 * p.f2, p.f3, 2.0 * p.f2 * p.f3 * p.f4),
        ),
        p.f4 * fma(
            p.d1 * f1mf2 * f1mf3 * f1mf4 * f2mf3 * f2mf42 * f3mf42,
            fma(
                p.f1 * (p.f2 + p.f3),
                p.f4,
                fma(2.0 * p.f1 * p.f2, p.f3, p.f2 * p.f3 * p.f4),
            ),
            p.f1 * acc0,
        ),
    ) / (f1mf22 * f1mf32 * f1mf43 * f2mf3 * f2mf42 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta2`, reconstruction 101.
pub(crate) fn delta2_101(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta2`, reconstruction 102.
pub(crate) fn delta2_102(p: &PolyInputs) -> f64 {
    let f1mf4 = p.f1 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    (fma(-p.d4, f1mf4, p.v1) - p.v4) / f1mf42
}

/// `IMRPhenomXHM_Intermediate_Amp_delta2`, reconstruction 1032.
pub(crate) fn delta2_1032(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f42 = p.f4 * p.f4;
    let f1mf4 = p.f1 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf43 = f1mf42 * f1mf4;
    fma(
        3.0 * (p.f1 + p.f4),
        p.v1 - p.v4,
        fma(
            -p.d1,
            fma(-2.0, f42, fma(p.f1, p.f4, f12)),
            p.d4 * (fma(-2.0, f12, p.f1 * p.f4) + f42),
        ),
    ) / f1mf43
}

/// `IMRPhenomXHM_Intermediate_Amp_delta2`, reconstruction 104.
pub(crate) fn delta2_104(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf4 = p.f1 - p.f4;
    let f2mf4 = p.f2 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 =
        p.d1 * f1mf2 * f1mf4 * f2mf42 * (fma(p.f1, p.f2, 2.0 * (p.f1 + p.f2) * p.f4) + f42);
    acc0 = fma(
        -(p.d4 * f1mf22 * f1mf4 * f2mf4),
        fma(2.0 * p.f1, p.f2 + p.f4, fma(p.f2, p.f4, f12)),
        acc0,
    );
    acc0 = fma(-(4.0 * f12 * f23), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f24, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f23 * p.f4), p.v1, acc0);
    acc0 = fma(3.0 * f24 * p.f4, p.v1, acc0);
    acc0 = fma(12.0 * f12 * p.f2 * f42, p.v1, acc0);
    acc0 = fma(-(4.0 * f23 * f42), p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * f43), p.v1, acc0);
    acc0 = fma(p.f1 * f44, p.v1, acc0);
    acc0 = fma(f45, p.v1, acc0);
    acc0 = fma(f15, p.v2, acc0);
    acc0 = fma(f14 * p.f4, p.v2, acc0);
    acc0 = fma(-(8.0 * f13 * f42), p.v2, acc0);
    acc0 = fma(8.0 * f12 * f43, p.v2, acc0);
    acc0 = fma(-(p.f1 * f44), p.v2, acc0);
    acc0 = fma(-f45, p.v2, acc0);
    acc0 = fma(
        -(f1mf22
            * fma(
                p.f1 * fma(3.0, p.f2, -(4.0 * p.f4)),
                fma(2.0, p.f4, p.f2),
                fma(
                    f12,
                    fma(2.0, p.f2, p.f4),
                    fma(p.f2 * fma(3.0, p.f2, -(4.0 * p.f4)), p.f4, f13),
                ),
            )),
        p.v4,
        acc0,
    );
    acc0 / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta2`, reconstruction 1042.
pub(crate) fn delta2_1042(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f1mf3 = p.f1 - p.f3;
    let f3mf4 = p.f3 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf32 = f1mf3 * f1mf3;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf42 * f1mf4;
    let mut acc0 =
        p.d1 * f1mf3 * f1mf4 * f3mf42 * (fma(p.f1, p.f3, 2.0 * (p.f1 + p.f3) * p.f4) + f42);
    acc0 = fma(
        -(p.d4 * f1mf32 * f1mf4 * f3mf4),
        fma(2.0 * p.f1, p.f3 + p.f4, fma(p.f3, p.f4, f12)),
        acc0,
    );
    acc0 = fma(-(4.0 * f12 * f33), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f34, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f33 * p.f4), p.v1, acc0);
    acc0 = fma(3.0 * f34 * p.f4, p.v1, acc0);
    acc0 = fma(12.0 * f12 * p.f3 * f42, p.v1, acc0);
    acc0 = fma(-(4.0 * f33 * f42), p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * f43), p.v1, acc0);
    acc0 = fma(p.f1 * f44, p.v1, acc0);
    acc0 = fma(f45, p.v1, acc0);
    acc0 = fma(f15, p.v3, acc0);
    acc0 = fma(f14 * p.f4, p.v3, acc0);
    acc0 = fma(-(8.0 * f13 * f42), p.v3, acc0);
    acc0 = fma(8.0 * f12 * f43, p.v3, acc0);
    acc0 = fma(-(p.f1 * f44), p.v3, acc0);
    acc0 = fma(-f45, p.v3, acc0);
    acc0 = fma(
        -(f1mf32
            * fma(
                p.f1 * fma(3.0, p.f3, -(4.0 * p.f4)),
                fma(2.0, p.f4, p.f3),
                fma(
                    f12,
                    fma(2.0, p.f3, p.f4),
                    fma(p.f3 * fma(3.0, p.f3, -(4.0 * p.f4)), p.f4, f13),
                ),
            )),
        p.v4,
        acc0,
    );
    acc0 / (f1mf32 * f1mf43 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta2`, reconstruction 105.
pub(crate) fn delta2_105(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f16 = f15 * p.f1;
    let f17 = f16 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f25 = f24 * p.f2;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f35 = f34 * p.f3;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f46 = f45 * p.f4;
    let f47 = f46 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf3 = p.f1 - p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f2mf3 = p.f2 - p.f3;
    let f2mf4 = p.f2 - p.f4;
    let f3mf4 = p.f3 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f1mf32 = f1mf3 * f1mf3;
    let f2mf42 = f2mf4 * f2mf4;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 = -(p.d1
        * f1mf2
        * f1mf3
        * f1mf4
        * f2mf3
        * f2mf42
        * f3mf42
        * fma(
            p.f1 + p.f2 + p.f3,
            f42,
            fma(
                p.f1 * p.f2,
                p.f3,
                2.0 * fma(p.f2, p.f3, p.f1 * (p.f2 + p.f3)) * p.f4,
            ),
        ));
    acc0 = fma(
        -(p.d4 * f1mf22 * f1mf32 * f1mf4 * f2mf3 * f2mf4 * f3mf4),
        fma(
            2.0 * p.f1,
            fma(p.f2, p.f3, (p.f2 + p.f3) * p.f4),
            fma(p.f2 * p.f3, p.f4, f12 * (p.f2 + p.f3 + p.f4)),
        ),
        acc0,
    );
    acc0 = fma(5.0 * f13 * f24 * f33, p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * f25 * f33), p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * f23 * f34), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f25 * f34, p.v1, acc0);
    acc0 = fma(4.0 * f12 * f23 * f35, p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f24 * f35), p.v1, acc0);
    acc0 = fma(5.0 * f12 * f24 * f33 * p.f4, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f25 * f33 * p.f4), p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * f23 * f34 * p.f4), p.v1, acc0);
    acc0 = fma(3.0 * f25 * f34 * p.f4, p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f23 * f35 * p.f4, p.v1, acc0);
    acc0 = fma(-(3.0 * f24 * f35 * p.f4), p.v1, acc0);
    acc0 = fma(-(15.0 * f13 * f24 * p.f3 * f42), p.v1, acc0);
    acc0 = fma(12.0 * f12 * f25 * p.f3 * f42, p.v1, acc0);
    acc0 = fma(5.0 * p.f1 * f24 * f33 * f42, p.v1, acc0);
    acc0 = fma(-(4.0 * f25 * f33 * f42), p.v1, acc0);
    acc0 = fma(15.0 * f13 * p.f2 * f34 * f42, p.v1, acc0);
    acc0 = fma(-(5.0 * p.f1 * f23 * f34 * f42), p.v1, acc0);
    acc0 = fma(-(12.0 * f12 * p.f2 * f35 * f42), p.v1, acc0);
    acc0 = fma(4.0 * f23 * f35 * f42, p.v1, acc0);
    acc0 = fma(10.0 * f13 * f24 * f43, p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * f25 * f43), p.v1, acc0);
    acc0 = fma(20.0 * f13 * f23 * p.f3 * f43, p.v1, acc0);
    acc0 = fma(-(15.0 * f12 * f24 * p.f3 * f43), p.v1, acc0);
    acc0 = fma(-(20.0 * f13 * p.f2 * f33 * f43), p.v1, acc0);
    acc0 = fma(5.0 * f24 * f33 * f43, p.v1, acc0);
    acc0 = fma(-(10.0 * f13 * f34 * f43), p.v1, acc0);
    acc0 = fma(15.0 * f12 * p.f2 * f34 * f43, p.v1, acc0);
    acc0 = fma(-(5.0 * f23 * f34 * f43), p.v1, acc0);
    acc0 = fma(8.0 * f12 * f35 * f43, p.v1, acc0);
    acc0 = fma(-(15.0 * f13 * f23 * f44), p.v1, acc0);
    acc0 = fma(10.0 * f12 * f24 * f44, p.v1, acc0);
    acc0 = fma(p.f1 * f25 * f44, p.v1, acc0);
    acc0 = fma(15.0 * f13 * f33 * f44, p.v1, acc0);
    acc0 = fma(-(10.0 * f12 * f34 * f44), p.v1, acc0);
    acc0 = fma(-(p.f1 * f35 * f44), p.v1, acc0);
    acc0 = fma(f12 * f23 * f45, p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * f24 * f45), p.v1, acc0);
    acc0 = fma(f25 * f45, p.v1, acc0);
    acc0 = fma(-(f12 * f33 * f45), p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * f34 * f45, p.v1, acc0);
    acc0 = fma(-(f35 * f45), p.v1, acc0);
    acc0 = fma(5.0 * f13 * p.f2 * f46, p.v1, acc0);
    acc0 = fma(p.f1 * f23 * f46, p.v1, acc0);
    acc0 = fma(-(2.0 * f24 * f46), p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * p.f3 * f46), p.v1, acc0);
    acc0 = fma(-(p.f1 * f33 * f46), p.v1, acc0);
    acc0 = fma(2.0 * f34 * f46, p.v1, acc0);
    acc0 = fma(-(3.0 * f12 * p.f2 * f47), p.v1, acc0);
    acc0 = fma(f23 * f47, p.v1, acc0);
    acc0 = fma(3.0 * f12 * p.f3 * f47, p.v1, acc0);
    acc0 = fma(-(f33 * f47), p.v1, acc0);
    acc0 = fma(-(f17 * f33), p.v2, acc0);
    acc0 = fma(2.0 * f16 * f34, p.v2, acc0);
    acc0 = fma(-(f15 * f35), p.v2, acc0);
    acc0 = fma(-(f16 * f33 * p.f4), p.v2, acc0);
    acc0 = fma(2.0 * f15 * f34 * p.f4, p.v2, acc0);
    acc0 = fma(-(f14 * f35 * p.f4), p.v2, acc0);
    acc0 = fma(3.0 * f17 * p.f3 * f42, p.v2, acc0);
    acc0 = fma(-(f15 * f33 * f42), p.v2, acc0);
    acc0 = fma(-(10.0 * f14 * f34 * f42), p.v2, acc0);
    acc0 = fma(8.0 * f13 * f35 * f42, p.v2, acc0);
    acc0 = fma(-(2.0 * f17 * f43), p.v2, acc0);
    acc0 = fma(-(5.0 * f16 * p.f3 * f43), p.v2, acc0);
    acc0 = fma(15.0 * f14 * f33 * f43, p.v2, acc0);
    acc0 = fma(-(8.0 * f12 * f35 * f43), p.v2, acc0);
    acc0 = fma(4.0 * f16 * f44, p.v2, acc0);
    acc0 = fma(-(15.0 * f13 * f33 * f44), p.v2, acc0);
    acc0 = fma(10.0 * f12 * f34 * f44, p.v2, acc0);
    acc0 = fma(p.f1 * f35 * f44, p.v2, acc0);
    acc0 = fma(f12 * f33 * f45, p.v2, acc0);
    acc0 = fma(-(2.0 * p.f1 * f34 * f45), p.v2, acc0);
    acc0 = fma(f35 * f45, p.v2, acc0);
    acc0 = fma(-(4.0 * f14 * f46), p.v2, acc0);
    acc0 = fma(5.0 * f13 * p.f3 * f46, p.v2, acc0);
    acc0 = fma(p.f1 * f33 * f46, p.v2, acc0);
    acc0 = fma(-(2.0 * f34 * f46), p.v2, acc0);
    acc0 = fma(2.0 * f13 * f47, p.v2, acc0);
    acc0 = fma(-(3.0 * f12 * p.f3 * f47), p.v2, acc0);
    acc0 = fma(f33 * f47, p.v2, acc0);
    acc0 = fma(f17 * f23, p.v3, acc0);
    acc0 = fma(-(2.0 * f16 * f24), p.v3, acc0);
    acc0 = fma(f15 * f25, p.v3, acc0);
    acc0 = fma(f16 * f23 * p.f4, p.v3, acc0);
    acc0 = fma(-(2.0 * f15 * f24 * p.f4), p.v3, acc0);
    acc0 = fma(f14 * f25 * p.f4, p.v3, acc0);
    acc0 = fma(-(3.0 * f17 * p.f2 * f42), p.v3, acc0);
    acc0 = fma(f15 * f23 * f42, p.v3, acc0);
    acc0 = fma(10.0 * f14 * f24 * f42, p.v3, acc0);
    acc0 = fma(-(8.0 * f13 * f25 * f42), p.v3, acc0);
    acc0 = fma(2.0 * f17 * f43, p.v3, acc0);
    acc0 = fma(5.0 * f16 * p.f2 * f43, p.v3, acc0);
    acc0 = fma(-(15.0 * f14 * f23 * f43), p.v3, acc0);
    acc0 = fma(8.0 * f12 * f25 * f43, p.v3, acc0);
    acc0 = fma(-(4.0 * f16 * f44), p.v3, acc0);
    acc0 = fma(15.0 * f13 * f23 * f44, p.v3, acc0);
    acc0 = fma(-(10.0 * f12 * f24 * f44), p.v3, acc0);
    acc0 = fma(-(p.f1 * f25 * f44), p.v3, acc0);
    acc0 = fma(-(f12 * f23 * f45), p.v3, acc0);
    acc0 = fma(2.0 * p.f1 * f24 * f45, p.v3, acc0);
    acc0 = fma(-(f25 * f45), p.v3, acc0);
    acc0 = fma(4.0 * f14 * f46, p.v3, acc0);
    acc0 = fma(-(5.0 * f13 * p.f2 * f46), p.v3, acc0);
    acc0 = fma(-(p.f1 * f23 * f46), p.v3, acc0);
    acc0 = fma(2.0 * f24 * f46, p.v3, acc0);
    acc0 = fma(-(2.0 * f13 * f47), p.v3, acc0);
    acc0 = fma(3.0 * f12 * p.f2 * f47, p.v3, acc0);
    acc0 = fma(-(f23 * f47), p.v3, acc0);
    acc0 = fma(
        -(f1mf22
            * f1mf32
            * f2mf3
            * fma(
                f12,
                fma(
                    5.0,
                    f43,
                    fma(
                        -(6.0 * (p.f2 + p.f3)),
                        f42,
                        fma(
                            2.0 * p.f2 * p.f3,
                            p.f2 + p.f3,
                            (fma(p.f2, p.f3, f22) + f32) * p.f4,
                        ),
                    ),
                ),
                fma(
                    p.f1 * fma(p.f2, p.f3, 2.0 * (p.f2 + p.f3) * p.f4),
                    fma(
                        5.0,
                        f42,
                        fma(3.0 * p.f2, p.f3, -(4.0 * (p.f2 + p.f3) * p.f4)),
                    ),
                    fma(
                        f13,
                        fma(-3.0, f42, fma(p.f2, p.f3, f22) + f32),
                        p.f2 * p.f3
                            * p.f4
                            * fma(
                                5.0,
                                f42,
                                fma(3.0 * p.f2, p.f3, -(4.0 * (p.f2 + p.f3) * p.f4)),
                            ),
                    ),
                ),
            )),
        p.v4,
        acc0,
    );
    acc0 / (f1mf22 * f1mf32 * f1mf43 * f2mf3 * f2mf42 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta3`, reconstruction 101.
pub(crate) fn delta3_101(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta3`, reconstruction 102.
pub(crate) fn delta3_102(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta3`, reconstruction 1032.
pub(crate) fn delta3_1032(p: &PolyInputs) -> f64 {
    let f1mf4 = p.f1 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf43 = f1mf42 * f1mf4;
    fma(2.0, p.v4, fma(-2.0, p.v1, fma(p.d1, f1mf4, p.d4 * f1mf4))) / f1mf43
}

/// `IMRPhenomXHM_Intermediate_Amp_delta3`, reconstruction 104.
pub(crate) fn delta3_104(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf4 = p.f1 - p.f4;
    let f2mf4 = p.f2 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 = 2.0 * f12 * f2mf42 * (p.v1 - p.v4);
    acc0 = fma(f44, -p.v1 + p.v2, acc0);
    acc0 = fma(2.0 * f22 * f42, p.v1 - p.v4, acc0);
    acc0 = fma(2.0 * f13 * p.f4, p.v2 - p.v4, acc0);
    acc0 = fma(f24, -p.v1 + p.v4, acc0);
    acc0 = fma(f14, -p.v2 + p.v4, acc0);
    acc0 = fma(
        2.0 * p.f1 * p.f4,
        fma(
            2.0 * p.f2 * p.f4,
            -p.v1 + p.v4,
            fma(f42, p.v1 - p.v2, f22 * (p.v1 - p.v4)),
        ),
        acc0,
    );
    fma(
        2.0,
        acc0,
        fma(
            p.d4 * f1mf22 * f1mf4 * f2mf4,
            fma(2.0, p.f1, p.f2) + p.f4,
            -(p.d1 * f1mf2 * f1mf4 * f2mf42 * fma(2.0, p.f4, p.f1 + p.f2)),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta3`, reconstruction 1042.
pub(crate) fn delta3_1042(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f1mf3 = p.f1 - p.f3;
    let f3mf4 = p.f3 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf32 = f1mf3 * f1mf3;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf42 * f1mf4;
    let mut acc0 = -(p.d1 * f13 * f32);
    acc0 = fma(2.0 * p.d4 * f14, p.f3, acc0);
    acc0 = fma(-(3.0 * p.d4 * f13), f32, acc0);
    acc0 = fma(p.d1 * p.f1, f34, acc0);
    acc0 = fma(p.d4 * p.f1, f34, acc0);
    acc0 = fma(-(2.0 * p.d4 * f14), p.f4, acc0);
    acc0 = fma(2.0 * p.d1 * f13 * p.f3, p.f4, acc0);
    acc0 = fma(2.0 * p.d4 * f13 * p.f3, p.f4, acc0);
    acc0 = fma(-(p.d1 * f12 * f32), p.f4, acc0);
    acc0 = fma(p.d4 * f12 * f32, p.f4, acc0);
    acc0 = fma(-(p.d1 * f34), p.f4, acc0);
    acc0 = fma(-(p.d4 * f34), p.f4, acc0);
    acc0 = fma(-(p.d1 * f13), f42, acc0);
    acc0 = fma(p.d4 * f13, f42, acc0);
    acc0 = fma(2.0 * p.d1 * f12 * p.f3, f42, acc0);
    acc0 = fma(-(2.0 * p.d4 * f12 * p.f3), f42, acc0);
    acc0 = fma(-(p.d1 * p.f1 * f32), f42, acc0);
    acc0 = fma(p.d4 * p.f1 * f32, f42, acc0);
    acc0 = fma(-(p.d1 * f12), f43, acc0);
    acc0 = fma(p.d4 * f12, f43, acc0);
    acc0 = fma(-(2.0 * p.d1 * p.f1 * p.f3), f43, acc0);
    acc0 = fma(-(2.0 * p.d4 * p.f1 * p.f3), f43, acc0);
    acc0 = fma(3.0 * p.d1 * f32, f43, acc0);
    acc0 = fma(p.d4 * f32, f43, acc0);
    acc0 = fma(2.0 * p.d1 * p.f1, f44, acc0);
    acc0 = fma(-(2.0 * p.d1 * p.f3), f44, acc0);
    acc0 = fma(4.0 * f12 * f32, p.v1, acc0);
    acc0 = fma(-(2.0 * f34), p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * p.f3 * p.f4), p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f32 * p.f4, p.v1, acc0);
    acc0 = fma(4.0 * f12 * f42, p.v1, acc0);
    acc0 = fma(-(8.0 * p.f1 * p.f3 * f42), p.v1, acc0);
    acc0 = fma(4.0 * f32 * f42, p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f43, p.v1, acc0);
    acc0 = fma(-(2.0 * f44), p.v1, acc0);
    acc0 = fma(-(2.0 * f14), p.v3, acc0);
    acc0 = fma(4.0 * f13 * p.f4, p.v3, acc0);
    acc0 = fma(-(4.0 * p.f1 * f43), p.v3, acc0);
    acc0 = fma(2.0 * f44, p.v3, acc0);
    acc0 = fma(2.0 * f14, p.v4, acc0);
    acc0 = fma(-(4.0 * f12 * f32), p.v4, acc0);
    acc0 = fma(2.0 * f34, p.v4, acc0);
    acc0 = fma(-(4.0 * f13 * p.f4), p.v4, acc0);
    acc0 = fma(8.0 * f12 * p.f3 * p.f4, p.v4, acc0);
    acc0 = fma(-(4.0 * p.f1 * f32 * p.f4), p.v4, acc0);
    acc0 = fma(-(4.0 * f12 * f42), p.v4, acc0);
    acc0 = fma(8.0 * p.f1 * p.f3 * f42, p.v4, acc0);
    acc0 = fma(-(4.0 * f32 * f42), p.v4, acc0);
    acc0 / (f1mf32 * f1mf43 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta3`, reconstruction 105.
pub(crate) fn delta3_105(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f16 = f15 * p.f1;
    let f17 = f16 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f25 = f24 * p.f2;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f35 = f34 * p.f3;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f46 = f45 * p.f4;
    let f47 = f46 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf3 = p.f1 - p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f2mf3 = p.f2 - p.f3;
    let f2mf4 = p.f2 - p.f4;
    let f3mf4 = p.f3 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f1mf32 = f1mf3 * f1mf3;
    let f2mf42 = f2mf4 * f2mf4;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 = p.d1
        * f1mf2
        * f1mf3
        * f1mf4
        * f2mf3
        * f2mf42
        * f3mf42
        * (fma(
            2.0 * (p.f1 + p.f2 + p.f3),
            p.f4,
            fma(p.f2, p.f3, fma(p.f1, p.f2, p.f1 * p.f3)),
        ) + f42);
    acc0 = fma(
        p.d4 * f1mf22 * f1mf32 * f1mf4 * f2mf3 * f2mf4 * f3mf4,
        fma(
            2.0 * p.f1,
            p.f2 + p.f3 + p.f4,
            fma(p.f2 + p.f3, p.f4, fma(p.f2, p.f3, f12)),
        ),
        acc0,
    );
    acc0 = fma(-(5.0 * f13 * f24 * f32), p.v1, acc0);
    acc0 = fma(4.0 * f12 * f25 * f32, p.v1, acc0);
    acc0 = fma(5.0 * f13 * f22 * f34, p.v1, acc0);
    acc0 = fma(-(2.0 * f25 * f34), p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * f22 * f35), p.v1, acc0);
    acc0 = fma(2.0 * f24 * f35, p.v1, acc0);
    acc0 = fma(10.0 * f13 * f24 * p.f3 * p.f4, p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * f25 * p.f3 * p.f4), p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * f24 * f32 * p.f4), p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f25 * f32 * p.f4, p.v1, acc0);
    acc0 = fma(-(10.0 * f13 * p.f2 * f34 * p.f4), p.v1, acc0);
    acc0 = fma(5.0 * f12 * f22 * f34 * p.f4, p.v1, acc0);
    acc0 = fma(8.0 * f12 * p.f2 * f35 * p.f4, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f22 * f35 * p.f4), p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * f24 * f42), p.v1, acc0);
    acc0 = fma(4.0 * f12 * f25 * f42, p.v1, acc0);
    acc0 = fma(10.0 * f12 * f24 * p.f3 * f42, p.v1, acc0);
    acc0 = fma(-(8.0 * p.f1 * f25 * p.f3 * f42), p.v1, acc0);
    acc0 = fma(-(5.0 * p.f1 * f24 * f32 * f42), p.v1, acc0);
    acc0 = fma(4.0 * f25 * f32 * f42, p.v1, acc0);
    acc0 = fma(5.0 * f13 * f34 * f42, p.v1, acc0);
    acc0 = fma(-(10.0 * f12 * p.f2 * f34 * f42), p.v1, acc0);
    acc0 = fma(5.0 * p.f1 * f22 * f34 * f42, p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * f35 * f42), p.v1, acc0);
    acc0 = fma(8.0 * p.f1 * p.f2 * f35 * f42, p.v1, acc0);
    acc0 = fma(-(4.0 * f22 * f35 * f42), p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * f24 * f43), p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f25 * f43, p.v1, acc0);
    acc0 = fma(-(20.0 * f13 * f22 * p.f3 * f43), p.v1, acc0);
    acc0 = fma(10.0 * p.f1 * f24 * p.f3 * f43, p.v1, acc0);
    acc0 = fma(20.0 * f13 * p.f2 * f32 * f43, p.v1, acc0);
    acc0 = fma(-(5.0 * f24 * f32 * f43), p.v1, acc0);
    acc0 = fma(5.0 * f12 * f34 * f43, p.v1, acc0);
    acc0 = fma(-(10.0 * p.f1 * p.f2 * f34 * f43), p.v1, acc0);
    acc0 = fma(5.0 * f22 * f34 * f43, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f35 * f43), p.v1, acc0);
    acc0 = fma(15.0 * f13 * f22 * f44, p.v1, acc0);
    acc0 = fma(-(5.0 * p.f1 * f24 * f44), p.v1, acc0);
    acc0 = fma(-(2.0 * f25 * f44), p.v1, acc0);
    acc0 = fma(-(15.0 * f13 * f32 * f44), p.v1, acc0);
    acc0 = fma(5.0 * p.f1 * f34 * f44, p.v1, acc0);
    acc0 = fma(2.0 * f35 * f44, p.v1, acc0);
    acc0 = fma(-(10.0 * f13 * p.f2 * f45), p.v1, acc0);
    acc0 = fma(-(f12 * f22 * f45), p.v1, acc0);
    acc0 = fma(3.0 * f24 * f45, p.v1, acc0);
    acc0 = fma(10.0 * f13 * p.f3 * f45, p.v1, acc0);
    acc0 = fma(f12 * f32 * f45, p.v1, acc0);
    acc0 = fma(-(3.0 * f34 * f45), p.v1, acc0);
    acc0 = fma(2.0 * f12 * p.f2 * f46, p.v1, acc0);
    acc0 = fma(-(p.f1 * f22 * f46), p.v1, acc0);
    acc0 = fma(-(2.0 * f12 * p.f3 * f46), p.v1, acc0);
    acc0 = fma(p.f1 * f32 * f46, p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * p.f2 * f47, p.v1, acc0);
    acc0 = fma(-(f22 * f47), p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * p.f3 * f47), p.v1, acc0);
    acc0 = fma(f32 * f47, p.v1, acc0);
    acc0 = fma(f17 * f32, p.v2, acc0);
    acc0 = fma(-(3.0 * f15 * f34), p.v2, acc0);
    acc0 = fma(2.0 * f14 * f35, p.v2, acc0);
    acc0 = fma(-(2.0 * f17 * p.f3 * p.f4), p.v2, acc0);
    acc0 = fma(f16 * f32 * p.f4, p.v2, acc0);
    acc0 = fma(5.0 * f14 * f34 * p.f4, p.v2, acc0);
    acc0 = fma(-(4.0 * f13 * f35 * p.f4), p.v2, acc0);
    acc0 = fma(f17 * f42, p.v2, acc0);
    acc0 = fma(-(2.0 * f16 * p.f3 * f42), p.v2, acc0);
    acc0 = fma(f15 * f32 * f42, p.v2, acc0);
    acc0 = fma(f16 * f43, p.v2, acc0);
    acc0 = fma(10.0 * f15 * p.f3 * f43, p.v2, acc0);
    acc0 = fma(-(15.0 * f14 * f32 * f43), p.v2, acc0);
    acc0 = fma(4.0 * p.f1 * f35 * f43, p.v2, acc0);
    acc0 = fma(-(8.0 * f15 * f44), p.v2, acc0);
    acc0 = fma(15.0 * f13 * f32 * f44, p.v2, acc0);
    acc0 = fma(-(5.0 * p.f1 * f34 * f44), p.v2, acc0);
    acc0 = fma(-(2.0 * f35 * f44), p.v2, acc0);
    acc0 = fma(8.0 * f14 * f45, p.v2, acc0);
    acc0 = fma(-(10.0 * f13 * p.f3 * f45), p.v2, acc0);
    acc0 = fma(-(f12 * f32 * f45), p.v2, acc0);
    acc0 = fma(3.0 * f34 * f45, p.v2, acc0);
    acc0 = fma(-(f13 * f46), p.v2, acc0);
    acc0 = fma(2.0 * f12 * p.f3 * f46, p.v2, acc0);
    acc0 = fma(-(p.f1 * f32 * f46), p.v2, acc0);
    acc0 = fma(-(f12 * f47), p.v2, acc0);
    acc0 = fma(2.0 * p.f1 * p.f3 * f47, p.v2, acc0);
    acc0 = fma(-(f32 * f47), p.v2, acc0);
    acc0 = fma(-(f17 * f22), p.v3, acc0);
    acc0 = fma(3.0 * f15 * f24, p.v3, acc0);
    acc0 = fma(-(2.0 * f14 * f25), p.v3, acc0);
    acc0 = fma(2.0 * f17 * p.f2 * p.f4, p.v3, acc0);
    acc0 = fma(-(f16 * f22 * p.f4), p.v3, acc0);
    acc0 = fma(-(5.0 * f14 * f24 * p.f4), p.v3, acc0);
    acc0 = fma(4.0 * f13 * f25 * p.f4, p.v3, acc0);
    acc0 = fma(-(f17 * f42), p.v3, acc0);
    acc0 = fma(2.0 * f16 * p.f2 * f42, p.v3, acc0);
    acc0 = fma(-(f15 * f22 * f42), p.v3, acc0);
    acc0 = fma(-(f16 * f43), p.v3, acc0);
    acc0 = fma(-(10.0 * f15 * p.f2 * f43), p.v3, acc0);
    acc0 = fma(15.0 * f14 * f22 * f43, p.v3, acc0);
    acc0 = fma(-(4.0 * p.f1 * f25 * f43), p.v3, acc0);
    acc0 = fma(8.0 * f15 * f44, p.v3, acc0);
    acc0 = fma(-(15.0 * f13 * f22 * f44), p.v3, acc0);
    acc0 = fma(5.0 * p.f1 * f24 * f44, p.v3, acc0);
    acc0 = fma(2.0 * f25 * f44, p.v3, acc0);
    acc0 = fma(-(8.0 * f14 * f45), p.v3, acc0);
    acc0 = fma(10.0 * f13 * p.f2 * f45, p.v3, acc0);
    acc0 = fma(f12 * f22 * f45, p.v3, acc0);
    acc0 = fma(-(3.0 * f24 * f45), p.v3, acc0);
    acc0 = fma(f13 * f46, p.v3, acc0);
    acc0 = fma(-(2.0 * f12 * p.f2 * f46), p.v3, acc0);
    acc0 = fma(p.f1 * f22 * f46, p.v3, acc0);
    acc0 = fma(f12 * f47, p.v3, acc0);
    acc0 = fma(-(2.0 * p.f1 * p.f2 * f47), p.v3, acc0);
    acc0 = fma(f22 * f47, p.v3, acc0);
    acc0 = fma(
        f1mf22
            * f1mf32
            * f2mf3
            * fma(
                p.f1,
                fma(
                    10.0,
                    f43,
                    fma(
                        -(3.0 * (p.f2 + p.f3)),
                        f42,
                        fma(
                            4.0 * p.f2 * p.f3,
                            p.f2 + p.f3,
                            -(4.0 * (fma(p.f2, p.f3, f22) + f32) * p.f4),
                        ),
                    ),
                ),
                fma(
                    5.0 * (p.f2 + p.f3),
                    f43,
                    fma(
                        -(4.0 * (fma(p.f2, p.f3, f22) + f32)),
                        f42,
                        fma(
                            f12 * fma(-2.0, p.f4, p.f2 + p.f3),
                            fma(2.0, p.f2 + p.f3, p.f4),
                            fma(2.0 * f22, f32, f13 * fma(-2.0, p.f4, p.f2 + p.f3)),
                        ),
                    ),
                ),
            ),
        p.v4,
        acc0,
    );
    acc0 / (f1mf22 * f1mf32 * f1mf43 * f2mf3 * f2mf42 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta4`, reconstruction 101.
pub(crate) fn delta4_101(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta4`, reconstruction 102.
pub(crate) fn delta4_102(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta4`, reconstruction 1032.
pub(crate) fn delta4_1032(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta4`, reconstruction 104.
pub(crate) fn delta4_104(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf4 = p.f1 - p.f4;
    let f2mf4 = p.f2 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 = p.d1 * f1mf2 * f1mf4 * f2mf42;
    acc0 = fma(-(p.d4 * f1mf22 * f1mf4), f2mf4, acc0);
    acc0 = fma(-(3.0 * p.f1 * f22), p.v1, acc0);
    acc0 = fma(2.0 * f23, p.v1, acc0);
    acc0 = fma(6.0 * p.f1 * p.f2 * p.f4, p.v1, acc0);
    acc0 = fma(-(3.0 * f22 * p.f4), p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f42), p.v1, acc0);
    acc0 = fma(f43, p.v1, acc0);
    acc0 = fma(f13, p.v2, acc0);
    acc0 = fma(-(3.0 * f12 * p.f4), p.v2, acc0);
    acc0 = fma(3.0 * p.f1 * f42, p.v2, acc0);
    acc0 = fma(-f43, p.v2, acc0);
    acc0 = fma(
        -(f1mf22 * fma(-3.0, p.f4, fma(2.0, p.f2, p.f1))),
        p.v4,
        acc0,
    );
    acc0 / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta4`, reconstruction 1042.
pub(crate) fn delta4_1042(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f1mf3 = p.f1 - p.f3;
    let f3mf4 = p.f3 - p.f4;
    let f1mf42 = f1mf4 * f1mf4;
    let f1mf32 = f1mf3 * f1mf3;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf42 * f1mf4;
    let mut acc0 = p.d1 * f1mf3 * f1mf4 * f3mf42;
    acc0 = fma(-(p.d4 * f1mf32 * f1mf4), f3mf4, acc0);
    acc0 = fma(-(3.0 * p.f1 * f32), p.v1, acc0);
    acc0 = fma(2.0 * f33, p.v1, acc0);
    acc0 = fma(6.0 * p.f1 * p.f3 * p.f4, p.v1, acc0);
    acc0 = fma(-(3.0 * f32 * p.f4), p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f42), p.v1, acc0);
    acc0 = fma(f43, p.v1, acc0);
    acc0 = fma(f13, p.v3, acc0);
    acc0 = fma(-(3.0 * f12 * p.f4), p.v3, acc0);
    acc0 = fma(3.0 * p.f1 * f42, p.v3, acc0);
    acc0 = fma(-f43, p.v3, acc0);
    acc0 = fma(
        -(f1mf32 * fma(-3.0, p.f4, fma(2.0, p.f3, p.f1))),
        p.v4,
        acc0,
    );
    acc0 / (f1mf32 * f1mf43 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta4`, reconstruction 105.
pub(crate) fn delta4_105(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f16 = f15 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f25 = f24 * p.f2;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f35 = f34 * p.f3;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f46 = f45 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf3 = p.f1 - p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f2mf3 = p.f2 - p.f3;
    let f2mf4 = p.f2 - p.f4;
    let f3mf4 = p.f3 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f1mf32 = f1mf3 * f1mf3;
    let f2mf42 = f2mf4 * f2mf4;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 = -(p.d1
        * f1mf2
        * f1mf3
        * f1mf4
        * f2mf3
        * f2mf42
        * f3mf42
        * fma(2.0, p.f4, p.f1 + p.f2 + p.f3));
    acc0 = fma(
        -(p.d4 * f1mf22 * f1mf32 * f1mf4 * f2mf3 * f2mf4 * f3mf4),
        fma(2.0, p.f1, p.f2) + p.f3 + p.f4,
        acc0,
    );
    acc0 = fma(5.0 * f13 * f23 * f32, p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f25 * f32), p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * f22 * f33), p.v1, acc0);
    acc0 = fma(2.0 * f25 * f33, p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f22 * f35, p.v1, acc0);
    acc0 = fma(-(2.0 * f23 * f35), p.v1, acc0);
    acc0 = fma(-(10.0 * f13 * f23 * p.f3 * p.f4), p.v1, acc0);
    acc0 = fma(6.0 * p.f1 * f25 * p.f3 * p.f4, p.v1, acc0);
    acc0 = fma(5.0 * f12 * f23 * f32 * p.f4, p.v1, acc0);
    acc0 = fma(-(3.0 * f25 * f32 * p.f4), p.v1, acc0);
    acc0 = fma(10.0 * f13 * p.f2 * f33 * p.f4, p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * f22 * f33 * p.f4), p.v1, acc0);
    acc0 = fma(-(6.0 * p.f1 * p.f2 * f35 * p.f4), p.v1, acc0);
    acc0 = fma(3.0 * f22 * f35 * p.f4, p.v1, acc0);
    acc0 = fma(5.0 * f13 * f23 * f42, p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f25 * f42), p.v1, acc0);
    acc0 = fma(15.0 * f13 * f22 * p.f3 * f42, p.v1, acc0);
    acc0 = fma(-(10.0 * f12 * f23 * p.f3 * f42), p.v1, acc0);
    acc0 = fma(-(15.0 * f13 * p.f2 * f32 * f42), p.v1, acc0);
    acc0 = fma(5.0 * p.f1 * f23 * f32 * f42, p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * f33 * f42), p.v1, acc0);
    acc0 = fma(10.0 * f12 * p.f2 * f33 * f42, p.v1, acc0);
    acc0 = fma(-(5.0 * p.f1 * f22 * f33 * f42), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f35 * f42, p.v1, acc0);
    acc0 = fma(-(10.0 * f13 * f22 * f43), p.v1, acc0);
    acc0 = fma(5.0 * f12 * f23 * f43, p.v1, acc0);
    acc0 = fma(f25 * f43, p.v1, acc0);
    acc0 = fma(15.0 * f12 * f22 * p.f3 * f43, p.v1, acc0);
    acc0 = fma(-(10.0 * p.f1 * f23 * p.f3 * f43), p.v1, acc0);
    acc0 = fma(10.0 * f13 * f32 * f43, p.v1, acc0);
    acc0 = fma(-(15.0 * f12 * p.f2 * f32 * f43), p.v1, acc0);
    acc0 = fma(5.0 * f23 * f32 * f43, p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * f33 * f43), p.v1, acc0);
    acc0 = fma(10.0 * p.f1 * p.f2 * f33 * f43, p.v1, acc0);
    acc0 = fma(-(5.0 * f22 * f33 * f43), p.v1, acc0);
    acc0 = fma(-(f35 * f43), p.v1, acc0);
    acc0 = fma(5.0 * f13 * p.f2 * f44, p.v1, acc0);
    acc0 = fma(-(10.0 * f12 * f22 * f44), p.v1, acc0);
    acc0 = fma(5.0 * p.f1 * f23 * f44, p.v1, acc0);
    acc0 = fma(-(5.0 * f13 * p.f3 * f44), p.v1, acc0);
    acc0 = fma(10.0 * f12 * f32 * f44, p.v1, acc0);
    acc0 = fma(-(5.0 * p.f1 * f33 * f44), p.v1, acc0);
    acc0 = fma(5.0 * f12 * p.f2 * f45, p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * f22 * f45, p.v1, acc0);
    acc0 = fma(-(3.0 * f23 * f45), p.v1, acc0);
    acc0 = fma(-(5.0 * f12 * p.f3 * f45), p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * f32 * f45), p.v1, acc0);
    acc0 = fma(3.0 * f33 * f45, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * p.f2 * f46), p.v1, acc0);
    acc0 = fma(2.0 * f22 * f46, p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * p.f3 * f46, p.v1, acc0);
    acc0 = fma(-(2.0 * f32 * f46), p.v1, acc0);
    acc0 = fma(-(2.0 * f16 * f32), p.v2, acc0);
    acc0 = fma(3.0 * f15 * f33, p.v2, acc0);
    acc0 = fma(-(f13 * f35), p.v2, acc0);
    acc0 = fma(4.0 * f16 * p.f3 * p.f4, p.v2, acc0);
    acc0 = fma(-(2.0 * f15 * f32 * p.f4), p.v2, acc0);
    acc0 = fma(-(5.0 * f14 * f33 * p.f4), p.v2, acc0);
    acc0 = fma(3.0 * f12 * f35 * p.f4, p.v2, acc0);
    acc0 = fma(-(2.0 * f16 * f42), p.v2, acc0);
    acc0 = fma(-(5.0 * f15 * p.f3 * f42), p.v2, acc0);
    acc0 = fma(10.0 * f14 * f32 * f42, p.v2, acc0);
    acc0 = fma(-(3.0 * p.f1 * f35 * f42), p.v2, acc0);
    acc0 = fma(4.0 * f15 * f43, p.v2, acc0);
    acc0 = fma(-(5.0 * f14 * p.f3 * f43), p.v2, acc0);
    acc0 = fma(f35 * f43, p.v2, acc0);
    acc0 = fma(5.0 * f13 * p.f3 * f44, p.v2, acc0);
    acc0 = fma(-(10.0 * f12 * f32 * f44), p.v2, acc0);
    acc0 = fma(5.0 * p.f1 * f33 * f44, p.v2, acc0);
    acc0 = fma(-(4.0 * f13 * f45), p.v2, acc0);
    acc0 = fma(5.0 * f12 * p.f3 * f45, p.v2, acc0);
    acc0 = fma(2.0 * p.f1 * f32 * f45, p.v2, acc0);
    acc0 = fma(-(3.0 * f33 * f45), p.v2, acc0);
    acc0 = fma(2.0 * f12 * f46, p.v2, acc0);
    acc0 = fma(-(4.0 * p.f1 * p.f3 * f46), p.v2, acc0);
    acc0 = fma(2.0 * f32 * f46, p.v2, acc0);
    acc0 = fma(2.0 * f16 * f22, p.v3, acc0);
    acc0 = fma(-(3.0 * f15 * f23), p.v3, acc0);
    acc0 = fma(f13 * f25, p.v3, acc0);
    acc0 = fma(-(4.0 * f16 * p.f2 * p.f4), p.v3, acc0);
    acc0 = fma(2.0 * f15 * f22 * p.f4, p.v3, acc0);
    acc0 = fma(5.0 * f14 * f23 * p.f4, p.v3, acc0);
    acc0 = fma(-(3.0 * f12 * f25 * p.f4), p.v3, acc0);
    acc0 = fma(2.0 * f16 * f42, p.v3, acc0);
    acc0 = fma(5.0 * f15 * p.f2 * f42, p.v3, acc0);
    acc0 = fma(-(10.0 * f14 * f22 * f42), p.v3, acc0);
    acc0 = fma(3.0 * p.f1 * f25 * f42, p.v3, acc0);
    acc0 = fma(-(4.0 * f15 * f43), p.v3, acc0);
    acc0 = fma(5.0 * f14 * p.f2 * f43, p.v3, acc0);
    acc0 = fma(-(f25 * f43), p.v3, acc0);
    acc0 = fma(-(5.0 * f13 * p.f2 * f44), p.v3, acc0);
    acc0 = fma(10.0 * f12 * f22 * f44, p.v3, acc0);
    acc0 = fma(-(5.0 * p.f1 * f23 * f44), p.v3, acc0);
    acc0 = fma(4.0 * f13 * f45, p.v3, acc0);
    acc0 = fma(-(5.0 * f12 * p.f2 * f45), p.v3, acc0);
    acc0 = fma(-(2.0 * p.f1 * f22 * f45), p.v3, acc0);
    acc0 = fma(3.0 * f23 * f45, p.v3, acc0);
    acc0 = fma(-(2.0 * f12 * f46), p.v3, acc0);
    acc0 = fma(4.0 * p.f1 * p.f2 * f46, p.v3, acc0);
    acc0 = fma(-(2.0 * f22 * f46), p.v3, acc0);
    acc0 = fma(
        -(f1mf22
            * f1mf32
            * f2mf3
            * fma(
                5.0,
                f43,
                fma(
                    p.f1,
                    fma(
                        5.0,
                        f42,
                        fma(
                            -(6.0 * (p.f2 + p.f3)),
                            p.f4,
                            fma(5.0 * p.f2, p.f3, f22) + f32,
                        ),
                    ),
                    fma(
                        -(3.0 * (fma(p.f2, p.f3, f22) + f32)),
                        p.f4,
                        fma(
                            2.0 * p.f2 * p.f3,
                            p.f2 + p.f3,
                            2.0 * f12 * fma(-2.0, p.f4, p.f2 + p.f3),
                        ),
                    ),
                ),
            )),
        p.v4,
        acc0,
    );
    acc0 / (f1mf22 * f1mf32 * f1mf43 * f2mf3 * f2mf42 * f3mf42)
}

/// `IMRPhenomXHM_Intermediate_Amp_delta5`, reconstruction 101.
pub(crate) fn delta5_101(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta5`, reconstruction 102.
pub(crate) fn delta5_102(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta5`, reconstruction 1032.
pub(crate) fn delta5_1032(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta5`, reconstruction 104.
pub(crate) fn delta5_104(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta5`, reconstruction 1042.
pub(crate) fn delta5_1042(_: &PolyInputs) -> f64 {
    0.0
}

/// `IMRPhenomXHM_Intermediate_Amp_delta5`, reconstruction 105.
pub(crate) fn delta5_105(p: &PolyInputs) -> f64 {
    let f12 = p.f1 * p.f1;
    let f13 = f12 * p.f1;
    let f14 = f13 * p.f1;
    let f15 = f14 * p.f1;
    let f22 = p.f2 * p.f2;
    let f23 = f22 * p.f2;
    let f24 = f23 * p.f2;
    let f32 = p.f3 * p.f3;
    let f33 = f32 * p.f3;
    let f34 = f33 * p.f3;
    let f42 = p.f4 * p.f4;
    let f43 = f42 * p.f4;
    let f44 = f43 * p.f4;
    let f45 = f44 * p.f4;
    let f1mf2 = p.f1 - p.f2;
    let f1mf3 = p.f1 - p.f3;
    let f1mf4 = p.f1 - p.f4;
    let f2mf3 = p.f2 - p.f3;
    let f2mf4 = p.f2 - p.f4;
    let f3mf4 = p.f3 - p.f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f1mf32 = f1mf3 * f1mf3;
    let f2mf42 = f2mf4 * f2mf4;
    let f3mf42 = f3mf4 * f3mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    let mut acc0 = p.d1 * f1mf2 * f1mf3 * f1mf4 * f2mf3 * f2mf42 * f3mf42;
    acc0 = fma(p.d4 * f1mf22 * f1mf32 * f1mf4 * f2mf3 * f2mf4, f3mf4, acc0);
    acc0 = fma(-(4.0 * f12 * f23 * f32), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f24 * f32, p.v1, acc0);
    acc0 = fma(4.0 * f12 * f22 * f33, p.v1, acc0);
    acc0 = fma(-(2.0 * f24 * f33), p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f22 * f34), p.v1, acc0);
    acc0 = fma(2.0 * f23 * f34, p.v1, acc0);
    acc0 = fma(8.0 * f12 * f23 * p.f3 * p.f4, p.v1, acc0);
    acc0 = fma(-(6.0 * p.f1 * f24 * p.f3 * p.f4), p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f23 * f32 * p.f4), p.v1, acc0);
    acc0 = fma(3.0 * f24 * f32 * p.f4, p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * p.f2 * f33 * p.f4), p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f22 * f33 * p.f4, p.v1, acc0);
    acc0 = fma(6.0 * p.f1 * p.f2 * f34 * p.f4, p.v1, acc0);
    acc0 = fma(-(3.0 * f22 * f34 * p.f4), p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * f23 * f42), p.v1, acc0);
    acc0 = fma(3.0 * p.f1 * f24 * f42, p.v1, acc0);
    acc0 = fma(-(12.0 * f12 * f22 * p.f3 * f42), p.v1, acc0);
    acc0 = fma(8.0 * p.f1 * f23 * p.f3 * f42, p.v1, acc0);
    acc0 = fma(12.0 * f12 * p.f2 * f32 * f42, p.v1, acc0);
    acc0 = fma(-(4.0 * f23 * f32 * f42), p.v1, acc0);
    acc0 = fma(4.0 * f12 * f33 * f42, p.v1, acc0);
    acc0 = fma(-(8.0 * p.f1 * p.f2 * f33 * f42), p.v1, acc0);
    acc0 = fma(4.0 * f22 * f33 * f42, p.v1, acc0);
    acc0 = fma(-(3.0 * p.f1 * f34 * f42), p.v1, acc0);
    acc0 = fma(8.0 * f12 * f22 * f43, p.v1, acc0);
    acc0 = fma(-(4.0 * p.f1 * f23 * f43), p.v1, acc0);
    acc0 = fma(-(f24 * f43), p.v1, acc0);
    acc0 = fma(-(8.0 * f12 * f32 * f43), p.v1, acc0);
    acc0 = fma(4.0 * p.f1 * f33 * f43, p.v1, acc0);
    acc0 = fma(f34 * f43, p.v1, acc0);
    acc0 = fma(-(4.0 * f12 * p.f2 * f44), p.v1, acc0);
    acc0 = fma(-(p.f1 * f22 * f44), p.v1, acc0);
    acc0 = fma(2.0 * f23 * f44, p.v1, acc0);
    acc0 = fma(4.0 * f12 * p.f3 * f44, p.v1, acc0);
    acc0 = fma(p.f1 * f32 * f44, p.v1, acc0);
    acc0 = fma(-(2.0 * f33 * f44), p.v1, acc0);
    acc0 = fma(2.0 * p.f1 * p.f2 * f45, p.v1, acc0);
    acc0 = fma(-(f22 * f45), p.v1, acc0);
    acc0 = fma(-(2.0 * p.f1 * p.f3 * f45), p.v1, acc0);
    acc0 = fma(f32 * f45, p.v1, acc0);
    acc0 = fma(f15 * f32, p.v2, acc0);
    acc0 = fma(-(2.0 * f14 * f33), p.v2, acc0);
    acc0 = fma(f13 * f34, p.v2, acc0);
    acc0 = fma(-(2.0 * f15 * p.f3 * p.f4), p.v2, acc0);
    acc0 = fma(f14 * f32 * p.f4, p.v2, acc0);
    acc0 = fma(4.0 * f13 * f33 * p.f4, p.v2, acc0);
    acc0 = fma(-(3.0 * f12 * f34 * p.f4), p.v2, acc0);
    acc0 = fma(f15 * f42, p.v2, acc0);
    acc0 = fma(4.0 * f14 * p.f3 * f42, p.v2, acc0);
    acc0 = fma(-(8.0 * f13 * f32 * f42), p.v2, acc0);
    acc0 = fma(3.0 * p.f1 * f34 * f42, p.v2, acc0);
    acc0 = fma(-(3.0 * f14 * f43), p.v2, acc0);
    acc0 = fma(8.0 * f12 * f32 * f43, p.v2, acc0);
    acc0 = fma(-(4.0 * p.f1 * f33 * f43), p.v2, acc0);
    acc0 = fma(-(f34 * f43), p.v2, acc0);
    acc0 = fma(3.0 * f13 * f44, p.v2, acc0);
    acc0 = fma(-(4.0 * f12 * p.f3 * f44), p.v2, acc0);
    acc0 = fma(-(p.f1 * f32 * f44), p.v2, acc0);
    acc0 = fma(2.0 * f33 * f44, p.v2, acc0);
    acc0 = fma(-(f12 * f45), p.v2, acc0);
    acc0 = fma(2.0 * p.f1 * p.f3 * f45, p.v2, acc0);
    acc0 = fma(-(f32 * f45), p.v2, acc0);
    acc0 = fma(-(f15 * f22), p.v3, acc0);
    acc0 = fma(2.0 * f14 * f23, p.v3, acc0);
    acc0 = fma(-(f13 * f24), p.v3, acc0);
    acc0 = fma(2.0 * f15 * p.f2 * p.f4, p.v3, acc0);
    acc0 = fma(-(f14 * f22 * p.f4), p.v3, acc0);
    acc0 = fma(-(4.0 * f13 * f23 * p.f4), p.v3, acc0);
    acc0 = fma(3.0 * f12 * f24 * p.f4, p.v3, acc0);
    acc0 = fma(-(f15 * f42), p.v3, acc0);
    acc0 = fma(-(4.0 * f14 * p.f2 * f42), p.v3, acc0);
    acc0 = fma(8.0 * f13 * f22 * f42, p.v3, acc0);
    acc0 = fma(-(3.0 * p.f1 * f24 * f42), p.v3, acc0);
    acc0 = fma(3.0 * f14 * f43, p.v3, acc0);
    acc0 = fma(-(8.0 * f12 * f22 * f43), p.v3, acc0);
    acc0 = fma(4.0 * p.f1 * f23 * f43, p.v3, acc0);
    acc0 = fma(f24 * f43, p.v3, acc0);
    acc0 = fma(-(3.0 * f13 * f44), p.v3, acc0);
    acc0 = fma(4.0 * f12 * p.f2 * f44, p.v3, acc0);
    acc0 = fma(p.f1 * f22 * f44, p.v3, acc0);
    acc0 = fma(-(2.0 * f23 * f44), p.v3, acc0);
    acc0 = fma(f12 * f45, p.v3, acc0);
    acc0 = fma(-(2.0 * p.f1 * p.f2 * f45), p.v3, acc0);
    acc0 = fma(f22 * f45, p.v3, acc0);
    acc0 = fma(
        f1mf22
            * f1mf32
            * f2mf3
            * fma(
                4.0,
                f42,
                fma(
                    -(3.0 * (p.f2 + p.f3)),
                    p.f4,
                    fma(2.0 * p.f2, p.f3, p.f1 * fma(-2.0, p.f4, p.f2 + p.f3)),
                ),
            ),
        p.v4,
        acc0,
    );
    acc0 / (f1mf22 * f1mf32 * f1mf43 * f2mf3 * f2mf42 * f3mf42)
}
