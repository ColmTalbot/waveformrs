// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Intermediate region of the (2,2) mode, `LALSimIMRPhenomX_intermediate.c`.
//!
//! Only the default versions are implemented: `IntAmpVersion = 104`
//! (fourth order polynomial for the inverse amplitude) and
//! `IntPhaseVersion = 105` (five coefficient phase derivative ansatz).
//! Collocation fits use the effective spin `STotR`.

use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::{fma, UsefulPowers};
use crate::imrphenomx::phase::PhaseCoefficients;

/* ******************************** Amplitude ******************************** */

/// `IMRPhenomX_Intermediate_Amp_22_vA`: value of the amplitude collocation
/// point at the centre of the intermediate region (IntAmpVersion 104).
pub(crate) fn amp_22_va(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let no_spin = fma(
        -19837.908020966777,
        eta3,
        fma(
            27563.641024162127,
            eta2,
            fma(1974.6112656679577, eta, 1.4873184918202145),
        ),
    ) / fma(458.4097306093354, eta2, fma(143.29004876335128, eta, 1.0));
    let eq_spin = s * fma(
        eta3,
        fma(
            1486.4658089990298,
            s2,
            fma(-6962.675551371755, s, 1612.2681322644232),
        ),
        fma(
            -19.38970173389662,
            s2,
            fma(
                3011.446602208493 * eta2,
                s,
                fma(
                    3.2646808851249016,
                    s,
                    fma(
                        eta,
                        fma(-260.3494489873286, s, -365.55631765202895),
                        27.952730865904343,
                    ),
                ),
            ),
        ),
    ) / fma(1.0, s2, fma(-10.540154508599963, s, 12.647425554323242));
    let uneq_spin = dchi * delta * fma(-296.473359655246, eta, -0.016404056649860943) * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Intermediate_Amp_22_delta0` (IntAmpVersion 104).
#[allow(clippy::too_many_arguments)]
pub(crate) fn amp_22_delta0(
    d1: f64,
    d4: f64,
    v1: f64,
    v2: f64,
    v4: f64,
    f1: f64,
    f2: f64,
    f4: f64,
) -> f64 {
    let f12 = f1 * f1;
    let f42 = f4 * f4;
    let f1mf2 = f1 - f2;
    let f1mf4 = f1 - f4;
    let f2mf4 = f2 - f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    fma(
        f12 * f1mf22
            * f2
            * fma(
                4.0,
                f42,
                fma(-(3.0 * f2), f4, fma(f1, f2, -(2.0 * f1 * f4))),
            ),
        v4,
        fma(
            f42,
            fma(
                f2 * f2mf42 * fma(-f2, f4, fma(2.0 * f1, f4, fma(-4.0, f12, 3.0 * f1 * f2))),
                v1,
                f12 * f1mf43 * v2,
            ),
            fma(
                -(d4 * f12 * f1mf22 * f1mf4 * f2 * f2mf4),
                f4,
                d1 * f1 * f1mf2 * f1mf4 * f2 * f2mf42 * f42,
            ),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomX_Intermediate_Amp_22_delta1` (IntAmpVersion 104).
#[allow(clippy::too_many_arguments)]
pub(crate) fn amp_22_delta1(
    d1: f64,
    d4: f64,
    v1: f64,
    v2: f64,
    v4: f64,
    f1: f64,
    f2: f64,
    f4: f64,
) -> f64 {
    let f12 = f1 * f1;
    let f13 = f12 * f1;
    let f14 = f13 * f1;
    let f22 = f2 * f2;
    let f23 = f22 * f2;
    let f24 = f23 * f2;
    let f42 = f4 * f4;
    let f43 = f42 * f4;
    let f44 = f43 * f4;
    let f1mf2 = f1 - f2;
    let f1mf4 = f1 - f4;
    let f2mf4 = f2 - f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    fma(
        d4 * f1 * f1mf22 * f1mf4 * f2mf4,
        fma(2.0 * f2, f4, f1 * (f2 + f4)),
        f4 * fma(
            -(d1 * f1mf2 * f1mf4 * f2mf42),
            fma(2.0 * f1, f2, (f1 + f2) * f4),
            -(2.0
                * f1
                * fma(
                    f1,
                    fma(
                        4.0 * f23,
                        -v1 + v4,
                        fma(2.0 * f43, -v1 + v2, 6.0 * f22 * f4 * (v1 - v4)),
                    ),
                    fma(
                        2.0 * f13 * f4,
                        -v2 + v4,
                        fma(
                            4.0 * f23 * f4,
                            -v1 + v4,
                            fma(f14, v2 - v4, fma(f44, v1 - v2, 3.0 * f24 * (v1 - v4))),
                        ),
                    ),
                )),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomX_Intermediate_Amp_22_delta2` (IntAmpVersion 104).
#[allow(clippy::too_many_arguments)]
pub(crate) fn amp_22_delta2(
    d1: f64,
    d4: f64,
    v1: f64,
    v2: f64,
    v4: f64,
    f1: f64,
    f2: f64,
    f4: f64,
) -> f64 {
    let f12 = f1 * f1;
    let f13 = f12 * f1;
    let f14 = f13 * f1;
    let f15 = f14 * f1;
    let f22 = f2 * f2;
    let f23 = f22 * f2;
    let f24 = f23 * f2;
    let f42 = f4 * f4;
    let f43 = f42 * f4;
    let f44 = f43 * f4;
    let f45 = f44 * f4;
    let f1mf2 = f1 - f2;
    let f1mf4 = f1 - f4;
    let f2mf4 = f2 - f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    fma(
        -(f1mf22
            * fma(
                f1 * fma(3.0, f2, -(4.0 * f4)),
                fma(2.0, f4, f2),
                fma(
                    f12,
                    fma(2.0, f2, f4),
                    fma(f2 * fma(3.0, f2, -(4.0 * f4)), f4, f13),
                ),
            )),
        v4,
        fma(
            -f45,
            v2,
            fma(
                -(f1 * f44),
                v2,
                fma(
                    8.0 * f12 * f43,
                    v2,
                    fma(
                        -(8.0 * f13 * f42),
                        v2,
                        fma(
                            f14 * f4,
                            v2,
                            fma(
                                f15,
                                v2,
                                fma(
                                    f45,
                                    v1,
                                    fma(
                                        f1 * f44,
                                        v1,
                                        fma(
                                            -(8.0 * f12 * f43),
                                            v1,
                                            fma(
                                                -(4.0 * f23 * f42),
                                                v1,
                                                fma(
                                                    12.0 * f12 * f2 * f42,
                                                    v1,
                                                    fma(
                                                        3.0 * f24 * f4,
                                                        v1,
                                                        fma(
                                                            -(4.0 * f1 * f23 * f4),
                                                            v1,
                                                            fma(
                                                                3.0 * f1 * f24,
                                                                v1,
                                                                fma(
                                                                    -(4.0 * f12 * f23),
                                                                    v1,
                                                                    fma(
                                                                        -(d4 * f1mf22
                                                                            * f1mf4
                                                                            * f2mf4),
                                                                        fma(
                                                                            2.0 * f1,
                                                                            f2 + f4,
                                                                            fma(f2, f4, f12),
                                                                        ),
                                                                        d1 * f1mf2
                                                                            * f1mf4
                                                                            * f2mf42
                                                                            * (fma(
                                                                                f1,
                                                                                f2,
                                                                                2.0 * (f1 + f2)
                                                                                    * f4,
                                                                            ) + f42),
                                                                    ),
                                                                ),
                                                            ),
                                                        ),
                                                    ),
                                                ),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomX_Intermediate_Amp_22_delta3` (IntAmpVersion 104).
#[allow(clippy::too_many_arguments)]
pub(crate) fn amp_22_delta3(
    d1: f64,
    d4: f64,
    v1: f64,
    v2: f64,
    v4: f64,
    f1: f64,
    f2: f64,
    f4: f64,
) -> f64 {
    let f12 = f1 * f1;
    let f13 = f12 * f1;
    let f14 = f13 * f1;
    let f22 = f2 * f2;
    let f23 = f22 * f2;
    let f24 = f23 * f2;
    let f42 = f4 * f4;
    let f43 = f42 * f4;
    let f44 = f43 * f4;
    let f1mf2 = f1 - f2;
    let f1mf4 = f1 - f4;
    let f2mf4 = f2 - f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    fma(
        2.0,
        fma(
            2.0 * f1 * f4,
            fma(2.0 * f2 * f4, -v1 + v4, fma(f42, v1 - v2, f22 * (v1 - v4))),
            fma(
                f14,
                -v2 + v4,
                fma(
                    f24,
                    -v1 + v4,
                    fma(
                        2.0 * f13 * f4,
                        v2 - v4,
                        fma(
                            2.0 * f22 * f42,
                            v1 - v4,
                            fma(f44, -v1 + v2, 2.0 * f12 * f2mf42 * (v1 - v4)),
                        ),
                    ),
                ),
            ),
        ),
        fma(
            d4 * f1mf22 * f1mf4 * f2mf4,
            fma(2.0, f1, f2) + f4,
            -(d1 * f1mf2 * f1mf4 * f2mf42 * fma(2.0, f4, f1 + f2)),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomX_Intermediate_Amp_22_delta4` (IntAmpVersion 104).
#[allow(clippy::too_many_arguments)]
pub(crate) fn amp_22_delta4(
    d1: f64,
    d4: f64,
    v1: f64,
    v2: f64,
    v4: f64,
    f1: f64,
    f2: f64,
    f4: f64,
) -> f64 {
    let f12 = f1 * f1;
    let f13 = f12 * f1;
    let f22 = f2 * f2;
    let f23 = f22 * f2;
    let f42 = f4 * f4;
    let f43 = f42 * f4;
    let f1mf2 = f1 - f2;
    let f1mf4 = f1 - f4;
    let f2mf4 = f2 - f4;
    let f1mf22 = f1mf2 * f1mf2;
    let f2mf42 = f2mf4 * f2mf4;
    let f1mf43 = f1mf4 * f1mf4 * f1mf4;
    fma(
        -(f1mf22 * fma(-3.0, f4, fma(2.0, f2, f1))),
        v4,
        fma(
            -f43,
            v2,
            fma(
                3.0 * f1 * f42,
                v2,
                fma(
                    -(3.0 * f12 * f4),
                    v2,
                    fma(
                        f13,
                        v2,
                        fma(
                            f43,
                            v1,
                            fma(
                                -(3.0 * f1 * f42),
                                v1,
                                fma(
                                    -(3.0 * f22 * f4),
                                    v1,
                                    fma(
                                        6.0 * f1 * f2 * f4,
                                        v1,
                                        fma(
                                            2.0 * f23,
                                            v1,
                                            fma(
                                                -(3.0 * f1 * f22),
                                                v1,
                                                fma(
                                                    -(d4 * f1mf22 * f1mf4),
                                                    f2mf4,
                                                    d1 * f1mf2 * f1mf4 * f2mf42,
                                                ),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    ) / (f1mf22 * f1mf43 * f2mf42)
}

/// `IMRPhenomX_Intermediate_Amp_22_delta0` to `_delta5` for IntAmpVersion
/// 104: coefficients of the fourth order polynomial for the inverse
/// amplitude given the boundary derivatives `d1`, `d4` and the inverse
/// amplitudes `v1`, `v2`, `v4` at `f1`, `f2`, `f4`.
///
/// The value `v3` at `f3` is unused for this version and `delta5 = 0`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn amp_22_deltas(
    d1: f64,
    d4: f64,
    v1: f64,
    v2: f64,
    v4: f64,
    f1: f64,
    f2: f64,
    f4: f64,
) -> [f64; 6] {
    [
        amp_22_delta0(d1, d4, v1, v2, v4, f1, f2, f4),
        amp_22_delta1(d1, d4, v1, v2, v4, f1, f2, f4),
        amp_22_delta2(d1, d4, v1, v2, v4, f1, f2, f4),
        amp_22_delta3(d1, d4, v1, v2, v4, f1, f2, f4),
        amp_22_delta4(d1, d4, v1, v2, v4, f1, f2, f4),
        0.0,
    ]
}

/// `IMRPhenomX_Intermediate_Amp_22_Ansatz` (IntAmpVersion 104):
/// f^(7/6) divided by a fourth order polynomial.
#[inline]
pub(crate) fn amp_22_ansatz(mf: f64, powers: &UsefulPowers, amp: &AmpCoefficients) -> f64 {
    let polynomial =
        amp.delta0 + mf * (amp.delta1 + mf * (amp.delta2 + mf * (amp.delta3 + mf * amp.delta4)));
    powers.seven_sixths / polynomial
}

/* ********************************** Phase ********************************** */

/// `IMRPhenomX_Intermediate_Phase_22_v2` (IntPhaseVersion 105): direct fit
/// to the second phase collocation point.
pub(crate) fn phase_22_v2(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let s4 = s3 * s;
    let no_spin = fma(
        1.4667258334378073e9,
        eta3,
        fma(
            -3.5225742421184325e8,
            eta2,
            fma(-5.58197349185435e6, eta, -82.54500000000004),
        ),
    ) / fma(
        2.5176585751772933e6,
        eta3,
        fma(5.385164380400193e6, eta2, fma(66757.12830903867, eta, 1.0)),
    );
    let eq_spin = s * fma(
        eta,
        fma(
            37.462217031512786,
            s4,
            fma(
                29.071788188638315,
                s3,
                fma(
                    -17.671713040498304,
                    s2,
                    fma(-132.88417400679026, s, 207.69898051583655),
                ),
            ),
        ),
        fma(
            4.984750041013893,
            s4,
            fma(
                eta3,
                fma(
                    -866.283566780576,
                    s4,
                    fma(
                        1897.5444343138167,
                        s3,
                        fma(
                            415.7134909564443,
                            s2,
                            fma(1152.425940764218, s, -1365.1499998427248),
                        ),
                    ),
                ),
                fma(
                    5.95010003393006,
                    s3,
                    fma(
                        eta2,
                        fma(
                            -647.8103976942541,
                            s3,
                            fma(-107.41099349364234, s, 170.97203068800542),
                        ),
                        fma(
                            -0.8612656616290079,
                            s2,
                            fma(-36.066611959079935, s, 19.416719811164853),
                        ),
                    ),
                ),
            ),
        ),
    ) / fma(1.0, s, -1.1492259468169692);
    let uneq_spin = dchi
        * delta
        * eta3
        * fma(
            515.9898508588834,
            s,
            fma(-20486.813161100774, eta, 7343.130973149263),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Intermediate_Phase_22_v2mRDv4` (IntPhaseVersion 105):
/// difference between the second intermediate collocation point and the
/// fourth ringdown collocation point.
pub(crate) fn phase_22_v2m_rdv4(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let no_spin =
        eta * fma(
            632.4731389009143,
            eta2,
            fma(101.21991715215253, eta, 0.9951733419499662),
        ) / fma(
            1.0,
            eta3,
            fma(
                1.8413983770369362,
                eta2,
                fma(0.11412314719189287, eta, 0.00016803066316882238),
            ),
        );
    let eq_spin = s * fma(
        eta3,
        fma(
            178.85850322465944,
            s2,
            fma(-14302.70838220423, s, 3531.552143264721),
        ),
        fma(
            0.3612417066833153,
            s2,
            fma(
                eta,
                fma(-147.53381808989846, s2, -697.6773920613674),
                fma(
                    4941.31613710257 * eta2,
                    s,
                    fma(16.89845522539974, s, 18.694178521101332),
                ),
            ),
        ),
    ) / fma(1.0, s2, fma(-2.7706595614504725, s, 2.965640445745779));
    let uneq_spin = dchi * delta * eta2 * fma(1693.326644293169 * eta2, s, 356.74395864902294);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Intermediate_Phase_22_v3mRDv4` (IntPhaseVersion 105):
/// difference between the third intermediate collocation point and the
/// fourth ringdown collocation point.
pub(crate) fn phase_22_v3m_rdv4(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta6 = eta3 * eta3;
    let s2 = s * s;
    let no_spin =
        eta * fma(
            -751.4184178636324,
            eta3,
            fma(
                688.3609087244353,
                eta2,
                fma(-227.46830225846668, eta, -5.126358906504587),
            ),
        ) / fma(
            1.0,
            eta2,
            fma(-0.7811680872741462, eta, -0.004551938711031158),
        );
    let eq_spin = s * fma(
        eta3,
        fma(
            224.77581754671272,
            s2,
            fma(1862.9097047992134, s, -936.3740515136005),
        ),
        fma(
            -2.84311102369862,
            s2,
            fma(
                eta,
                fma(-8.13494176717772, s2, 73.79645135116367),
                fma(
                    -(539.4071941841604 * eta2),
                    s,
                    fma(-0.9539250460041732, s, 0.1549280856660919),
                ),
            ),
        ),
    ) / fma(1.0, s, -1.5308507364054487);
    let uneq_spin = 2993.3598520496153 * dchi * delta * eta6;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Intermediate_Phase_22_d43` (IntPhaseVersion 105):
/// difference between the fourth and third intermediate collocation points.
pub(crate) fn phase_22_d43(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s2 = s * s;
    let no_spin = fma(
        670077.5414916876,
        eta4,
        fma(
            -967049.2793750163,
            eta3,
            fma(
                -282820.39946006844,
                eta2,
                fma(-906.746595921514, eta, 0.4248820426833804),
            ),
        ),
    ) / fma(19783.077247023448, eta2, fma(1670.9440812294847, eta, 1.0));
    let eq_spin = s * fma(
        eta4,
        fma(
            8745.201037897836,
            s2,
            fma(-227.7672020783411, s, -4939.642457025497),
        ),
        fma(
            eta2,
            fma(
                941.6974723887743,
                s2,
                fma(269.9234918621958, s, -711.8532052499075),
            ),
            fma(
                0.4818323187946999,
                s2,
                fma(
                    eta,
                    fma(-61.17248283184154, s2, 12.840649528989287),
                    fma(
                        eta3,
                        fma(
                            -4974.189172654984,
                            s2,
                            fma(-877.8811492839261, s, 3499.432393555856),
                        ),
                        fma(1.1366593671801855, s, 0.22814271667259703),
                    ),
                ),
            ),
        ),
    ) / fma(1.0, s, -1.2442293719740283);
    let uneq_spin = dchi * delta * fma(1493.3851099678195, eta, -514.8494071830514) * eta3;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Intermediate_Phase_22_Ansatz` evaluated operation by
/// operation as in LALSimulation (including its fused multiply-adds), for
/// the coefficient computations that must agree with LAL to the last bit.
pub(crate) fn phase_22_ansatz_exact(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    let fda = phase.f_damp;
    let dfr = mf - phase.f_ring;
    let lorentzian = (4.0 * phase.c_l_gr) / fma(4.0 * fda, fda, dfr * dfr);
    let mut phase_out = fma(phase.b1, powers.m_one, phase.b0);
    phase_out = fma(phase.b2, powers.m_two, phase_out);
    phase_out = fma(phase.b3, powers.m_three, phase_out);
    phase_out = fma(phase.b4, powers.m_four, phase_out);
    phase_out + lorentzian
}

/// `IMRPhenomX_Intermediate_Phase_22_AnsatzInt` evaluated operation by
/// operation as in LALSimulation, see [`phase_22_ansatz_exact`]; also
/// used for the phase at each frequency.
#[inline]
pub(crate) fn phase_22_ansatz_int_exact(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    let fda = phase.f_damp;
    let mut phase_out = fma(phase.b0, mf, phase.b1 * powers.log);
    phase_out = fma(-phase.b2, powers.m_one, phase_out);
    phase_out -= phase.b3 * powers.m_two / 2.0;
    phase_out -= phase.b4 * powers.m_three / 3.0;
    phase_out + (2.0 * phase.c_l_gr * ((mf - phase.f_ring) / (2.0 * fda)).atan()) / fda
}
