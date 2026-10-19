// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! Merger-ringdown region of the (2,2) mode, `LALSimIMRPhenomX_ringdown.c`.
//!
//! Only the default versions are implemented: `RDAmpVersion = 103`
//! (deformed Lorentzian with three coefficients) and `RDPhaseVersion = 105`
//! (five coefficient phase derivative ansatz). Collocation fits use the
//! effective spin `STotR`.

use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::{fma, UsefulPowers};
use crate::imrphenomx::phase::PhaseCoefficients;

/* ******************************** Amplitude ******************************** */

/// `IMRPhenomX_Ringdown_Amp_22_gamma2`: ringdown amplitude coefficient
/// gamma_2 (lambda in arXiv:2001.11412).
pub(crate) fn amp_22_gamma2(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let s2 = s * s;
    let eta2 = eta * eta;
    let no_spin = fma(
        -18.256121237800397,
        eta2,
        fma(7.480371544268765, eta, 0.8312293675316895),
    ) / fma(-30.578409433912874, eta2, fma(10.915453595496611, eta, 1.0));
    let eq_spin = s * fma(
        eta2,
        fma(4.6849496672664594, s, 0.25295441250444334),
        fma(
            0.031852563636196894,
            s,
            fma(
                eta,
                fma(-2.8489481072076472, s, -0.1467158405070222),
                0.5869408584532747,
            ),
        ),
    ) / fma(1.0, s2, fma(-3.41755361841226, s, 3.8775263105069953));
    let uneq_spin = -0.00548054788508203 * dchi * delta * eta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Amp_22_gamma3`: ringdown amplitude coefficient
/// gamma_3 (sigma in arXiv:2001.11412).
pub(crate) fn amp_22_gamma3(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let no_spin = fma(
        -4.222259944408823,
        eta3,
        fma(
            2.109081209912545,
            eta2,
            fma(-4.091333144596439, eta, 1.3666000000000007),
        ),
    ) / fma(-2.7440263888207594, eta, 1.0);
    let eq_spin = fma(
        -0.05633734476062242,
        s,
        fma(
            eta,
            fma(0.4168560229353532, s, -0.8752427297525086),
            fma(
                eta2,
                fma(-0.6330998412809531, s, 2.331724812782498),
                0.07179105336478316,
            ),
        ),
    ) * s;
    let uneq_spin = 0.0 * delta * dchi;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Amp_22_v1`: amplitude collocation point at the
/// peak frequency.
pub(crate) fn amp_22_v1(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let s3 = s * s2;
    let s4 = s * s3;
    let no_spin = fma(
        162.52904393600332,
        eta2,
        fma(25.417967754401182, eta, 0.03689164742964719),
    ) / fma(-29.628854485544874, eta2, fma(61.19874463331437, eta, 1.0));
    let eq_spin = s * fma(
        eta2,
        fma(
            72.08174136362386,
            s4,
            fma(
                -5.104257870393138,
                s3,
                fma(
                    -47.99096500250743,
                    s2,
                    fma(-2.398327419614959, s, 27.675454081988036),
                ),
            ),
        ),
        fma(
            -0.18147275151697131,
            s4,
            fma(
                eta,
                fma(
                    -8.990840289951514,
                    s4,
                    fma(
                        -0.2874540719094058,
                        s3,
                        fma(
                            6.176053843938542,
                            s2,
                            fma(1.11147906765112, s, -4.805034453745424),
                        ),
                    ),
                ),
                fma(
                    eta3,
                    fma(
                        -132.77982622925845,
                        s4,
                        fma(
                            30.491948143930266,
                            s3,
                            fma(
                                81.9290740950083,
                                s2,
                                fma(-3.751501972663298, s, -48.31945248941757),
                            ),
                        ),
                    ),
                    fma(
                        -0.05292913111731128,
                        s3,
                        fma(
                            0.19967405175523437,
                            s2,
                            fma(0.026356911108320547, s, -0.14352506969368556),
                        ),
                    ),
                ),
            ),
        ),
    ) / fma(1.0, s, -1.4160870461211452);
    let uneq_spin = -0.04426571511345366 * dchi * delta * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Amp_22_Ansatz`: deformed Lorentzian ringdown
/// amplitude (without the f^(-7/6) and normalisation prefactors).
#[inline]
pub(crate) fn amp_22_ansatz(mf: f64, amp: &AmpCoefficients) -> f64 {
    let dfr = mf - amp.f_ring;
    (-dfr * amp.gamma_r).exp() * amp.gamma_d13 / (dfr * dfr + amp.gamma_d2)
}

/// `XLALSimIMRPhenomXRingdownAmplitude22AnsatzAnalytical`: the ringdown
/// amplitude ansatz evaluated operation by operation as in LALSimulation
/// (including its fused multiply-add), from the coefficients `gamma1..3`.
pub(crate) fn amp_22_ansatz_exact(mf: f64, f_ring: f64, f_damp: f64, amp: &AmpCoefficients) -> f64 {
    let gamma_d13 = f_damp * amp.gamma1 * amp.gamma3;
    let gamma_r = amp.gamma2 / (amp.gamma3 * f_damp);
    let gamma_d2 = (f_damp * amp.gamma3) * (f_damp * amp.gamma3);
    let dfr = mf - f_ring;
    (-dfr * gamma_r).exp() * gamma_d13 / fma(dfr, dfr, gamma_d2)
}

/// `IMRPhenomX_Ringdown_Amp_22_DAnsatz`: derivative of the ringdown
/// amplitude ansatz.
pub(crate) fn amp_22_dansatz(mf: f64, f_ring: f64, f_damp: f64, amp: &AmpCoefficients) -> f64 {
    let g1 = amp.gamma1;
    let g2 = amp.gamma2;
    let g3 = amp.gamma3;
    let frd = f_ring;
    let fda = f_damp;
    let dfr = mf - frd;
    let dfd = fda * g3;
    let prefactor = -(-g2 * dfr / dfd).exp() * g1;
    let numerator = fma(
        fda * fda * g2 * g3,
        g3,
        fma(dfr * dfr, g2, 2.0 * fda * dfr * g3),
    );
    let denominator = fma(dfr, dfr, dfd * dfd) * fma(dfr, dfr, dfd * dfd);
    prefactor * numerator / denominator
}

/// `IMRPhenomX_Ringdown_Amp_22_PeakFrequency`: frequency of the peak of
/// the ringdown amplitude, Eq. 5.14 of arXiv:2001.11412.
pub(crate) fn amp_22_peak_frequency(gamma2: f64, gamma3: f64, frd: f64, fda: f64) -> f64 {
    if gamma2 <= 1.0 {
        // If gamma2 > 1 the square root becomes imaginary and is dropped.
        (frd + fda * gamma3 * (fma(-gamma2, gamma2, 1.0).sqrt() - 1.0) / gamma2).abs()
    } else {
        (frd + fda * (-1.0) * gamma3 / gamma2).abs()
    }
}

/* ********************************** Phase ********************************** */

/// `IMRPhenomX_Ringdown_Phase_22_v4`: ringdown phase collocation point at
/// the ringdown frequency.
pub(crate) fn phase_22_v4(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta4 * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let s4 = s3 * s;
    let no_spin = fma(
        -39783.51698102803,
        eta5,
        fma(
            12828.269960300782,
            eta4,
            fma(
                7732.064464348168,
                eta3,
                fma(
                    -4925.756920247186,
                    eta2,
                    fma(-4616.740713893726, eta, -85.86062966719405),
                ),
            ),
        ),
    ) / fma(50.206318806624004, eta, 1.0);
    let eq_spin = s * fma(
        eta4,
        fma(
            5064.173605639933,
            s4,
            fma(
                -608.2901586790335,
                s3,
                fma(
                    -4253.974688619423,
                    s2,
                    fma(-3443.0233614187396, s, 1075.8686153198323),
                ),
            ),
        ),
        fma(
            eta,
            fma(
                97.69027029734173,
                s4,
                fma(
                    31.231253209893488,
                    s3,
                    fma(
                        -15.399770764623746,
                        s2,
                        fma(26.580344399838758, s, -69.19412903018717),
                    ),
                ),
            ),
            fma(
                1.6533417657003922,
                s4,
                fma(
                    eta2,
                    fma(
                        -719.8775484010988,
                        s4,
                        fma(
                            -104.36120236420928,
                            s3,
                            fma(
                                423.48863373726243,
                                s2,
                                fma(-18.184492163348665, s, 93.64156367505917),
                            ),
                        ),
                    ),
                    fma(
                        2.302712009652155,
                        s3,
                        fma(
                            -3.835967351280833,
                            s2,
                            fma(
                                eta3 * fma(-101.72731770500685, s, 1497.3545918387515),
                                s,
                                fma(-36.49019206094966, s, 33.335857451144356),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    ) / fma(1.0, s, -1.3705601055555852);
    let uneq_spin = dchi * delta * eta * fma(156.08206945239374, eta, 22.363215261437862);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Phase_22_d12`: difference v1 - v2 of the ringdown
/// phase collocation points.
pub(crate) fn phase_22_d12(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let s4 = s3 * s;
    let s5 = s4 * s;
    let no_spin =
        eta * fma(
            6.086871214811216,
            eta2,
            fma(-1.237332073800276, eta, 0.7207992174994245),
        ) / fma(
            1.0,
            eta3,
            fma(
                -0.15500218299268662,
                eta2,
                fma(0.06099184229137391, eta, 0.006851189888541745),
            ),
        );
    let eq_spin = fma(
        fma(
            -11.61488280763592,
            eta2,
            fma(
                -22.384949087140086,
                eta4,
                fma(4.917070939324979, eta, -0.2697933899920511),
            ),
        ),
        s5,
        fma(
            eta * fma(
                338.7263666984089,
                eta3,
                fma(-43.97332874253772, eta, 5.7756520242745735),
            ),
            s4,
            fma(
                fma(
                    -47.56042058800358,
                    eta4,
                    fma(-0.28745205203100666, eta, 0.20035146870472367),
                ),
                s3,
                fma(
                    fma(
                        160.14971486043524,
                        eta3,
                        fma(
                            58.59408241189781,
                            eta2,
                            fma(
                                -308.62513664956975,
                                eta4,
                                fma(-25.25397971063995, eta, 0.06519048552628343),
                            ),
                        ),
                    ),
                    s,
                    eta * fma(
                        835.1725103648205,
                        eta3,
                        fma(
                            -693.0504179144295,
                            eta2,
                            fma(153.95945758807616, eta, -5.215945111216946),
                        ),
                    ) * s2,
                ),
            ),
        ),
    ) / fma(-0.6628745847248266, s, 1.0);
    let uneq_spin = -23.504907495268824 * dchi * delta * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Phase_22_d24`: difference v2 - v4 of the ringdown
/// phase collocation points.
pub(crate) fn phase_22_d24(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let no_spin =
        eta * fma(
            64.69109972468395,
            eta2,
            fma(9.429314399633007, eta, -9.460253118496386),
        ) / fma(
            1.0,
            eta2,
            fma(-0.09987544893382533, eta, -0.0670554310666559),
        );
    let eq_spin = fma(
        0.04497628581617564,
        s3,
        fma(
            eta2 * s,
            fma(
                64.42947340363101,
                s2,
                fma(-62.997389062600035, s, -191.00932194869588),
            ),
            fma(
                eta4 * s,
                fma(
                    -1031.846477275069,
                    s2,
                    fma(-2177.554979351284, s, -774.3633787391745),
                ),
                fma(
                    17.36495157980372 * eta,
                    s,
                    eta3 * s * fma(808.457330742532, s, 930.3458437154668),
                ),
            ),
        ),
    ) / fma(-0.7267610313751913, s, 1.0);
    let uneq_spin = dchi * delta * fma(91.60477826830407, eta, -36.66374091965371) * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Phase_22_d34`: difference v3 - v4 of the ringdown
/// phase collocation points.
pub(crate) fn phase_22_d34(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta5 = eta3 * eta2;
    let s2 = s * s;
    let s3 = s2 * s;
    let s4 = s3 * s;
    let no_spin = eta * fma(13.936621412517798, eta, -8.506898502692536)
        / fma(1.0, eta, -0.40919671232073945);
    let eq_spin = fma(
        eta3 * s,
        fma(
            109.04824706217418,
            s2,
            fma(
                -633.1966645925428,
                s3,
                fma(6.404060932334562, s, -27.976989112929353),
            ),
        ),
        fma(
            0.046849371468156265,
            s2,
            fma(
                eta,
                fma(
                    -13.743271480938104,
                    s4,
                    fma(1.7280582989361533, s, 18.41570325463385 * s3),
                ),
                eta2 * fma(
                    215.78111099820157,
                    s4,
                    fma(73.8367329022058, s, -(95.57802408341716 * s3)),
                ),
            ),
        ),
    ) / fma(-0.6862449113932192, s, 1.0);
    let uneq_spin = 641.8965762829259 * dchi * delta * eta5;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Phase_22_d54`: difference v5 - v4 of the ringdown
/// phase collocation points.
pub(crate) fn phase_22_d54(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let no_spin =
        eta * fma(
            119.43820622871043,
            eta2,
            fma(22.455288821807095, eta, 7.05731400277692),
        ) / fma(1.0, eta, 0.26026709603623255);
    let eq_spin = fma(
        eta3 * s,
        fma(90.31815139272628, s, -316.26970506215554),
        fma(
            eta2 * fma(-56.05992404859163, s, 134.88158268621922),
            s,
            eta * s * fma(9.486783128047414, s, -7.9407123129681425),
        ),
    ) / fma(-0.7162058321905909, s, 1.0);
    let uneq_spin = 43.82713604567481 * dchi * delta * eta3;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Phase_22_Ansatz` evaluated operation by operation
/// as in LALSimulation (including its fused multiply-adds), for the
/// coefficient computations that must agree with LAL to the last bit.
pub(crate) fn phase_22_ansatz_exact(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    let fda = phase.f_damp;
    let dfr = mf - phase.f_ring;
    let mut phase_rd = fma(phase.c1, powers.m_one_third, phase.c0);
    phase_rd = fma(phase.c2, powers.m_two, phase_rd);
    phase_rd = fma(phase.c4, powers.m_four, phase_rd);
    phase_rd + phase.c_l / fma(fda, fda, dfr * dfr)
}

/// `IMRPhenomX_Ringdown_Phase_22_AnsatzInt` evaluated operation by
/// operation as in LALSimulation, see [`phase_22_ansatz_exact`]; also
/// used for the phase at each frequency.
#[inline]
pub(crate) fn phase_22_ansatz_int_exact(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    let mut phase_rd = fma(phase.c0, mf, 1.5 * phase.c1 * powers.two_thirds);
    phase_rd = fma(-phase.c2, powers.m_one, phase_rd);
    phase_rd = fma(-phase.c4_ov3, powers.m_three, phase_rd);
    fma(
        phase.c_l_ov_fda,
        ((mf - phase.f_ring) / phase.f_damp).atan(),
        phase_rd,
    )
}
