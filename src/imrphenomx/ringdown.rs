//! Merger-ringdown region of the (2,2) mode, `LALSimIMRPhenomX_ringdown.c`.
//!
//! Only the default versions are implemented: `RDAmpVersion = 103`
//! (deformed Lorentzian with three coefficients) and `RDPhaseVersion = 105`
//! (five coefficient phase derivative ansatz). Collocation fits use the
//! effective spin `STotR`.

use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::UsefulPowers;
use crate::imrphenomx::phase::PhaseCoefficients;

/* ******************************** Amplitude ******************************** */

/// `IMRPhenomX_Ringdown_Amp_22_gamma2`: ringdown amplitude coefficient
/// gamma_2 (lambda in arXiv:2001.11412).
pub(crate) fn amp_22_gamma2(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let s2 = s * s;
    let eta2 = eta * eta;
    let no_spin = (0.8312293675316895 + 7.480371544268765 * eta - 18.256121237800397 * eta2)
        / (1. + 10.915453595496611 * eta - 30.578409433912874 * eta2);
    let eq_spin = (s
        * (0.5869408584532747
            + eta * (-0.1467158405070222 - 2.8489481072076472 * s)
            + 0.031852563636196894 * s
            + eta2 * (0.25295441250444334 + 4.6849496672664594 * s)))
        / (3.8775263105069953 - 3.41755361841226 * s + 1. * s2);
    let uneq_spin = -0.00548054788508203 * dchi * delta * eta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Amp_22_gamma3`: ringdown amplitude coefficient
/// gamma_3 (sigma in arXiv:2001.11412).
pub(crate) fn amp_22_gamma3(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let no_spin = (1.3666000000000007 - 4.091333144596439 * eta + 2.109081209912545 * eta2
        - 4.222259944408823 * eta3)
        / (1. - 2.7440263888207594 * eta);
    let eq_spin = (0.07179105336478316
        + eta2 * (2.331724812782498 - 0.6330998412809531 * s)
        + eta * (-0.8752427297525086 + 0.4168560229353532 * s)
        - 0.05633734476062242 * s)
        * s;
    let uneq_spin = 0. * delta * dchi;
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
    let no_spin = (0.03689164742964719 + 25.417967754401182 * eta + 162.52904393600332 * eta2)
        / (1. + 61.19874463331437 * eta - 29.628854485544874 * eta2);
    let eq_spin = (s
        * (-0.14352506969368556 + 0.026356911108320547 * s + 0.19967405175523437 * s2
            - 0.05292913111731128 * s3
            + eta3
                * (-48.31945248941757 - 3.751501972663298 * s
                    + 81.9290740950083 * s2
                    + 30.491948143930266 * s3
                    - 132.77982622925845 * s4)
            + eta
                * (-4.805034453745424 + 1.11147906765112 * s + 6.176053843938542 * s2
                    - 0.2874540719094058 * s3
                    - 8.990840289951514 * s4)
            - 0.18147275151697131 * s4
            + eta2
                * (27.675454081988036
                    - 2.398327419614959 * s
                    - 47.99096500250743 * s2
                    - 5.104257870393138 * s3
                    + 72.08174136362386 * s4)))
        / (-1.4160870461211452 + 1. * s);
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
    let numerator = dfr * dfr * g2 + 2.0 * fda * dfr * g3 + fda * fda * g2 * g3 * g3;
    let denominator = (dfr * dfr + dfd * dfd) * (dfr * dfr + dfd * dfd);
    prefactor * numerator / denominator
}

/// `IMRPhenomX_Ringdown_Amp_22_PeakFrequency`: frequency of the peak of
/// the ringdown amplitude, Eq. 5.14 of arXiv:2001.11412.
pub(crate) fn amp_22_peak_frequency(gamma2: f64, gamma3: f64, frd: f64, fda: f64) -> f64 {
    if gamma2 <= 1.0 {
        // If gamma2 > 1 the square root becomes imaginary and is dropped.
        (frd + fda * gamma3 * ((1.0 - gamma2 * gamma2).sqrt() - 1.0) / gamma2).abs()
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
    let no_spin = (-85.86062966719405 - 4616.740713893726 * eta - 4925.756920247186 * eta2
        + 7732.064464348168 * eta3
        + 12828.269960300782 * eta4
        - 39783.51698102803 * eta5)
        / (1. + 50.206318806624004 * eta);
    let eq_spin = (s
        * (33.335857451144356 - 36.49019206094966 * s
            + eta3 * (1497.3545918387515 - 101.72731770500685 * s) * s
            - 3.835967351280833 * s2
            + 2.302712009652155 * s3
            + eta2
                * (93.64156367505917 - 18.184492163348665 * s + 423.48863373726243 * s2
                    - 104.36120236420928 * s3
                    - 719.8775484010988 * s4)
            + 1.6533417657003922 * s4
            + eta
                * (-69.19412903018717 + 26.580344399838758 * s - 15.399770764623746 * s2
                    + 31.231253209893488 * s3
                    + 97.69027029734173 * s4)
            + eta4
                * (1075.8686153198323
                    - 3443.0233614187396 * s
                    - 4253.974688619423 * s2
                    - 608.2901586790335 * s3
                    + 5064.173605639933 * s4)))
        / (-1.3705601055555852 + 1. * s);
    let uneq_spin = dchi * delta * eta * (22.363215261437862 + 156.08206945239374 * eta);
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
    let no_spin = (eta * (0.7207992174994245 - 1.237332073800276 * eta + 6.086871214811216 * eta2))
        / (0.006851189888541745 + 0.06099184229137391 * eta - 0.15500218299268662 * eta2
            + 1. * eta3);
    let eq_spin = ((0.06519048552628343 - 25.25397971063995 * eta - 308.62513664956975 * eta4
        + 58.59408241189781 * eta2
        + 160.14971486043524 * eta3)
        * s
        + eta
            * (-5.215945111216946 + 153.95945758807616 * eta - 693.0504179144295 * eta2
                + 835.1725103648205 * eta3)
            * s2
        + (0.20035146870472367 - 0.28745205203100666 * eta - 47.56042058800358 * eta4) * s3
        + eta * (5.7756520242745735 - 43.97332874253772 * eta + 338.7263666984089 * eta3) * s4
        + (-0.2697933899920511 + 4.917070939324979 * eta
            - 22.384949087140086 * eta4
            - 11.61488280763592 * eta2)
            * s5)
        / (1. - 0.6628745847248266 * s);
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
    let no_spin = (eta * (-9.460253118496386 + 9.429314399633007 * eta + 64.69109972468395 * eta2))
        / (-0.0670554310666559 - 0.09987544893382533 * eta + 1. * eta2);
    let eq_spin = (17.36495157980372 * eta * s
        + eta3 * s * (930.3458437154668 + 808.457330742532 * s)
        + eta4 * s * (-774.3633787391745 - 2177.554979351284 * s - 1031.846477275069 * s2)
        + eta2 * s * (-191.00932194869588 - 62.997389062600035 * s + 64.42947340363101 * s2)
        + 0.04497628581617564 * s3)
        / (1. - 0.7267610313751913 * s);
    let uneq_spin = dchi * delta * (-36.66374091965371 + 91.60477826830407 * eta) * eta2;
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
    let no_spin =
        (eta * (-8.506898502692536 + 13.936621412517798 * eta)) / (-0.40919671232073945 + 1. * eta);
    let eq_spin = (eta
        * (1.7280582989361533 * s + 18.41570325463385 * s3 - 13.743271480938104 * s4)
        + eta2 * (73.8367329022058 * s - 95.57802408341716 * s3 + 215.78111099820157 * s4)
        + 0.046849371468156265 * s2
        + eta3
            * s
            * (-27.976989112929353 + 6.404060932334562 * s - 633.1966645925428 * s3
                + 109.04824706217418 * s2))
        / (1. - 0.6862449113932192 * s);
    let uneq_spin = 641.8965762829259 * dchi * delta * eta5;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Phase_22_d54`: difference v5 - v4 of the ringdown
/// phase collocation points.
pub(crate) fn phase_22_d54(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let no_spin = (eta * (7.05731400277692 + 22.455288821807095 * eta + 119.43820622871043 * eta2))
        / (0.26026709603623255 + 1. * eta);
    let eq_spin = (eta2 * (134.88158268621922 - 56.05992404859163 * s) * s
        + eta * s * (-7.9407123129681425 + 9.486783128047414 * s)
        + eta3 * s * (-316.26970506215554 + 90.31815139272628 * s))
        / (1. - 0.7162058321905909 * s);
    let uneq_spin = 43.82713604567481 * dchi * delta * eta3;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Ringdown_Phase_22_Ansatz` (RDPhaseVersion 105): phase
/// derivative in the merger-ringdown (to be scaled by 1/eta).
#[inline]
pub(crate) fn phase_22_ansatz(mf: f64, powers: &UsefulPowers, phase: &PhaseCoefficients) -> f64 {
    let dfr = mf - phase.f_ring;
    phase.c0
        + phase.c1 * powers.m_one_third
        + phase.c2 * powers.m_two
        + phase.c4 * powers.m_four
        + phase.c_l / (phase.f_damp2 + dfr * dfr)
}

/// `IMRPhenomX_Ringdown_Phase_22_AnsatzInt` (RDPhaseVersion 105): phase in
/// the merger-ringdown (to be scaled by 1/eta).
#[inline]
pub(crate) fn phase_22_ansatz_int(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    phase.c0 * mf + phase.c1_3ov2 * powers.two_thirds
        - phase.c2 * powers.m_one
        - phase.c4_ov3 * powers.m_three
        + phase.c_l_ov_fda * ((mf - phase.f_ring) * phase.inv_f_damp).atan()
}
