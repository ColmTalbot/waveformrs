//! Intermediate region of the (2,2) mode, `LALSimIMRPhenomX_intermediate.c`.
//!
//! Only the default versions are implemented: `IntAmpVersion = 104`
//! (fourth order polynomial for the inverse amplitude) and
//! `IntPhaseVersion = 105` (five coefficient phase derivative ansatz).
//! Collocation fits use the effective spin `STotR`.

use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::UsefulPowers;
use crate::imrphenomx::phase::PhaseCoefficients;

/* ******************************** Amplitude ******************************** */

/// `IMRPhenomX_Intermediate_Amp_22_vA`: value of the amplitude collocation
/// point at the centre of the intermediate region (IntAmpVersion 104).
pub(crate) fn amp_22_va(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let no_spin = (1.4873184918202145 + 1974.6112656679577 * eta + 27563.641024162127 * eta2
        - 19837.908020966777 * eta3)
        / (1. + 143.29004876335128 * eta + 458.4097306093354 * eta2);
    let eq_spin = (s
        * (27.952730865904343
            + eta * (-365.55631765202895 - 260.3494489873286 * s)
            + 3.2646808851249016 * s
            + 3011.446602208493 * eta2 * s
            - 19.38970173389662 * s2
            + eta3 * (1612.2681322644232 - 6962.675551371755 * s + 1486.4658089990298 * s2)))
        / (12.647425554323242 - 10.540154508599963 * s + 1. * s2);
    let uneq_spin = dchi * delta * (-0.016404056649860943 - 296.473359655246 * eta) * eta2;
    no_spin + eq_spin + uneq_spin
}

/// Powers and differences of the collocation frequencies shared by the
/// `IMRPhenomX_Intermediate_Amp_22_delta*` functions.
struct DeltaTerms {
    f12: f64,
    f13: f64,
    f14: f64,
    f15: f64,
    f22: f64,
    f23: f64,
    f24: f64,
    f42: f64,
    f43: f64,
    f44: f64,
    f45: f64,
    f1mf2: f64,
    f1mf4: f64,
    f2mf4: f64,
    f1mf22: f64,
    f2mf42: f64,
    f1mf43: f64,
}

impl DeltaTerms {
    fn new(f1: f64, f2: f64, f4: f64) -> Self {
        let f12 = f1 * f1;
        let f13 = f12 * f1;
        let f14 = f13 * f1;
        let f22 = f2 * f2;
        let f23 = f22 * f2;
        let f42 = f4 * f4;
        let f43 = f42 * f4;
        let f44 = f43 * f4;
        let f1mf2 = f1 - f2;
        let f1mf4 = f1 - f4;
        let f2mf4 = f2 - f4;
        Self {
            f12,
            f13,
            f14,
            f15: f14 * f1,
            f22,
            f23,
            f24: f23 * f2,
            f42,
            f43,
            f44,
            f45: f44 * f4,
            f1mf2,
            f1mf4,
            f2mf4,
            f1mf22: f1mf2 * f1mf2,
            f2mf42: f2mf4 * f2mf4,
            f1mf43: f1mf4 * f1mf4 * f1mf4,
        }
    }
}

/// `IMRPhenomX_Intermediate_Amp_22_delta0` to `_delta5` for IntAmpVersion
/// 104: coefficients of the fourth order polynomial for the inverse
/// amplitude given the boundary derivatives `d1`, `d4` and the inverse
/// amplitudes `v1`, `v2`, `v4` at `f1`, `f2`, `f4`.
///
/// The value `v3` at `f3` is unused for this version, and `delta5 = 0`.
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
    let DeltaTerms {
        f12,
        f13,
        f14,
        f15,
        f22,
        f23,
        f24,
        f42,
        f43,
        f44,
        f45,
        f1mf2,
        f1mf4,
        f2mf4,
        f1mf22,
        f2mf42,
        f1mf43,
    } = DeltaTerms::new(f1, f2, f4);
    let denominator = f1mf22 * f1mf43 * f2mf42;

    let delta0 = (-(d4 * f12 * f1mf22 * f1mf4 * f2 * f2mf4 * f4)
        + d1 * f1 * f1mf2 * f1mf4 * f2 * f2mf42 * f42
        + f42
            * (f2 * f2mf42 * (-4. * f12 + 3. * f1 * f2 + 2. * f1 * f4 - f2 * f4) * v1
                + f12 * f1mf43 * v2)
        + f12 * f1mf22 * f2 * (f1 * f2 - 2. * f1 * f4 - 3. * f2 * f4 + 4. * f42) * v4)
        / denominator;

    let delta1 = (d4 * f1 * f1mf22 * f1mf4 * f2mf4 * (2. * f2 * f4 + f1 * (f2 + f4))
        + f4 * (-(d1 * f1mf2 * f1mf4 * f2mf42 * (2. * f1 * f2 + (f1 + f2) * f4))
            - 2. * f1
                * (f44 * (v1 - v2)
                    + 3. * f24 * (v1 - v4)
                    + f14 * (v2 - v4)
                    + 4. * f23 * f4 * (-v1 + v4)
                    + 2. * f13 * f4 * (-v2 + v4)
                    + f1 * (2. * f43 * (-v1 + v2)
                        + 6. * f22 * f4 * (v1 - v4)
                        + 4. * f23 * (-v1 + v4)))))
        / denominator;

    let delta2 = (-(d4 * f1mf22 * f1mf4 * f2mf4 * (f12 + f2 * f4 + 2. * f1 * (f2 + f4)))
        + d1 * f1mf2 * f1mf4 * f2mf42 * (f1 * f2 + 2. * (f1 + f2) * f4 + f42)
        - 4. * f12 * f23 * v1
        + 3. * f1 * f24 * v1
        - 4. * f1 * f23 * f4 * v1
        + 3. * f24 * f4 * v1
        + 12. * f12 * f2 * f42 * v1
        - 4. * f23 * f42 * v1
        - 8. * f12 * f43 * v1
        + f1 * f44 * v1
        + f45 * v1
        + f15 * v2
        + f14 * f4 * v2
        - 8. * f13 * f42 * v2
        + 8. * f12 * f43 * v2
        - f1 * f44 * v2
        - f45 * v2
        - f1mf22
            * (f13
                + f2 * (3. * f2 - 4. * f4) * f4
                + f12 * (2. * f2 + f4)
                + f1 * (3. * f2 - 4. * f4) * (f2 + 2. * f4))
            * v4)
        / denominator;

    let delta3 = (d4 * f1mf22 * f1mf4 * f2mf4 * (2. * f1 + f2 + f4)
        - d1 * f1mf2 * f1mf4 * f2mf42 * (f1 + f2 + 2. * f4)
        + 2. * (f44 * (-v1 + v2)
            + 2. * f12 * f2mf42 * (v1 - v4)
            + 2. * f22 * f42 * (v1 - v4)
            + 2. * f13 * f4 * (v2 - v4)
            + f24 * (-v1 + v4)
            + f14 * (-v2 + v4)
            + 2. * f1 * f4 * (f42 * (v1 - v2) + f22 * (v1 - v4) + 2. * f2 * f4 * (-v1 + v4))))
        / denominator;

    let delta4 = (-(d4 * f1mf22 * f1mf4 * f2mf4) + d1 * f1mf2 * f1mf4 * f2mf42
        - 3. * f1 * f22 * v1
        + 2. * f23 * v1
        + 6. * f1 * f2 * f4 * v1
        - 3. * f22 * f4 * v1
        - 3. * f1 * f42 * v1
        + f43 * v1
        + f13 * v2
        - 3. * f12 * f4 * v2
        + 3. * f1 * f42 * v2
        - f43 * v2
        - f1mf22 * (f1 + 2. * f2 - 3. * f4) * v4)
        / denominator;

    [delta0, delta1, delta2, delta3, delta4, 0.0]
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
    let no_spin = (-82.54500000000004 - 5.58197349185435e6 * eta - 3.5225742421184325e8 * eta2
        + 1.4667258334378073e9 * eta3)
        / (1. + 66757.12830903867 * eta + 5.385164380400193e6 * eta2 + 2.5176585751772933e6 * eta3);
    let eq_spin = (s
        * (19.416719811164853 - 36.066611959079935 * s - 0.8612656616290079 * s2
            + eta2 * (170.97203068800542 - 107.41099349364234 * s - 647.8103976942541 * s3)
            + 5.95010003393006 * s3
            + eta3
                * (-1365.1499998427248
                    + 1152.425940764218 * s
                    + 415.7134909564443 * s2
                    + 1897.5444343138167 * s3
                    - 866.283566780576 * s4)
            + 4.984750041013893 * s4
            + eta
                * (207.69898051583655 - 132.88417400679026 * s - 17.671713040498304 * s2
                    + 29.071788188638315 * s3
                    + 37.462217031512786 * s4)))
        / (-1.1492259468169692 + 1. * s);
    let uneq_spin = dchi
        * delta
        * eta3
        * (7343.130973149263 - 20486.813161100774 * eta + 515.9898508588834 * s);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Intermediate_Phase_22_v2mRDv4` (IntPhaseVersion 105):
/// difference between the second intermediate collocation point and the
/// fourth ringdown collocation point.
pub(crate) fn phase_22_v2m_rdv4(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let no_spin = (eta
        * (0.9951733419499662 + 101.21991715215253 * eta + 632.4731389009143 * eta2))
        / (0.00016803066316882238
            + 0.11412314719189287 * eta
            + 1.8413983770369362 * eta2
            + 1. * eta3);
    let eq_spin = (s
        * (18.694178521101332
            + 16.89845522539974 * s
            + 4941.31613710257 * eta2 * s
            + eta * (-697.6773920613674 - 147.53381808989846 * s2)
            + 0.3612417066833153 * s2
            + eta3 * (3531.552143264721 - 14302.70838220423 * s + 178.85850322465944 * s2)))
        / (2.965640445745779 - 2.7706595614504725 * s + 1. * s2);
    let uneq_spin = dchi * delta * eta2 * (356.74395864902294 + 1693.326644293169 * eta2 * s);
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
    let no_spin = (eta
        * (-5.126358906504587 - 227.46830225846668 * eta + 688.3609087244353 * eta2
            - 751.4184178636324 * eta3))
        / (-0.004551938711031158 - 0.7811680872741462 * eta + 1. * eta2);
    let eq_spin = (s
        * (0.1549280856660919 - 0.9539250460041732 * s - 539.4071941841604 * eta2 * s
            + eta * (73.79645135116367 - 8.13494176717772 * s2)
            - 2.84311102369862 * s2
            + eta3 * (-936.3740515136005 + 1862.9097047992134 * s + 224.77581754671272 * s2)))
        / (-1.5308507364054487 + 1. * s);
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
    let no_spin = (0.4248820426833804
        - 906.746595921514 * eta
        - 282820.39946006844 * eta2
        - 967049.2793750163 * eta3
        + 670077.5414916876 * eta4)
        / (1. + 1670.9440812294847 * eta + 19783.077247023448 * eta2);
    let eq_spin = (s
        * (0.22814271667259703
            + 1.1366593671801855 * s
            + eta3 * (3499.432393555856 - 877.8811492839261 * s - 4974.189172654984 * s2)
            + eta * (12.840649528989287 - 61.17248283184154 * s2)
            + 0.4818323187946999 * s2
            + eta2 * (-711.8532052499075 + 269.9234918621958 * s + 941.6974723887743 * s2)
            + eta4 * (-4939.642457025497 - 227.7672020783411 * s + 8745.201037897836 * s2)))
        / (-1.2442293719740283 + 1. * s);
    let uneq_spin = dchi * delta * (-514.8494071830514 + 1493.3851099678195 * eta) * eta3;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Intermediate_Phase_22_Ansatz` (IntPhaseVersion 105): phase
/// derivative in the intermediate region (to be scaled by 1/eta).
#[inline]
pub(crate) fn phase_22_ansatz(mf: f64, powers: &UsefulPowers, phase: &PhaseCoefficients) -> f64 {
    let dfr = mf - phase.f_ring;
    // Lorentzian term with the GR coefficient cL = - a_RD dphase0.
    let lorentzian = phase.int_lorentzian_numerator / (phase.int_lorentzian_width2 + dfr * dfr);
    phase.b0
        + phase.b1 * powers.m_one
        + phase.b2 * powers.m_two
        + phase.b3 * powers.m_three
        + phase.b4 * powers.m_four
        + lorentzian
}

/// `IMRPhenomX_Intermediate_Phase_22_AnsatzInt` (IntPhaseVersion 105):
/// phase in the intermediate region (to be scaled by 1/eta).
#[inline]
pub(crate) fn phase_22_ansatz_int(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    phase.b0 * mf + phase.b1 * powers.log
        - phase.b2 * powers.m_one
        - phase.b3_ov2 * powers.m_two
        - phase.b4_ov3 * powers.m_three
        + phase.int_atan_prefactor * ((mf - phase.f_ring) * phase.inv_2_f_damp).atan()
}
