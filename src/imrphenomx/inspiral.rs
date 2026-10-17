//! Inspiral region of the (2,2) mode, `LALSimIMRPhenomX_inspiral.c`.
//!
//! Only the default versions are implemented: `InsAmpVersion = 103`
//! (canonical TaylorF2 amplitude with three pseudo-PN coefficients) and
//! `InsPhaseVersion = 104` (canonical 3.5PN TaylorF2 phase with four
//! pseudo-PN coefficients). Collocation fits use the effective spin
//! `chiPNHat`.

use std::f64::consts::PI;

use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::{UsefulPowers, WaveformStruct, POWERS_OF_LALPI};
use crate::imrphenomx::phase::PhaseCoefficients;

/* ******************************** Amplitude ******************************** */

/// `IMRPhenomX_Inspiral_Amp_22_v2`: collocation point at 2/4 of the
/// inspiral amplitude region.
pub(crate) fn amp_22_v2(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let no_spin = (-0.015178276424448592 - 0.06098548699809163 * eta + 0.4845148547154606 * eta2)
        / (1. + 0.09799277215675059 * eta);
    let eq_spin = ((0.02300153747158323 + 0.10495263104245876 * eta2) * s
        + (0.04834642258922544 - 0.14189350657140673 * eta) * eta * s3
        + (0.01761591799745109 - 0.14404522791467844 * eta2) * s2)
        / (1. - 0.7340448493183307 * s);
    let uneq_spin = dchi * delta * eta4 * (0.0018724905795891192 + 34.90874132485147 * eta);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Amp_22_v3`: collocation point at 3/4 of the
/// inspiral amplitude region.
pub(crate) fn amp_22_v3(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let no_spin = (-0.058572000924124644 - 1.1970535595488723 * eta + 8.4630293045015 * eta2)
        / (1. + 15.430818840453686 * eta);
    let eq_spin = ((-0.08746408292050666
        + eta * (-0.20646621646484237 - 0.21291764491897636 * s)
        + eta2 * (0.788717372588848 + 0.8282888482429105 * s)
        - 0.018924013869130434 * s)
        * s)
        / (-1.332123330797879 + 1. * s);
    let uneq_spin = dchi * delta * eta4 * (0.004389995099201855 + 105.84553997647659 * eta);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Amp_22_v4`: collocation point at the end of the
/// inspiral amplitude region.
pub(crate) fn amp_22_v4(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s2 = s * s;
    let no_spin = (-0.16212854591357853 + 1.617404703616985 * eta - 3.186012733446088 * eta2
        + 5.629598195000046 * eta3)
        / (1. + 0.04507019231274476 * eta);
    let eq_spin = (s
        * (1.0055835408962206 + eta2 * (18.353433894421833 - 18.80590889704093 * s)
            - 0.31443470118113853 * s
            + eta * (-4.127597118865669 + 5.215501942120774 * s)
            + eta3 * (-41.0378120175805 + 19.099315016873643 * s)))
        / (5.852706459485663 - 5.717874483424523 * s + 1. * s2);
    let uneq_spin = dchi * delta * eta4 * (0.05575955418803233 + 208.92352600701068 * eta);
    no_spin + eq_spin + uneq_spin
}

/// `pow_7_of` from `LALSimIMRPhenomXUtilities.c`.
fn pow_7_of(x: f64) -> f64 {
    let pow2 = x * x;
    pow2 * pow2 * pow2 * x
}

/// Powers f^(1/3), f^(7/3), f^(8/3) of the three collocation frequencies
/// as used by the pseudo-PN amplitude coefficient solutions.
struct RhoPowers {
    p7o3: [f64; 3],
    p8o3: [f64; 3],
    cubes: [f64; 3],
    denominator: f64,
}

impl RhoPowers {
    fn new(f: [f64; 3]) -> Self {
        let p1o3 = [f[0].cbrt(), f[1].cbrt(), f[2].cbrt()];
        let p7o3 = [pow_7_of(p1o3[0]), pow_7_of(p1o3[1]), pow_7_of(p1o3[2])];
        let p8o3 = [p7o3[0] * p1o3[0], p7o3[1] * p1o3[1], p7o3[2] * p1o3[2]];
        let cubes = [f[0] * f[0] * f[0], f[1] * f[1] * f[1], f[2] * f[2] * f[2]];
        let denominator = p7o3[0]
            * (p1o3[0] - p1o3[1])
            * p7o3[1]
            * (p1o3[0] - p1o3[2])
            * (p1o3[1] - p1o3[2])
            * p7o3[2];
        Self {
            p7o3,
            p8o3,
            cubes,
            denominator,
        }
    }
}

/// `IMRPhenomX_Inspiral_Amp_22_rho1`, `_rho2` and `_rho3`: pseudo-PN
/// amplitude coefficients of f^(7/3), f^(8/3) and f^3 from the values `v`
/// of the collocation points at frequencies `f`.
pub(crate) fn amp_22_rho(v: [f64; 3], f: [f64; 3]) -> [f64; 3] {
    let RhoPowers {
        p7o3,
        p8o3,
        cubes,
        denominator,
    } = RhoPowers::new(f);
    let [v1, v2, v3] = v;
    let [f1p7o3, f2p7o3, f3p7o3] = p7o3;
    let [f1p8o3, f2p8o3, f3p8o3] = p8o3;
    // F^3 written as F*F*F in LAL; `cubes` holds the same products.
    let [f1c, f2c, f3c] = cubes;
    let rho1 = (-(f2p8o3 * f3c * v1) + f2c * f3p8o3 * v1 + f1p8o3 * f3c * v2
        - f1c * f3p8o3 * v2
        - f1p8o3 * f2c * v3
        + f1c * f2p8o3 * v3)
        / denominator;
    let rho2 = (f2p7o3 * f3c * v1 - f2c * f3p7o3 * v1 - f1p7o3 * f3c * v2
        + f1c * f3p7o3 * v2
        + f1p7o3 * f2c * v3
        - f1c * f2p7o3 * v3)
        / denominator;
    let rho3 = (f2p8o3 * f3p7o3 * v1 - f2p7o3 * f3p8o3 * v1 - f1p8o3 * f3p7o3 * v2
        + f1p7o3 * f3p8o3 * v2
        + f1p8o3 * f2p7o3 * v3
        - f1p7o3 * f2p8o3 * v3)
        / denominator;
    [rho1, rho2, rho3]
}

/// `IMRPhenomX_Inspiral_Amp_22_Ansatz`: TaylorF2 PN amplitude plus
/// pseudo-PN terms (without the f^(-7/6) and normalisation prefactors).
#[inline]
pub(crate) fn amp_22_ansatz(mf: f64, powers: &UsefulPowers, amp: &AmpCoefficients) -> f64 {
    // The f^(1/3) coefficient vanishes and is omitted.
    amp.pn_initial
        + powers.two_thirds * amp.pn_two_thirds
        + mf * amp.pn_three_thirds
        + mf * (powers.one_third * amp.pn_four_thirds
            + powers.two_thirds * amp.pn_five_thirds
            + mf * amp.pn_six_thirds
            + mf * (powers.one_third * amp.rho1 + powers.two_thirds * amp.rho2 + mf * amp.rho3))
}

/// `IMRPhenomX_Inspiral_Amp_22_DAnsatz`: derivative of the inspiral
/// amplitude ansatz with respect to Mf.
pub(crate) fn amp_22_dansatz(mf: f64, wf: &WaveformStruct, amp: &AmpCoefficients) -> f64 {
    let eta = wf.eta;
    let chi1l = wf.chi1l;
    let chi2l = wf.chi2l;
    let rho1 = amp.rho1;
    let rho2 = amp.rho2;
    let rho3 = amp.rho3;

    let chi1l2 = chi1l * chi1l;
    let chi1l3 = chi1l2 * chi1l;
    let chi2l2 = chi2l * chi2l;

    let eta2 = eta * eta;
    let eta3 = eta * eta2;
    let mf2 = mf * mf;
    let pi = PI;
    let delta = wf.delta;

    ((chi2l * (81. - 81. * delta - 44. * eta) + chi1l * (81. * (1. + delta) - 44. * eta)) * pi)
        / 48.
        + ((-969. + 1804. * eta) * pi.powf(2. / 3.)) / (1008. * mf.powf(1. / 3.))
        + ((-27312085. - 10287648. * chi2l2 + 10287648. * chi2l2 * delta
            - 10287648. * chi1l2 * (1. + delta)
            + 24.
                * (-1975055. + 857304. * chi1l2 - 994896. * chi1l * chi2l + 857304. * chi2l2)
                * eta
            + 35371056. * eta2)
            * pi.powf(4. / 3.)
            * mf.powf(1. / 3.))
            / 6.096384e6
        + (5.
            * pi.powf(5. / 3.)
            * (-6048. * chi1l3 * (-1. - delta + (3. + delta) * eta)
                + chi1l
                    * (287213. * (1. + delta)
                        - 4. * (93414. + 2083. * delta) * eta
                        - 35632. * eta2)
                + chi2l
                    * (-((287213. + 6048. * chi2l2) * (-1. + delta))
                        + 4. * (-93414. + 1512. * chi2l2 * (-3. + delta) + 2083. * delta) * eta
                        - 35632. * eta2)
                + 42840. * (-1. + 4. * eta) * pi)
            * mf.powf(2. / 3.))
            / 96768.
        - (pi.powf(2.0)
            * (-336.
                * (-3248849057. + 1809550512. * chi1l2 - 2954929824. * chi1l * chi2l
                    + 1809550512. * chi2l2)
                * eta2
                - 324322727232. * eta3
                + 7. * (177520268561. + 29362199328. * chi2l2 - 29362199328. * chi2l2 * delta
                    + 29362199328. * chi1l2 * (1. + delta)
                    + 12160253952. * (chi1l + chi2l + chi1l * delta - chi2l * delta) * pi)
                + 12.
                    * eta
                    * (-545384828789. + 49568837472. * chi1l * chi2l
                        - 12312458928. * chi2l2
                        - 21943440288. * chi2l2 * delta
                        + 77616. * chi1l2 * (-158633. + 282718. * delta)
                        - 8345272320. * (chi1l + chi2l) * pi
                        + 21384760320. * pi.powf(2.0)))
            * mf)
            / 3.0042980352e10
        + (7.0 / 3.0) * mf.powf(4.0 / 3.0) * rho1
        + (8.0 / 3.0) * mf.powf(5.0 / 3.0) * rho2
        + 3. * mf2 * rho3
}

/* ********************************** Phase ********************************** */

/// `IMRPhenomX_Inspiral_Phase_22_v3`: value of the phase collocation point
/// v3 (InsPhaseVersion 104).
pub(crate) fn phase_22_v3(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let dchi2 = dchi * dchi;
    let no_spin = (15415.000000000007 + 873401.6255736464 * eta + 376665.64637025696 * eta2
        - 3.9719980569125614e6 * eta3
        + 8.913612508054944e6 * eta4)
        / (1. + 46.83697749859996 * eta);
    let eq_spin = (s
        * (397951.95299014193 - 207180.42746987 * s
            + eta3 * (4.662143741417853e6 - 584728.050612325 * s - 1.6894189124921719e6 * s2)
            + eta * (-1.0053073129700898e6 + 1.235279439281927e6 * s - 174952.69161683554 * s2)
            - 130668.37221912303 * s2
            + eta2 * (-1.9826323844247842e6 + 208349.45742548333 * s + 895372.155565861 * s2)))
        / (-9.675704197652225 + 3.5804521763363075 * s + 2.5298346636273306 * s2 + 1. * s3);
    let uneq_spin = -1296.9289110696955 * dchi2 * eta
        + dchi
            * delta
            * eta
            * (-24708.109411857182 + 24703.28267342699 * eta + 47752.17032707405 * s);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Phase_22_d13`: difference v1 - v3 of the phase
/// collocation points (InsPhaseVersion 104).
pub(crate) fn phase_22_d13(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let s4 = s3 * s;
    let no_spin = (-17294.000000000007 - 19943.076428555978 * eta + 483033.0998073767 * eta2)
        / (1. + 4.460294035404433 * eta);
    let eq_spin = (s
        * (68384.62786426462 + 67663.42759836042 * s - 2179.3505885609297 * s2
            + eta
                * (-58475.33302037833 + 62190.404951852535 * s + 18298.307770807573 * s2
                    - 303141.1945565486 * s3)
            + 19703.894135534803 * s3
            + eta2
                * (-148368.4954044637 - 758386.5685734496 * s - 137991.37032619823 * s2
                    + 1.0765877367729193e6 * s3)
            + 32614.091002011017 * s4))
        / (2.0412979553629143 + 1. * s);
    let uneq_spin = 12017.062595934838 * dchi * delta * eta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Phase_22_d23`: difference v2 - v3 of the phase
/// collocation points (InsPhaseVersion 104).
pub(crate) fn phase_22_d23(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let no_spin = (-7579.300000000004 - 120297.86185566607 * eta + 1.1694356931282217e6 * eta2
        - 557253.0066989232 * eta3)
        / (1. + 18.53018618227582 * eta);
    let eq_spin = (s
        * (-27089.36915061857 - 66228.9369155027 * s
            + eta2 * (150022.21343386435 - 50166.382087278434 * s - 399712.22891153296 * s2)
            - 44331.41741405198 * s2
            + eta * (50644.13475990821 + 157036.45676788126 * s + 126736.43159783827 * s2)
            + eta3 * (-593633.5370110178 - 325423.99477314285 * s + 847483.2999508682 * s2)))
        / (-1.5232497464826662 - 3.062957826830017 * s - 1.130185486082531 * s2 + 1. * s3);
    let uneq_spin = 3843.083992827935 * dchi * delta * eta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Phase_22_d43`: difference v4 - v3 of the phase
/// collocation points (InsPhaseVersion 104).
pub(crate) fn phase_22_d43(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let no_spin = (2439.000000000001 - 31133.52170083207 * eta + 28867.73328134167 * eta2)
        / (1. + 0.41143032589262585 * eta);
    let eq_spin = (s
        * (16116.057657391262
            + eta3 * (-375818.0132734753 - 386247.80765802023 * s)
            + eta * (-82355.86732027541 - 25843.06175439942 * s)
            + 9861.635308837876 * s
            + eta2 * (229284.04542668918 + 117410.37432997991 * s)))
        / (-3.7385208695213668 + 0.25294420589064653 * s + 1. * s2);
    let uneq_spin = 194.5554531509207 * dchi * delta * eta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Phase_22_Ansatz`: phase derivative in the inspiral,
/// TaylorF2 plus pseudo-PN terms (to be scaled by 1/eta).
///
/// Terms whose coefficients vanish for the default (GR, version 104) model
/// (`dphi_minus2`, `dphi_minus1`, `dphi1`, `dphi9`, `dphi9L`, `a4`) are omitted.
#[inline]
pub(crate) fn phase_22_ansatz(mf: f64, powers: &UsefulPowers, phase: &PhaseCoefficients) -> f64 {
    let mut phase_in = phase.dphi0; // f^{0/3}
    phase_in += phase.dphi2 * powers.two_thirds; // f^{2/3}
    phase_in += phase.dphi3 * mf; // f^{3/3}
    phase_in += phase.dphi4 * powers.four_thirds; // f^{4/3}
    phase_in += phase.dphi5 * powers.five_thirds; // f^{5/3}
    phase_in += phase.dphi6 * powers.two; // f^{6/3}
    phase_in += phase.dphi6l * powers.two * powers.log; // f^{6/3}, Log[f]
    phase_in += phase.dphi7 * powers.seven_thirds; // f^{7/3}
    phase_in += phase.dphi8 * powers.eight_thirds; // f^{8/3}
    phase_in += phase.dphi8l * powers.eight_thirds * powers.log; // f^{8/3}, Log[f]

    // Pseudo-PN coefficients
    phase_in += phase.a0 * powers.eight_thirds
        + phase.a1 * powers.three
        + phase.a2 * powers.eight_thirds * powers.two_thirds
        + phase.a3 * powers.eight_thirds * powers.itself;

    phase_in * powers.m_eight_thirds * phase.dphase_norm
}

/// `IMRPhenomX_Inspiral_Phase_22_AnsatzInt`: inspiral phase, TaylorF2 plus
/// pseudo-PN terms (to be scaled by 1/eta).
///
/// Terms whose coefficients vanish for the default (GR, version 104) model
/// (`phi_minus2`, `phi_minus1`, `phi1`, `phi5`, `phi9`, `phi9L`, `sigma5`)
/// are omitted.
#[inline]
pub(crate) fn phase_22_ansatz_int(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    let mut phasing = phase.phi0; // f^{-5/3}, Newtonian
    phasing += phase.phi2 * powers.two_thirds; // 1.0PN
    phasing += phase.phi3 * mf; // 1.5PN
    phasing += phase.phi4 * powers.four_thirds; // 2.0PN
    phasing += phase.phi5l * powers.five_thirds * powers.log; // 2.5PN, log
    phasing += phase.phi6 * powers.two; // 3.0PN
    phasing += phase.phi6l * powers.two * powers.log; // 3.0PN, log
    phasing += phase.phi7 * powers.seven_thirds; // 3.5PN
    phasing += phase.phi8 * powers.eight_thirds; // 4.0PN
    phasing += phase.phi8l * powers.eight_thirds * powers.log; // 4.0PN, log

    // Pseudo-PN coefficients
    phasing += phase.sigma1 * powers.eight_thirds
        + phase.sigma2 * powers.three
        + phase.sigma3 * powers.one_third * powers.three
        + phase.sigma4 * powers.two_thirds * powers.three;

    phasing * phase.phi_norm * powers.m_five_thirds
}

/// Leading-order normalisation of the inspiral phase derivative,
/// 5 / (128 pi^(5/3)) with pi^(5/3) as in `powers_of_lalpi`.
pub(crate) fn inspiral_dphase_norm() -> f64 {
    5.0 / (128.0 * POWERS_OF_LALPI.five_thirds)
}
