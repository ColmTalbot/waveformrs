//! Inspiral region of the (2,2) mode, `LALSimIMRPhenomX_inspiral.c`.
//!
//! Only the default versions are implemented: `InsAmpVersion = 103`
//! (canonical TaylorF2 amplitude with three pseudo-PN coefficients) and
//! `InsPhaseVersion = 104` (canonical 3.5PN TaylorF2 phase with four
//! pseudo-PN coefficients). Collocation fits use the effective spin
//! `chiPNHat`.

use std::f64::consts::PI;

use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct, POWERS_OF_LALPI};
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
    let no_spin = fma(
        0.4845148547154606,
        eta2,
        fma(-0.06098548699809163, eta, -0.015178276424448592),
    ) / fma(0.09799277215675059, eta, 1.0);
    let eq_spin = fma(
        fma(-0.14404522791467844, eta2, 0.01761591799745109),
        s2,
        fma(
            fma(0.10495263104245876, eta2, 0.02300153747158323),
            s,
            fma(-0.14189350657140673, eta, 0.04834642258922544) * eta * s3,
        ),
    ) / fma(-0.7340448493183307, s, 1.0);
    let uneq_spin = dchi * delta * eta4 * fma(34.90874132485147, eta, 0.0018724905795891192);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Amp_22_v3`: collocation point at 3/4 of the
/// inspiral amplitude region.
pub(crate) fn amp_22_v3(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let no_spin = fma(
        8.4630293045015,
        eta2,
        fma(-1.1970535595488723, eta, -0.058572000924124644),
    ) / fma(15.430818840453686, eta, 1.0);
    let eq_spin = fma(
        -0.018924013869130434,
        s,
        fma(
            eta2,
            fma(0.8282888482429105, s, 0.788717372588848),
            fma(
                eta,
                fma(-0.21291764491897636, s, -0.20646621646484237),
                -0.08746408292050666,
            ),
        ),
    ) * s
        / fma(1.0, s, -1.332123330797879);
    let uneq_spin = dchi * delta * eta4 * fma(105.84553997647659, eta, 0.004389995099201855);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Amp_22_v4`: collocation point at the end of the
/// inspiral amplitude region.
pub(crate) fn amp_22_v4(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s2 = s * s;
    let no_spin = fma(
        5.629598195000046,
        eta3,
        fma(
            -3.186012733446088,
            eta2,
            fma(1.617404703616985, eta, -0.16212854591357853),
        ),
    ) / fma(0.04507019231274476, eta, 1.0);
    let eq_spin = s * fma(
        eta3,
        fma(19.099315016873643, s, -41.0378120175805),
        fma(
            eta,
            fma(5.215501942120774, s, -4.127597118865669),
            fma(
                -0.31443470118113853,
                s,
                fma(
                    eta2,
                    fma(-18.80590889704093, s, 18.353433894421833),
                    1.0055835408962206,
                ),
            ),
        ),
    ) / fma(1.0, s2, fma(-5.717874483424523, s, 5.852706459485663));
    let uneq_spin = dchi * delta * eta4 * fma(208.92352600701068, eta, 0.05575955418803233);
    no_spin + eq_spin + uneq_spin
}

/// `pow_7_of` from `LALSimIMRPhenomXUtilities.c`.
pub(crate) fn pow_7_of(x: f64) -> f64 {
    let pow2 = x * x;
    pow2 * pow2 * pow2 * x
}

/// `IMRPhenomX_Inspiral_Amp_22_rho1`: pseudo-PN amplitude coefficient of
/// f^(7/3) from the collocation values `v1..v3` at `f1..f3`.
pub(crate) fn amp_22_rho1(v1: f64, v2: f64, v3: f64, f1: f64, f2: f64, f3: f64) -> f64 {
    let f1p1o3 = f1.cbrt();
    let f2p1o3 = f2.cbrt();
    let f3p1o3 = f3.cbrt();
    let f1p7o3 = pow_7_of(f1p1o3);
    let f1p8o3 = f1p7o3 * f1p1o3;
    let f2p7o3 = pow_7_of(f2p1o3);
    let f2p8o3 = f2p7o3 * f2p1o3;
    let f3p7o3 = pow_7_of(f3p1o3);
    let f3p8o3 = f3p7o3 * f3p1o3;
    fma(
        f1 * f1 * f1 * f2p8o3,
        v3,
        fma(
            -(f1p8o3 * (f2 * f2 * f2)),
            v3,
            fma(
                -(f1 * f1 * f1 * f3p8o3),
                v2,
                fma(
                    f1p8o3 * (f3 * f3 * f3),
                    v2,
                    fma(-(f2p8o3 * (f3 * f3 * f3)), v1, f2 * f2 * f2 * f3p8o3 * v1),
                ),
            ),
        ),
    ) / (f1p7o3 * (f1p1o3 - f2p1o3) * f2p7o3 * (f1p1o3 - f3p1o3) * (f2p1o3 - f3p1o3) * f3p7o3)
}

/// `IMRPhenomX_Inspiral_Amp_22_rho2`: pseudo-PN amplitude coefficient of
/// f^(8/3).
pub(crate) fn amp_22_rho2(v1: f64, v2: f64, v3: f64, f1: f64, f2: f64, f3: f64) -> f64 {
    let f1p1o3 = f1.cbrt();
    let f2p1o3 = f2.cbrt();
    let f3p1o3 = f3.cbrt();
    let f1p7o3 = pow_7_of(f1p1o3);
    let f2p7o3 = pow_7_of(f2p1o3);
    let f3p7o3 = pow_7_of(f3p1o3);
    fma(
        -(f1 * f1 * f1 * f2p7o3),
        v3,
        fma(
            f1p7o3 * (f2 * f2 * f2),
            v3,
            fma(
                f1 * f1 * f1 * f3p7o3,
                v2,
                fma(
                    -(f1p7o3 * (f3 * f3 * f3)),
                    v2,
                    fma(f2p7o3 * (f3 * f3 * f3), v1, -(f2 * f2 * f2 * f3p7o3 * v1)),
                ),
            ),
        ),
    ) / (f1p7o3 * (f1p1o3 - f2p1o3) * f2p7o3 * (f1p1o3 - f3p1o3) * (f2p1o3 - f3p1o3) * f3p7o3)
}

/// `IMRPhenomX_Inspiral_Amp_22_rho3`: pseudo-PN amplitude coefficient of
/// f^3.
pub(crate) fn amp_22_rho3(v1: f64, v2: f64, v3: f64, f1: f64, f2: f64, f3: f64) -> f64 {
    let f1p1o3 = f1.cbrt();
    let f2p1o3 = f2.cbrt();
    let f3p1o3 = f3.cbrt();
    let f1p7o3 = pow_7_of(f1p1o3);
    let f1p8o3 = f1p7o3 * f1p1o3;
    let f2p7o3 = pow_7_of(f2p1o3);
    let f2p8o3 = f2p7o3 * f2p1o3;
    let f3p7o3 = pow_7_of(f3p1o3);
    let f3p8o3 = f3p7o3 * f3p1o3;
    fma(
        -(f1p7o3 * f2p8o3),
        v3,
        fma(
            f1p8o3 * f2p7o3,
            v3,
            fma(
                f1p7o3 * f3p8o3,
                v2,
                fma(
                    -(f1p8o3 * f3p7o3),
                    v2,
                    fma(f2p8o3 * f3p7o3, v1, -(f2p7o3 * f3p8o3 * v1)),
                ),
            ),
        ),
    ) / (f1p7o3 * (f1p1o3 - f2p1o3) * f2p7o3 * (f1p1o3 - f3p1o3) * (f2p1o3 - f3p1o3) * f3p7o3)
}

/// `IMRPhenomX_Inspiral_Amp_22_rho1`, `_rho2` and `_rho3`: pseudo-PN
/// amplitude coefficients of f^(7/3), f^(8/3) and f^3 from the values `v`
/// of the collocation points at frequencies `f`.
pub(crate) fn amp_22_rho(v: [f64; 3], f: [f64; 3]) -> [f64; 3] {
    let [v1, v2, v3] = v;
    let [f1, f2, f3] = f;
    [
        amp_22_rho1(v1, v2, v3, f1, f2, f3),
        amp_22_rho2(v1, v2, v3, f1, f2, f3),
        amp_22_rho3(v1, v2, v3, f1, f2, f3),
    ]
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

/// `IMRPhenomX_Inspiral_Amp_22_Ansatz` evaluated operation by operation as
/// in LALSimulation (including its fused multiply-adds), for the
/// coefficient computations that must agree with LAL to the last bit.
pub(crate) fn amp_22_ansatz_exact(mf: f64, powers: &UsefulPowers, amp: &AmpCoefficients) -> f64 {
    let mut inner = fma(powers.one_third, amp.rho1, powers.two_thirds * amp.rho2);
    inner = fma(mf, amp.rho3, inner);
    let mut middle = fma(
        powers.one_third,
        amp.pn_four_thirds,
        powers.two_thirds * amp.pn_five_thirds,
    );
    middle = fma(mf, amp.pn_six_thirds, middle);
    middle = fma(mf, inner, middle);
    let mut pn_amp = fma(powers.one_third, amp.pn_one_third, amp.pn_initial);
    pn_amp = fma(powers.two_thirds, amp.pn_two_thirds, pn_amp);
    pn_amp = fma(mf, amp.pn_three_thirds, pn_amp);
    fma(mf, middle, pn_amp)
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
    let lalpi = PI;
    let delta = wf.delta;
    fma(
        3.0 * mf2,
        rho3,
        fma(
            8.0 / 3.0 * mf.powf(5.0 / 3.0),
            rho2,
            fma(
                7.0 / 3.0 * mf.powf(4.0 / 3.0),
                rho1,
                fma(
                    chi2l,
                    fma(-44.0, eta, fma(-81.0, delta, 81.0)),
                    chi1l * fma(81.0, 1.0 + delta, -(44.0 * eta)),
                ) * lalpi
                    / 48.0
                    + fma(1804.0, eta, -969.0) * lalpi.powf(2.0 / 3.0)
                        / (1008.0 * mf.powf(1.0 / 3.0))
                    + fma(
                        35371056.0,
                        eta2,
                        fma(
                            24.0 * fma(
                                857304.0,
                                chi2l2,
                                fma(
                                    -(994896.0 * chi1l),
                                    chi2l,
                                    fma(857304.0, chi1l2, -1975055.0),
                                ),
                            ),
                            eta,
                            fma(
                                -(10287648.0 * chi1l2),
                                1.0 + delta,
                                fma(
                                    10287648.0 * chi2l2,
                                    delta,
                                    fma(-10287648.0, chi2l2, -27312085.0),
                                ),
                            ),
                        ),
                    ) * lalpi.powf(4.0 / 3.0)
                        * mf.powf(1.0 / 3.0)
                        / 6.096384e6
                    + 5.0
                        * lalpi.powf(5.0 / 3.0)
                        * fma(
                            42840.0 * fma(4.0, eta, -1.0),
                            lalpi,
                            fma(
                                chi2l,
                                fma(
                                    -35632.0,
                                    eta2,
                                    fma(
                                        -fma(6048.0, chi2l2, 287213.0),
                                        -1.0 + delta,
                                        4.0 * fma(
                                            2083.0,
                                            delta,
                                            fma(1512.0 * chi2l2, -3.0 + delta, -93414.0),
                                        ) * eta,
                                    ),
                                ),
                                fma(
                                    -6048.0 * chi1l3,
                                    fma(3.0 + delta, eta, -1.0 - delta),
                                    chi1l
                                        * fma(
                                            -35632.0,
                                            eta2,
                                            fma(
                                                287213.0,
                                                1.0 + delta,
                                                -(4.0 * fma(2083.0, delta, 93414.0) * eta),
                                            ),
                                        ),
                                ),
                            ),
                        )
                        * mf.powf(2.0 / 3.0)
                        / 96768.0
                    - lalpi
                        * lalpi
                        * fma(
                            12.0 * eta,
                            fma(
                                21384760320.0,
                                lalpi * lalpi,
                                fma(
                                    -(8345272320.0 * (chi1l + chi2l)),
                                    lalpi,
                                    fma(
                                        77616.0 * chi1l2,
                                        fma(282718.0, delta, -158633.0),
                                        fma(
                                            -(21943440288.0 * chi2l2),
                                            delta,
                                            fma(
                                                -12312458928.0,
                                                chi2l2,
                                                fma(49568837472.0 * chi1l, chi2l, -545384828789.0),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                            fma(
                                7.0,
                                fma(
                                    12160253952.0
                                        * fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)),
                                    lalpi,
                                    fma(
                                        29362199328.0 * chi1l2,
                                        1.0 + delta,
                                        fma(
                                            -(29362199328.0 * chi2l2),
                                            delta,
                                            fma(29362199328.0, chi2l2, 177520268561.0),
                                        ),
                                    ),
                                ),
                                fma(
                                    -336.0
                                        * fma(
                                            1809550512.0,
                                            chi2l2,
                                            fma(
                                                -(2954929824.0 * chi1l),
                                                chi2l,
                                                fma(1809550512.0, chi1l2, -3248849057.0),
                                            ),
                                        ),
                                    eta2,
                                    -(324322727232.0 * eta3),
                                ),
                            ),
                        )
                        * mf
                        / 3.0042980352e10,
            ),
        ),
    )
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
    let no_spin = fma(
        8.913612508054944e6,
        eta4,
        fma(
            -3.9719980569125614e6,
            eta3,
            fma(
                376665.64637025696,
                eta2,
                fma(873401.6255736464, eta, 15415.000000000007),
            ),
        ),
    ) / fma(46.83697749859996, eta, 1.0);
    let eq_spin = s * fma(
        eta2,
        fma(
            895372.155565861,
            s2,
            fma(208349.45742548333, s, -1.9826323844247842e6),
        ),
        fma(
            -130668.37221912303,
            s2,
            fma(
                eta,
                fma(
                    -174952.69161683554,
                    s2,
                    fma(1.235279439281927e6, s, -1.0053073129700898e6),
                ),
                fma(
                    eta3,
                    fma(
                        -1.6894189124921719e6,
                        s2,
                        fma(-584728.050612325, s, 4.662143741417853e6),
                    ),
                    fma(-207180.42746987, s, 397951.95299014193),
                ),
            ),
        ),
    ) / fma(
        1.0,
        s3,
        fma(
            2.5298346636273306,
            s2,
            fma(3.5804521763363075, s, -9.675704197652225),
        ),
    );
    let uneq_spin = fma(
        -1296.9289110696955 * dchi2,
        eta,
        dchi * delta
            * eta
            * fma(
                47752.17032707405,
                s,
                fma(24703.28267342699, eta, -24708.109411857182),
            ),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Phase_22_d13`: difference v1 - v3 of the phase
/// collocation points (InsPhaseVersion 104).
pub(crate) fn phase_22_d13(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let s2 = s * s;
    let s3 = s2 * s;
    let s4 = s3 * s;
    let no_spin = fma(
        483033.0998073767,
        eta2,
        fma(-19943.076428555978, eta, -17294.000000000007),
    ) / fma(4.460294035404433, eta, 1.0);
    let eq_spin = s * fma(
        32614.091002011017,
        s4,
        fma(
            eta2,
            fma(
                1.0765877367729193e6,
                s3,
                fma(
                    -137991.37032619823,
                    s2,
                    fma(-758386.5685734496, s, -148368.4954044637),
                ),
            ),
            fma(
                19703.894135534803,
                s3,
                fma(
                    eta,
                    fma(
                        -303141.1945565486,
                        s3,
                        fma(
                            18298.307770807573,
                            s2,
                            fma(62190.404951852535, s, -58475.33302037833),
                        ),
                    ),
                    fma(
                        -2179.3505885609297,
                        s2,
                        fma(67663.42759836042, s, 68384.62786426462),
                    ),
                ),
            ),
        ),
    ) / fma(1.0, s, 2.0412979553629143);
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
    let no_spin = fma(
        -557253.0066989232,
        eta3,
        fma(
            1.1694356931282217e6,
            eta2,
            fma(-120297.86185566607, eta, -7579.300000000004),
        ),
    ) / fma(18.53018618227582, eta, 1.0);
    let eq_spin = s * fma(
        eta3,
        fma(
            847483.2999508682,
            s2,
            fma(-325423.99477314285, s, -593633.5370110178),
        ),
        fma(
            eta,
            fma(
                126736.43159783827,
                s2,
                fma(157036.45676788126, s, 50644.13475990821),
            ),
            fma(
                -44331.41741405198,
                s2,
                fma(
                    eta2,
                    fma(
                        -399712.22891153296,
                        s2,
                        fma(-50166.382087278434, s, 150022.21343386435),
                    ),
                    fma(-66228.9369155027, s, -27089.36915061857),
                ),
            ),
        ),
    ) / fma(
        1.0,
        s3,
        fma(
            -1.130185486082531,
            s2,
            fma(-3.062957826830017, s, -1.5232497464826662),
        ),
    );
    let uneq_spin = 3843.083992827935 * dchi * delta * eta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomX_Inspiral_Phase_22_d43`: difference v4 - v3 of the phase
/// collocation points (InsPhaseVersion 104).
pub(crate) fn phase_22_d43(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s2 = s * s;
    let no_spin = fma(
        28867.73328134167,
        eta2,
        fma(-31133.52170083207, eta, 2439.000000000001),
    ) / fma(0.41143032589262585, eta, 1.0);
    let eq_spin = s * fma(
        eta2,
        fma(117410.37432997991, s, 229284.04542668918),
        fma(
            9861.635308837876,
            s,
            fma(
                eta,
                fma(-25843.06175439942, s, -82355.86732027541),
                fma(
                    eta3,
                    fma(-386247.80765802023, s, -375818.0132734753),
                    16116.057657391262,
                ),
            ),
        ),
    ) / fma(1.0, s2, fma(0.25294420589064653, s, -3.7385208695213668));
    let uneq_spin = 194.5554531509207 * dchi * delta * eta;
    no_spin + eq_spin + uneq_spin
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

/// `IMRPhenomX_Inspiral_Phase_22_Ansatz` evaluated operation by operation as
/// in LALSimulation (including its fused multiply-adds), for the
/// coefficient computations that must agree with LAL to the last bit.
///
/// Terms with vanishing coefficients are skipped: `fma(0, x, y) == y`.
pub(crate) fn phase_22_ansatz_exact(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    let mut phase_in = phase.dphi0;
    phase_in = fma(phase.dphi2, powers.two_thirds, phase_in);
    phase_in = fma(phase.dphi3, mf, phase_in);
    phase_in = fma(phase.dphi4, powers.four_thirds, phase_in);
    phase_in = fma(phase.dphi5, powers.five_thirds, phase_in);
    phase_in = fma(phase.dphi6, powers.two, phase_in);
    phase_in = fma(phase.dphi6l * powers.two, powers.log, phase_in);
    phase_in = fma(phase.dphi7, powers.seven_thirds, phase_in);
    phase_in = fma(phase.dphi8, powers.eight_thirds, phase_in);
    phase_in = fma(phase.dphi8l * powers.eight_thirds, powers.log, phase_in);
    let mut pseudo = fma(phase.a0, powers.eight_thirds, phase.a1 * powers.three);
    pseudo = fma(phase.a2 * powers.eight_thirds, powers.two_thirds, pseudo);
    pseudo = fma(phase.a3 * powers.eight_thirds, powers.itself, pseudo);
    phase_in += pseudo;
    phase_in * powers.m_eight_thirds * phase.dphase_norm
}

/// `IMRPhenomX_Inspiral_Phase_22_AnsatzInt` evaluated operation by
/// operation as in LALSimulation, see [`phase_22_ansatz_exact`].
pub(crate) fn phase_22_ansatz_int_exact(
    mf: f64,
    powers: &UsefulPowers,
    phase: &PhaseCoefficients,
) -> f64 {
    let mut phasing = phase.phi0;
    phasing = fma(phase.phi2, powers.two_thirds, phasing);
    phasing = fma(phase.phi3, mf, phasing);
    phasing = fma(phase.phi4, powers.four_thirds, phasing);
    phasing = fma(phase.phi5l * powers.five_thirds, powers.log, phasing);
    phasing = fma(phase.phi6, powers.two, phasing);
    phasing = fma(phase.phi6l * powers.two, powers.log, phasing);
    phasing = fma(phase.phi7, powers.seven_thirds, phasing);
    phasing = fma(phase.phi8, powers.eight_thirds, phasing);
    phasing = fma(phase.phi8l * powers.eight_thirds, powers.log, phasing);
    let mut pseudo = fma(
        phase.sigma1,
        powers.eight_thirds,
        phase.sigma2 * powers.three,
    );
    pseudo = fma(phase.sigma3 * powers.one_third, powers.three, pseudo);
    pseudo = fma(phase.sigma4 * powers.two_thirds, powers.three, pseudo);
    phasing += pseudo;
    phasing * phase.phi_norm * powers.m_five_thirds
}
