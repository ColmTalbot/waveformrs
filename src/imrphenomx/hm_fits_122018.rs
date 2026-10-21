// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Parameter-space fits of the IMRPhenomXHM amplitude for the 122019
//! release (fits version 122018): the inspiral and intermediate collocation
//! points and the ringdown coefficients, from
//! `LALSimIMRPhenomXHM_inspiral.c`, `LALSimIMRPhenomXHM_intermediate.c` and
//! `LALSimIMRPhenomXHM_ringdown.c`.
//!
//! The fits are generated from the LAL sources with LAL's fused
//! multiply-adds (see [`fma`]); long sums are written as accumulations, which
//! does not change the operations.

use crate::imrphenomx::hm_internals::Fit;
use crate::imrphenomx::internals::{fma, WaveformStruct};

/// `IMRPhenomXHM_Insp_Amp_21_iv1` (version 122018): inspiral amplitude
/// collocation point at fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_21_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let s2 = s * s;
    let no_spin = fma(-4.0, eta, 1.0).sqrt()
        * fma(
            -98.3062466823662,
            eta5,
            fma(
                69.75910808095745,
                eta4,
                fma(
                    -16.706455229589558,
                    eta3,
                    fma(
                        1.963812986867654,
                        eta2,
                        fma(0.10740090317702103, eta, 0.037868557189995156),
                    ),
                ),
            ),
        );
    let eq_spin = fma(-4.0, eta, 1.0).sqrt()
        * s
        * fma(
            eta2,
            fma(
                -0.08217216660522082,
                s2,
                fma(-0.5985436493302649, s, -0.4735861262934258),
            ),
            fma(
                -0.008044970210401218,
                s,
                fma(0.10627108779259965, eta, -0.007963757232702219),
            ),
        );
    let uneq_spin = fma(
        -0.257787704938017 * wf.dchi * eta2,
        fma(8.75928187268504, eta2, 1.0),
        -(0.2597503605427412 * wf.dchi * eta2 * s),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_21_iv2` (version 122018): inspiral amplitude
/// collocation point at 0.75 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_21_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let s2 = s * s;
    let no_spin = fma(-4.0, eta, 1.0).sqrt()
        * fma(
            28.3307320191161,
            eta4,
            fma(
                -14.27268643447161,
                eta3,
                fma(
                    2.831411618302815,
                    eta2,
                    fma(-0.12579599745414977, eta, 0.05511628628738656),
                ),
            ),
        );
    let eq_spin = fma(-4.0, eta, 1.0).sqrt()
        * s
        * fma(
            -0.006383232900211555,
            s2,
            fma(
                eta2,
                fma(
                    -0.13099477081654226,
                    s2,
                    fma(-0.8321002363767336, s, -0.3409769288480959),
                ),
                fma(
                    -0.009520793625590234,
                    s,
                    fma(
                        eta,
                        fma(0.116470975986383, s, 0.09512553997347649),
                        -0.008692738851491525,
                    ),
                ),
            ),
        );
    let uneq_spin = fma(
        -0.2962753588645467 * wf.dchi * eta2,
        fma(1.3993978458830476, eta2, 1.0),
        -(0.17100612756133535 * wf.dchi * eta2 * s * fma(18.974303741922743 * eta2, delta, 1.0)),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_21_iv3` (version 122018): inspiral amplitude
/// collocation point at 0.5 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_21_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let s2 = s * s;
    let no_spin = fma(-4.0, eta, 1.0).sqrt()
        * fma(
            0.2428578654261086,
            eta2,
            fma(-0.0024538774422098405, eta, 0.059110044024271766),
        );
    let eq_spin = fma(-4.0, eta, 1.0).sqrt()
        * s
        * fma(
            -0.007685497720848461,
            s2,
            fma(
                eta2,
                fma(
                    0.004623467175906347,
                    s2,
                    fma(-0.12702579620537421, s, -0.016643018304732624),
                ),
                fma(-0.006952154764487417, s, -0.007044339356171243),
            ),
        );
    let uneq_spin = fma(
        -0.3172310538516028 * wf.dchi * fma(-2.9155919835488024, eta2, 1.0),
        eta2,
        -(0.11975485688200693 * wf.dchi * eta2 * s * fma(17.27626751837825 * eta2, delta, 1.0)),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_33_iv1` (version 122018): inspiral amplitude
/// collocation point at fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_33_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let eta2 = eta * eta;
    let s2 = s * s;
    let no_spin = fma(-4.0, eta, 1.0).sqrt()
        * fma(
            0.5584776628959615,
            eta2,
            fma(-0.14374841547279146, eta, -0.056586690934283326),
        )
        / (-0.3996185676368123 + eta);
    let eq_spin = fma(-4.0, eta, 1.0).sqrt()
        * s
        * fma(
            eta2,
            fma(
                5.343277927456605,
                s2,
                fma(-0.07646730296478217, s, -7.657635515668849),
            ),
            fma(
                fma(0.12482426029674777, s, 0.056042044149691175),
                s,
                eta * fma(-1.7827773156978863, s2, 2.1108074577110343),
            ),
        );
    let uneq_spin = 0.45866449225302536 * wf.dchi * fma(-9.603750707244906, eta2, 1.0) * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_33_iv2` (version 122018): inspiral amplitude
/// collocation point at 0.75 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_33_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let mut acc0 = 0.2137734510411439;
    acc0 = fma(-0.7692194209223682, eta, acc0);
    acc0 = fma(26.10570221351058, eta2, acc0);
    acc0 = fma(-316.0643979123107, eta3, acc0);
    acc0 = fma(2090.9063511488234, eta4, acc0);
    acc0 = fma(-6897.3285171507105, eta5, acc0);
    acc0 = fma(8968.893362362503, eta6, acc0);
    let no_spin = fma(-4.0, eta, 1.0).sqrt() * acc0;
    let eq_spin = fma(-4.0, eta, 1.0).sqrt()
        * s
        * fma(
            eta2,
            fma(
                8.351789152096782,
                s2,
                fma(1.6692205620001157, s, -6.51575737684721),
            ),
            fma(
                0.10701786057882816,
                s2,
                fma(
                    eta,
                    fma(
                        -2.110311135456494,
                        s2,
                        fma(-0.4683932646001618, s, 1.6484440612224325),
                    ),
                    fma(0.05924304311104228, s, 0.018546836505210842),
                ),
            ),
        );
    let uneq_spin = 0.3929315188124088 * wf.dchi * fma(-11.289452844364227, eta2, 1.0) * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_33_iv3` (version 122018): inspiral amplitude
/// collocation point at 0.5 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_33_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let mut acc0 = 0.2363760327127446;
    acc0 = fma(0.2855410252403732, eta, acc0);
    acc0 = fma(-10.159877125359897, eta2, acc0);
    acc0 = fma(162.65372389693505, eta3, acc0);
    acc0 = fma(-1154.7315106095564, eta4, acc0);
    acc0 = fma(3952.61320206691, eta5, acc0);
    acc0 = fma(-5207.67472857814, eta6, acc0);
    let no_spin = fma(-4.0, eta, 1.0).sqrt() * acc0;
    let eq_spin = fma(-4.0, eta, 1.0).sqrt()
        * s
        * fma(
            0.052963737236081304,
            s2,
            fma(
                eta,
                fma(
                    -0.2834348500565196,
                    s2,
                    fma(-0.1837289613228469, s, 0.15922377052827502),
                ),
                fma(0.048249943132325494, s, 0.04573095188775319),
            ),
        );
    let uneq_spin = 0.25187274502769835 * wf.dchi * fma(-12.172961866410864, eta2, 1.0) * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_32_iv1` (version 122018): inspiral amplitude
/// collocation point at fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_32_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let eta8 = eta.powf(8.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let mut acc0 = 0.019069933430190773;
    acc0 = fma(-0.19396651989685837, eta, acc0);
    acc0 = fma(11.95224600241255, eta2, acc0);
    acc0 = fma(-158.90113442757382, eta3, acc0);
    acc0 = fma(1046.65239329071, eta4, acc0);
    acc0 = fma(-3476.940285294999, eta5, acc0);
    acc0 = fma(4707.249209858949, eta6, acc0);
    let no_spin = fma(-3.0, eta, 1.0).sqrt() * acc0;
    let mut acc1 = 0.0046910348789512895;
    acc1 = fma(0.40231360805609434, eta, acc1);
    acc1 = fma(-0.0038263656140933152, s, acc1);
    acc1 = fma(0.018963579407636953, s2, acc1);
    acc1 = fma(
        eta2,
        fma(
            -0.9085620866763245,
            s3,
            fma(2.3753413452420133, s, -1.955352354930108),
        ),
        acc1,
    );
    acc1 = fma(0.02738043801805805, s3, acc1);
    acc1 = fma(
        eta3,
        fma(
            5.2255665027119145,
            s3,
            fma(
                0.49784942656123987,
                s2,
                fma(-7.9259853291789515, s, 7.977057990568723),
            ),
        ),
        acc1,
    );
    let eq_spin = fma(-3.0, eta, 1.0).sqrt() * s * acc1;
    let uneq_spin = fma(
        0.2735219358839411 * wf.dchi * eta2 * s,
        delta,
        fma(
            0.058560321425018165 * (wf.dchi * wf.dchi) * fma(-19.936477485971217, eta2, 1.0),
            eta2,
            1635.4240644598524 * wf.dchi * eta8 * delta,
        ),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_32_iv2` (version 122018): inspiral amplitude
/// collocation point at 0.75 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_32_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta8 = eta.powf(8.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(-3.0, eta, 1.0).sqrt()
        * fma(
            32.930430889650836,
            eta4,
            fma(
                -16.160563094841066,
                eta3,
                fma(
                    2.7200998230836158,
                    eta2,
                    fma(-0.09692699636236377, eta, 0.024621376891809633),
                ),
            ),
        );
    let mut acc0 = 0.008522695567479373;
    acc0 = fma(-1.1104639098529456, eta2, acc0);
    acc0 = fma(-0.00362963820787208, s, acc0);
    acc0 = fma(0.016978054142418417, s2, acc0);
    acc0 = fma(
        eta,
        fma(
            -0.1470288177047577,
            s3,
            fma(0.15878436411950506, s, 0.24280554040831698),
        ),
        acc0,
    );
    acc0 = fma(0.029465887557447824, s3, acc0);
    acc0 = fma(
        eta3,
        fma(
            2.5663386135613093,
            s3,
            fma(
                0.3381436950547799,
                s2,
                fma(-0.7550771176087877, s, 4.649438233164449),
            ),
        ),
        acc0,
    );
    let eq_spin = fma(-3.0, eta, 1.0).sqrt() * s * acc0;
    let uneq_spin = fma(
        0.1465612311350642 * wf.dchi * eta2 * s,
        delta,
        fma(
            -0.007061187955941243 * (wf.dchi * wf.dchi) * fma(-2.024701925508361, eta2, 1.0),
            eta2,
            215.06940561269835 * wf.dchi * eta8 * delta,
        ),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_32_iv3` (version 122018): inspiral amplitude
/// collocation point at 0.5 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_32_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta8 = eta.powf(8.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(-3.0, eta, 1.0).sqrt()
        * fma(
            -0.019043090896351363,
            eta3,
            fma(
                0.02620962593739105,
                eta2,
                fma(0.017454430190035, eta, -0.006150151041614737),
            ),
        )
        / (-0.2655505633361449 + eta);
    let eq_spin = fma(-3.0, eta, 1.0).sqrt()
        * s
        * fma(
            0.016323101357296865,
            s3,
            fma(
                eta2,
                fma(
                    0.6703688097531089,
                    s3,
                    fma(
                        0.04615964369925028,
                        s2,
                        fma(-0.008730335213434078, s, 0.5440307361144313),
                    ),
                ),
                fma(
                    0.012079324401547036,
                    s2,
                    fma(
                        eta * s,
                        fma(-0.15666140197050316, s2, 0.05592992411391443),
                        fma(0.00347699923233349, s, 0.011073381681404716),
                    ),
                ),
            ),
        );
    let uneq_spin = fma(
        0.10464021067936538 * wf.dchi * eta2 * s,
        delta,
        fma(
            -0.020140175824954427 * (wf.dchi * wf.dchi) * fma(-12.675522774051249, eta2, 1.0),
            eta2,
            -(417.3604094454253 * wf.dchi * eta8 * delta),
        ),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_44_iv1` (version 122018): inspiral amplitude
/// collocation point at fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_44_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let s2 = s * s;
    let no_spin = fma(-3.0, eta, 1.0).sqrt()
        * fma(
            35.21461767354364,
            eta4,
            fma(
                -15.988716302668415,
                eta3,
                fma(
                    1.9024723168424225,
                    eta2,
                    fma(0.1928897813606222, eta, 0.06190013067931406),
                ),
            ),
        );
    let eq_spin = fma(-3.0, eta, 1.0).sqrt()
        * s
        * fma(
            eta2,
            fma(
                1.1237234702682313,
                s2,
                fma(-0.6512813450832168, s, -2.1705970511116486),
            ),
            fma(
                0.07327810908370004,
                s2,
                fma(
                    eta,
                    fma(
                        -0.4482406547006759,
                        s2,
                        fma(0.12149520289658673, s, 0.6600413908621988),
                    ),
                    fma(0.044702230915643903, s, 0.011454874900772544),
                ),
            ),
        );
    let uneq_spin = fma(
        0.4766851579723911 * wf.dchi * fma(-15.950025762198988, eta2, 1.0),
        eta2,
        0.127900699645338 * (wf.dchi * wf.dchi) * fma(-15.79329306044842, eta2, 1.0) * eta2,
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_44_iv2` (version 122018): inspiral amplitude
/// collocation point at 0.75 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_44_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -1.2910560244510723,
        eta3,
        fma(
            0.2997223283799925,
            eta2,
            fma(-0.1469952725049322, eta, 0.08406011695496626),
        ),
    );
    let eq_spin = fma(
        fma(
            1.1090751649419717,
            eta2,
            fma(-0.4759412992529712, eta, 0.06028104440131858),
        ),
        s3,
        fma(
            fma(
                -1.1536009170220438,
                eta2,
                fma(0.26110236039648027, eta, 0.023924074703897662),
            ),
            s,
            fma(
                0.05736387085230215,
                eta2,
                fma(-0.1439868858871802, eta, 0.04479727299752669),
            ) * s2,
        ),
    );
    let uneq_spin = fma(
        0.10346324686812074 * (wf.dchi * wf.dchi) * fma(-16.135903382018213, eta2, 1.0),
        eta2,
        0.2648241309154185 * wf.dchi * eta2 * delta,
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Amp_44_iv3` (version 122018): inspiral amplitude
/// collocation point at 0.5 fcutInsp, rescaled by the leading order of the 22.
pub(crate) fn insp_amp_44_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.chi_pn_hat;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let s2 = s * s;
    let no_spin = fma(
        157.05954559045156,
        eta5,
        fma(
            -111.36325359782991,
            eta4,
            fma(
                28.404235115663706,
                eta3,
                fma(
                    -3.2466088293309885,
                    eta2,
                    fma(-0.025332770704783136, eta, 0.08212436946985402),
                ),
            ),
        ),
    );
    let eq_spin = s * fma(
        0.034468577523793176,
        s2,
        fma(
            eta,
            fma(
                -0.11199012099701576,
                s2,
                fma(-0.12754920943544915, s, -0.08968833480313292),
            ),
            fma(0.039491331923244756, s, 0.03488890057062679),
        ),
    );
    let uneq_spin = 0.2062291124580944 * wf.dchi * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_21_int1` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_21_int1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let s2 = s * s;
    let no_spin = fma(-4.0, eta2, eta).sqrt()
        * fma(
            30.14761650482866,
            eta2,
            fma(-25.594352690383847, eta, 21.256776327599113),
        );
    let eq_spin = fma(-4.0, eta2, eta).sqrt()
        * s
        * fma(
            eta,
            fma(
                -4.881427482271166,
                s2,
                fma(6.344825546437461, s, -1.1798437990445079),
            ),
            fma(-1.8167045597937677, s, -11.262044985632757),
        );
    let uneq_spin = fma(
        -3.6366100759176696 * (wf.dchi * wf.dchi) * fma(-4.0, eta, 1.0),
        eta,
        -(31.60048733143782 * wf.dchi * eta2 * fma(2.1502870640831855, eta2, 1.0)),
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_21_int2` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + 2 (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_21_int2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let no_spin = fma(-4.0, eta2, eta).sqrt()
        * fma(
            29.742565944285772,
            eta2,
            fma(-21.13596229438309, eta, 19.15445065708005),
        );
    let eq_spin = fma(-4.0, eta2, eta).sqrt()
        * s
        * fma(
            eta,
            fma(6.006571549661901, s, -2.913184982025043),
            fma(-2.123816950673979, s, -12.766814596085734),
        );
    let uneq_spin = -25.856046423804255 * wf.dchi * eta2 * fma(5.7871199275552, eta2, 1.0);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_21_int0` (version 122018): intermediate
/// amplitude collocation point of the first of two intermediate regions (extreme mass ratios).
pub(crate) fn inter_amp_21_int0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let no_spin = fma(
        518054.9485981792,
        eta3,
        fma(
            -24617.068739152524,
            eta2,
            fma(441.76285124642845, eta, 0.872895771366973),
        ),
    );
    let eq_spin = s * fma(
        eta2,
        fma(2640.919925429635, s, 3571.0393588230045),
        fma(
            -0.2687302438646897,
            s,
            fma(
                eta,
                fma(-113.29725582509889, s, -173.67847091983123),
                -0.0720494539485585,
            ),
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_21_dint0` (version 122018): intermediate
/// amplitude collocation point derivative at the collocation point of `int0`.
pub(crate) fn inter_amp_21_dint0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let no_spin = fma(
        -337621.44851304166,
        eta3,
        fma(
            13641.071903017495,
            eta2,
            fma(-93.1876950411214, eta, -0.8535048463050732),
        ),
    );
    let eq_spin = s * fma(
        eta,
        fma(189.63365484152857, s, 77.26350785961219),
        fma(
            -0.26539816223182355,
            s,
            fma(
                eta2,
                fma(-8172.057025783849, s, -1972.284151572111),
                -1.2067842398131878,
            ),
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_33_int1` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_33_int1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta6 = eta.powf(6.0);
    let no_spin = fma(-4.0, eta2, eta).sqrt()
        * fma(
            5621.897260910284,
            eta4,
            fma(
                -3744.785831952632,
                eta3,
                fma(
                    974.8550901501316,
                    eta2,
                    fma(-133.56611389260297, eta, 27.927652424857733),
                ),
            ),
        );
    let eq_spin = fma(-4.0, eta2, eta).sqrt()
        * s
        * fma(
            eta2,
            fma(83.33807316873204, s, 159.68630712802727),
            fma(
                5.059236579431119,
                s,
                fma(
                    eta,
                    fma(-37.07212326362276, s, -60.248696675045565),
                    7.348313807306079,
                ),
            ),
        );
    let uneq_spin = 1412.367880056888 * wf.dchi * eta6;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_33_int2` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + 2 (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_33_int2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta6 = eta.powf(6.0);
    let no_spin = fma(-4.0, eta2, eta).sqrt()
        * fma(
            53.04107631052987,
            eta2,
            fma(-18.666422946967764, eta, 20.162169689041903),
        );
    let eq_spin = fma(-4.0, eta2, eta).sqrt()
        * s
        * fma(
            eta2,
            fma(196.40535070402336, s, 91.31681057861915),
            fma(
                4.878506403725656,
                s,
                fma(
                    eta,
                    fma(-61.1244771077077, s, -33.707998325000965),
                    3.896260108714186,
                ),
            ),
        );
    let uneq_spin = 1637.4256048973248 * wf.dchi * eta6;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_33_int0` (version 122018): intermediate
/// amplitude collocation point of the first of two intermediate regions (extreme mass ratios).
pub(crate) fn inter_amp_33_int0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let no_spin = fma(
        743142.8286902909,
        eta3,
        fma(
            -34257.76380246282,
            eta2,
            fma(549.5183711492834, eta, 1.5852399637975103),
        ),
    );
    let eq_spin = s * fma(
        eta2,
        fma(6510.983041407191, s, 1185.1654325913717),
        fma(
            0.4253024979725725,
            s,
            fma(
                eta,
                fma(-174.5730646548662, s, -89.49451655594787),
                0.7436306553052219,
            ),
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_33_dint0` (version 122018): intermediate
/// amplitude collocation point derivative at the collocation point of `int0`.
pub(crate) fn inter_amp_33_dint0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let no_spin = fma(
        -310993.1309846956,
        eta3,
        fma(
            9262.994550540048,
            eta2,
            fma(101.4338937535679, eta, -4.691600252198376),
        ),
    );
    let eq_spin = s * fma(
        eta,
        fma(292.04566260701574, s, 734.7118618746317),
        fma(
            -0.40986595512314733,
            s,
            fma(
                eta2,
                fma(-5100.09336069277, s, -28714.904192060643),
                -4.198232394219111,
            ),
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_32_int1` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_32_int1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let mut acc0 = 6.523612598187996;
    acc0 = fma(-56.93956111746338, eta, acc0);
    acc0 = fma(1021.6414686597869, eta2, acc0);
    acc0 = fma(-12107.114370361525, eta3, acc0);
    acc0 = fma(76320.90587515048, eta4, acc0);
    acc0 = fma(-244144.92645448362, eta5, acc0);
    acc0 = fma(321790.55131499085, eta6, acc0);
    let no_spin = fma(-3.0, eta2, eta).sqrt() * acc0;
    let eq_spin = fma(-3.0, eta2, eta).sqrt()
        * s
        * fma(
            eta4,
            fma(7369.386205324284, s, -5200.659417644607),
            fma(
                eta2,
                fma(1341.2619677741804, s, 48.39743086535961),
                fma(
                    4.063724538613828,
                    s,
                    fma(
                        eta,
                        fma(-126.1263286618178, s, -37.87005271181108),
                        fma(
                            eta3,
                            fma(-5438.911035114849, s, 1790.8363334078751),
                            2.9649243713119895,
                        ),
                    ),
                ),
            ),
        );
    let uneq_spin = eta2
        * fma(
            -0.4386152975075188,
            fma(-(2.0 * wf.chi1l), wf.chi2l, wf.chi1l * wf.chi1l) + wf.chi2l * wf.chi2l,
            fma(
                wf.chi2l,
                fma(-7.324266404418883, s, 3.6527252109313233),
                wf.chi1l * fma(7.324266404418883, s, -3.6527252109313233),
            ) * delta,
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_32_int2` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + 2 (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_32_int2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let s2 = s * s;
    let no_spin = fma(-3.0, eta2, eta).sqrt()
        * fma(
            6214.038393898584,
            eta4,
            fma(
                -2113.9652334868965,
                eta3,
                fma(
                    271.105632998832,
                    eta2,
                    fma(-31.905244419036794, eta, 5.941845842405418),
                ),
            ),
        );
    let mut acc0 = -2.726472456645038;
    acc0 = fma(2.9454485454761827, s, acc0);
    acc0 = fma(
        eta3,
        fma(
            -11680.937129551317,
            s2,
            fma(-8474.190197512324, s, 10581.664858726683),
        ),
        acc0,
    );
    acc0 = fma(
        eta,
        fma(
            -145.5079981415436,
            s2,
            fma(-119.88112323140916, s, 98.08119212251981),
        ),
        acc0,
    );
    acc0 = fma(3.5684571473795095, s2, acc0);
    acc0 = fma(
        eta2,
        fma(
            2139.8290160628144,
            s2,
            fma(1686.7137359336039, s, -1595.8027347570667),
        ),
        acc0,
    );
    acc0 = fma(
        eta4,
        fma(
            20863.270079587106,
            s2,
            fma(13866.428366595079, s, -21488.25117198268),
        ),
        acc0,
    );
    let eq_spin = fma(-3.0, eta2, eta).sqrt() * s * acc0;
    let uneq_spin = 0.0038732029045487884 * wf.dchi * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_32_int0` (version 122018): intermediate
/// amplitude collocation point of the first of two intermediate regions (extreme mass ratios).
pub(crate) fn inter_amp_32_int0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        141004.29332593994,
        eta3,
        fma(
            -6626.167995915723,
            eta2,
            fma(115.81823862983131, eta, 0.24794156582503746),
        ),
    );
    let eq_spin = fma(
        fma(
            2118.213562899934,
            eta2,
            fma(-63.98011009365723, eta, 0.23428222791962147),
        ),
        s3,
        fma(
            fma(
                -1794.2301585086836,
                eta2,
                fma(35.10041265469983, eta, 0.21144389781375486),
            ),
            s,
            fma(
                1258.628375238807,
                eta2,
                fma(-37.038950686633, eta, 0.2781735549493081),
            ) * s2,
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_32_dint0` (version 122018): intermediate
/// amplitude collocation point derivative at the collocation point of `int0`.
pub(crate) fn inter_amp_32_dint0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let no_spin = fma(
        -95482.02951271653,
        eta3,
        fma(
            3694.1706648870427,
            eta2,
            fma(-14.604141885467747, eta, -0.3391808620221253),
        ),
    );
    let eq_spin = s * fma(
        eta,
        fma(188.2052849646003, s, 149.04865679660233),
        fma(
            -1.2356159239385598,
            s,
            fma(
                eta2,
                fma(-6332.389157828062, s, -5018.762853306415),
                -1.2844502090793946,
            ),
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_44_int1` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_44_int1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let no_spin = fma(-3.0, eta2, eta).sqrt()
        * fma(
            4779.096001663427,
            eta4,
            fma(
                -2612.999845214264,
                eta3,
                fma(
                    540.0541240482852,
                    eta2,
                    fma(-72.3834734399584, eta, 10.804555518381166),
                ),
            ),
        );
    let eq_spin = fma(-3.0, eta2, eta).sqrt()
        * s
        * fma(
            eta2,
            fma(106.91295627237112, s, 119.70401520575753),
            fma(
                3.0973959822174297,
                s,
                fma(
                    eta,
                    fma(-39.31284390368824, s, -47.94914754514519),
                    4.26336253142121,
                ),
            ),
        );
    let uneq_spin = fma(
        0.7262636326998003 * (wf.dchi * wf.dchi) * fma(-4.0, eta, 1.0),
        eta,
        3.001401833124412 * wf.dchi * eta2 * delta,
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_44_int2` (version 122018): intermediate
/// amplitude collocation point at fcutInsp + 2 (fcutRD - fcutInsp) / 3.
pub(crate) fn inter_amp_44_int2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let no_spin = fma(-3.0, eta2, eta).sqrt()
        * fma(
            6769.9274392345915,
            eta4,
            fma(
                -3194.0620894511508,
                eta3,
                fma(
                    508.07176447172264,
                    eta2,
                    fma(-53.221883492311235, eta, 9.020721305469884),
                ),
            ),
        );
    let eq_spin = fma(-3.0, eta2, eta).sqrt()
        * s
        * fma(
            eta2,
            fma(64.74544118094362, s, 96.41777041220982),
            fma(
                2.374434219852751,
                s,
                fma(
                    eta,
                    fma(-25.286684856422735, s, -38.38922554651356),
                    3.256591670091969,
                ),
            ),
        );
    let uneq_spin = 3.2337593375595417 * wf.dchi * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_44_int0` (version 122018): intermediate
/// amplitude collocation point of the first of two intermediate regions (extreme mass ratios).
pub(crate) fn inter_amp_44_int0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        249386.7511724409,
        eta3,
        fma(
            -11458.768824989507,
            eta2,
            fma(185.58965113823874, eta, 0.5664660641971224),
        ),
    );
    let eq_spin = fma(
        fma(
            -994.9241585715828,
            eta2,
            fma(-23.60535149579996, eta, 0.4423803798742513),
        ),
        s3,
        fma(
            fma(
                703.2622732011035,
                eta2,
                fma(-9.365114803167128, eta, 0.1741768776210781),
            ),
            s,
            fma(
                2833.5738711424974,
                eta2,
                fma(-62.13147149352512, eta, 0.20169229783048184),
            ) * s2,
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Amp_44_dint0` (version 122018): intermediate
/// amplitude collocation point derivative at the collocation point of `int0`.
pub(crate) fn inter_amp_44_dint0(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let no_spin = fma(
        -1728.7493675776548,
        eta2,
        fma(111.51170611049032, eta, -1.796444922382065),
    );
    let eq_spin = s * fma(
        eta,
        fma(128.64488103364167, s, 312.34531117524097),
        fma(
            -0.36655273031432567,
            s,
            fma(
                eta2,
                fma(-2927.019210835522, s, -11235.484645624338),
                -1.842119860613924,
            ),
        ),
    );
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_21_alambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_21_alambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let s2 = s * s;
    let no_spin = fma(-4.0, eta2, eta).sqrt()
        * fma(
            0.01042318959002753,
            eta2,
            fma(-0.0012619735607202085, eta, 0.00734983387668636),
        );
    let eq_spin = fma(-4.0, eta2, eta).sqrt()
        * s
        * fma(
            eta2,
            fma(
                0.020360826886107204,
                s2,
                fma(0.025956604949552205, s, -0.054621206928483663),
            ),
            fma(-0.0013927779195756036, s, -0.004839645742570202),
        );
    let uneq_spin =
        -0.018115657394753674 * wf.dchi * eta2 * fma(-(10.539795474715346 * eta2), delta, 1.0);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_21_lambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_21_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let no_spin = fma(
        1.8084545267208734,
        eta2,
        fma(0.12651770333481904, eta, 0.5566284518926176),
    );
    let eq_spin = fma(
        0.27420263462336675,
        s,
        fma(
            eta,
            fma(-0.31924748117518226, s, 0.059317243606471406),
            fma(
                eta2,
                fma(-3.4969956644617946, s, -2.101111399437034),
                0.29074922226651545,
            ),
        ),
    ) * s;
    let uneq_spin = 1.0122975748481835 * wf.dchi * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_21_sigma` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_21_sigma(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let no_spin = fma(0.0019318405961363861, eta, 1.2922261617161441);
    let eq_spin = fma(
        2.6625014134659772,
        eta2,
        fma(-0.6703778360948937, eta, 0.04927982551108649),
    ) * s;
    let uneq_spin = 1.2001101665670462
        * fma(
            -1.0 + wf.chi2l,
            wf.chi2l,
            fma(-(2.0 * wf.chi1l), wf.chi2l, wf.chi1l + wf.chi1l * wf.chi1l),
        )
        * eta2
        * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_33_alambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_33_alambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta4 = eta.powf(4.0);
    let no_spin = fma(-4.0, eta2, eta).sqrt()
        * fma(
            0.0898095508889557,
            eta2,
            fma(0.01202732427321774, eta, 0.013700854227665184),
        );
    let eq_spin = fma(-4.0, eta2, eta).sqrt()
        * fma(
            0.0035617441651710473,
            s,
            fma(
                eta,
                fma(-0.018186317026076343, s, -0.013132320758494439),
                0.0075858980586079065,
            ),
        )
        * s;
    let uneq_spin = eta4
        * fma(
            eta2 * fma(
                wf.chi1l,
                fma(-11.877399902871485, s, -4.2679864481479886),
                wf.chi2l * fma(11.877399902871485, s, 4.2679864481479886),
            ),
            delta,
            fma(
                wf.chi2l,
                fma(-0.05745949361626237, s, -0.09802218411554885),
                wf.chi1l * fma(0.05745949361626237, s, 0.09802218411554885),
            ),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_33_lambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_33_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let s2 = s * s;
    let no_spin = fma(
        1.471989765837694,
        eta2,
        fma(-0.06688558533374556, eta, 0.7435306475478924),
    );
    let eq_spin = s * fma(
        eta,
        fma(
            -0.1988991467974821,
            s2,
            fma(0.17261430318577403, s, -0.4809350398289311),
        ),
        fma(0.07564220573555203, s, 0.19457194111990656),
    );
    let uneq_spin = 1.8881959341735146 * wf.dchi * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_33_sigma` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_33_sigma(_wf: &WaveformStruct) -> f64 {
    let no_spin = 1.3;
    let eq_spin = 0.0;
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_32_alambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_32_alambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let mut acc0 = 0.00012587900257140724;
    acc0 = fma(0.03927886286971654, eta, acc0);
    acc0 = fma(-0.8109309606583066, eta2, acc0);
    acc0 = fma(8.820604907164254, eta3, acc0);
    acc0 = fma(-51.43344812454074, eta4, acc0);
    acc0 = fma(141.81940900657446, eta5, acc0);
    acc0 = fma(-140.0426973304466, eta6, acc0);
    let no_spin = acc0;
    let eq_spin = s * fma(
        1.0670798489407887 * eta3,
        s,
        fma(
            -0.00009363591928190229,
            s,
            fma(
                eta,
                fma(0.00568899459837252, s, -0.01083175709906557),
                fma(
                    eta2,
                    fma(-0.15569578955822236, s, 0.08311497969032984),
                    fma(
                        eta4,
                        fma(-2.09188976953315, s, -0.7849112300598181),
                        -0.00006001471234796344,
                    ),
                ),
            ),
        ),
    );
    let uneq_spin = fma(
        -0.04537308968659669 * (wf.dchi * wf.dchi) * eta2,
        fma(18.362371966229926, eta2, fma(-8.711096029480697, eta, 1.0)),
        wf.dchi
            * fma(
                9386.734883473799,
                eta3,
                fma(
                    -10001.774055779177,
                    eta2,
                    fma(3103.2516759087644, eta, -297.36978685672733),
                ),
            )
            * eta6,
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_32_lambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_32_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let no_spin = fma(-3.0, eta, 1.0).sqrt()
        * fma(
            0.7689553234961991,
            eta2,
            fma(-0.3197209728114808, eta, 0.0341611244787871),
        )
        / (fma(-0.43758296068790314, eta, 0.048429644168112324) + eta2);
    let eq_spin = fma(-3.0, eta, 1.0).sqrt()
        * s
        * fma(
            eta3,
            fma(175.47259563543165, s, -56.96407763839491),
            fma(
                9.790509295728649 * eta,
                s,
                fma(
                    eta2,
                    fma(-71.18182757443142, s, 25.536336676250748),
                    0.11057199932233873,
                ),
            ),
        );
    let uneq_spin = -5.002106168893265 * (wf.dchi * wf.dchi) * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_32_sigma` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_32_sigma(_wf: &WaveformStruct) -> f64 {
    let no_spin = 1.33;
    let eq_spin = 0.0;
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_44_alambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_44_alambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let s2 = s * s;
    let no_spin = fma(-3.0, eta2, eta).sqrt()
        * fma(
            233.23167528016833,
            eta5,
            fma(
                -136.10446022757958,
                eta4,
                fma(
                    28.298192768381554,
                    eta3,
                    fma(
                        -2.663803397359775,
                        eta2,
                        fma(0.09558474985614368, eta, 0.007904587819112173),
                    ),
                ),
            ),
        );
    let eq_spin = fma(-3.0, eta2, eta).sqrt()
        * s
        * fma(
            0.000657528128497184,
            s2,
            fma(
                eta3,
                fma(
                    0.19708289555425246,
                    s2,
                    fma(-0.9656098473922222, s, 0.9225758392294605),
                ),
                fma(
                    eta2,
                    fma(
                        -0.1389170612199015,
                        s2,
                        fma(0.4139026619690879, s, -0.2945455034809188),
                    ),
                    fma(
                        eta * fma(0.014107365722576927, s, -0.06166686913913691),
                        s,
                        fma(0.004122811292229324, s, 0.0049703757209330025),
                    ),
                ),
            ),
        );
    let uneq_spin = 0.00659873279539475 * wf.dchi * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_44_lambda` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_44_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let eta7 = eta.powf(7.0);
    let mut acc0 = 0.7702864948772887;
    acc0 = fma(32.81532373698395, eta, acc0);
    acc0 = fma(-1625.1795901450212, eta2, acc0);
    acc0 = fma(31305.458876573215, eta3, acc0);
    acc0 = fma(-297375.5347399236, eta4, acc0);
    acc0 = fma(1.4726521941846698e6, eta5, acc0);
    acc0 = fma(-3.616582470072637e6, eta6, acc0);
    acc0 = fma(3.4585865843680725e6, eta7, acc0);
    let no_spin = acc0;
    let eq_spin = fma(
        -(5.621635317494836 * eta3),
        s,
        fma(
            1.8738784391649446 * eta2,
            s,
            fma(-0.03011582009308575, s, 0.09034746027925727 * eta * s),
        ),
    ) / (-1.1340218677260014 + s);
    let uneq_spin = fma(
        0.959943270591552 * (wf.dchi * wf.dchi),
        eta2,
        0.853573071529436 * wf.dchi * eta2 * delta,
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Amp_44_sigma` (version 122018): ringdown amplitude
/// coefficient.
pub(crate) fn rd_amp_44_sigma(_wf: &WaveformStruct) -> f64 {
    let no_spin = 1.33;
    let eq_spin = 0.0;
    let uneq_spin = 0.0;
    no_spin + eq_spin + uneq_spin
}

/// The inspiral collocation points of mode `mode_int` (index in
/// `HIGHER_MODES`), at `fcutInsp`, `0.75 fcutInsp` and `0.5 fcutInsp`
/// (`pAmp->InspiralAmpFits`).
pub(crate) fn inspiral_fits(mode_int: usize) -> [Fit; 3] {
    const FITS: [[Fit; 3]; 4] = [
        [insp_amp_21_iv1, insp_amp_21_iv2, insp_amp_21_iv3],
        [insp_amp_33_iv1, insp_amp_33_iv2, insp_amp_33_iv3],
        [insp_amp_32_iv1, insp_amp_32_iv2, insp_amp_32_iv3],
        [insp_amp_44_iv1, insp_amp_44_iv2, insp_amp_44_iv3],
    ];
    FITS[mode_int]
}

/// The intermediate collocation points of mode `mode_int`: `int1`, `int2`
/// and, for the first of two intermediate regions, `int0` and `dint0`
/// (`pAmp->IntermediateAmpFits`).
pub(crate) fn intermediate_fits(mode_int: usize) -> [Fit; 4] {
    const FITS: [[Fit; 4]; 4] = [
        [
            inter_amp_21_int1,
            inter_amp_21_int2,
            inter_amp_21_int0,
            inter_amp_21_dint0,
        ],
        [
            inter_amp_33_int1,
            inter_amp_33_int2,
            inter_amp_33_int0,
            inter_amp_33_dint0,
        ],
        [
            inter_amp_32_int1,
            inter_amp_32_int2,
            inter_amp_32_int0,
            inter_amp_32_dint0,
        ],
        [
            inter_amp_44_int1,
            inter_amp_44_int2,
            inter_amp_44_int0,
            inter_amp_44_dint0,
        ],
    ];
    FITS[mode_int]
}

/// The ringdown coefficients `alambda`, `lambda` and `sigma` of mode
/// `mode_int` (`pAmp->RingdownAmpFits`).
pub(crate) fn ringdown_fits(mode_int: usize) -> [Fit; 3] {
    const FITS: [[Fit; 3]; 4] = [
        [rd_amp_21_alambda, rd_amp_21_lambda, rd_amp_21_sigma],
        [rd_amp_33_alambda, rd_amp_33_lambda, rd_amp_33_sigma],
        [rd_amp_32_alambda, rd_amp_32_lambda, rd_amp_32_sigma],
        [rd_amp_44_alambda, rd_amp_44_lambda, rd_amp_44_sigma],
    ];
    FITS[mode_int]
}
