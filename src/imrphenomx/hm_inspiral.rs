// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Inspiral region of the IMRPhenomXHM higher modes,
//! `LALSimIMRPhenomXHM_inspiral.c`: parameter-space fits of the amplitude
//! collocation points (version 122022) and of the extreme-mass-ratio phase
//! correction (version 122019), the PN amplitude coefficients, the
//! amplitude and phase ansatzes and the pseudo-PN amplitude coefficients.
//!
//! The fits are generated from the LAL sources with LAL's fused
//! multiply-adds (see [`fma`]); the fit functions take the 22-mode waveform
//! variables as LAL does.

use std::f64::consts::LN_2;

use num_complex::Complex;

use crate::imrphenomx::hm_internals::{
    Fit, HMAmpCoefficients, HMPhaseCoefficients, HMWaveformStruct,
};
use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct, POWERS_OF_LALPI};
use crate::imrphenomx::jet::Jet;

/// `IMRPhenomXHM_Insp_Amp_21_iv1` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_21_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    fma(
        delta
            * s1
            * fma(
                -16.75726972301224,
                fma(
                    876.4704270866478,
                    eta4,
                    fma(
                        -500.4821414428368,
                        eta3,
                        fma(
                            99.47071002831267,
                            eta2,
                            fma(-7.812073811917883, eta1, 1.1787350890261943),
                        ),
                    ),
                ),
                2.3439955698372663
                    * fma(
                        -4768.212511142035,
                        eta4,
                        fma(
                            2178.375177755584,
                            eta3,
                            fma(
                                -279.6409723479635,
                                eta2,
                                fma(7.176140122833879, eta1, 0.9373952326655807),
                            ),
                        ),
                    )
                    * s1,
            ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    347.1621871204769,
                    eta3,
                    fma(32.374226994179054, eta1, -(202.86279451816662 * eta2)),
                )
                * s1,
            sqroot,
            fma(
                chidiff1
                    * delta
                    * fma(
                        2638.3077554662905,
                        eta3,
                        fma(31.364509907424765, eta1, -(843.6414532232126 * eta2)),
                    ),
                sqroot,
                fma(
                    chidiff1 * eta5,
                    fma(
                        -134.98527058315528,
                        s2,
                        fma(987.635855365408, s1, -3962.5020052272976),
                    ),
                    delta
                        * fma(
                            220.77233521626252,
                            eta3,
                            fma(
                                -120.58166037019478,
                                eta2,
                                fma(16.6640319856064, eta1, 19.30531354642419),
                            ),
                        )
                        * sqroot,
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_Insp_Amp_21_iv2` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_21_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        delta
            * fma(
                chidiff1
                    * fma(
                        748.6659668096737,
                        eta3,
                        fma(36.88578491943111, eta1, -(321.2569602623214 * eta2)),
                    ),
                s1,
                chidiff1
                    * fma(
                        -3656.192371615541,
                        eta3,
                        fma(-95.42418611585117, eta1, 1217.338674959742 * eta2),
                    )
                    * s2,
            ),
        sqroot,
        fma(
            delta
                * s1
                * fma(
                    -10.809007068469844,
                    fma(
                        -16075.676528787526,
                        eta5,
                        fma(
                            11501.848278358668,
                            eta4,
                            fma(
                                -3064.7291458207815,
                                eta3,
                                fma(
                                    368.25127088896744,
                                    eta2,
                                    fma(-18.140470190766937, eta1, 1.1408749895922659),
                                ),
                            ),
                        ),
                    ),
                    1.0088254664333147
                        * fma(
                            -177907.92440833704,
                            eta5,
                            fma(
                                130764.22485304279,
                                eta4,
                                fma(
                                    -35561.24587952242,
                                    eta3,
                                    fma(
                                        4257.760834055382,
                                        eta2,
                                        fma(-192.2461213084741, eta1, 1.2322739396680107),
                                    ),
                                ),
                            ),
                        )
                        * s1,
                ),
            sqroot,
            fma(
                delta
                    * fma(
                        -1556.9475493549746,
                        eta4,
                        fma(
                            1046.3187137392433,
                            eta3,
                            fma(
                                -260.9979123967537,
                                eta2,
                                fma(30.883118528233638, eta1, 11.414934585404561),
                            ),
                        ),
                    ),
                sqroot,
                fma(
                    chidiff1 * eta5,
                    fma(
                        22.251142639924076,
                        s2,
                        fma(580.9465034962822, s1, -2898.9172078672705),
                    ),
                    delta
                        * fma(
                            chidiff2,
                            fma(
                                -417.5186332459383,
                                eta3,
                                fma(-18.541685007214625, eta1, 166.7427445020744 * eta2),
                            ),
                            chidiff1
                                * fma(
                                    2308.6520892707795,
                                    eta3,
                                    fma(41.61457952037761, eta1, -(779.9151607638761 * eta2)),
                                ),
                        )
                        * sqroot,
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_Insp_Amp_21_iv3` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_21_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    fma(
        delta
            * s1
            * fma(
                -(0.28384547935698246
                    * fma(
                        -24507.800226675925,
                        eta4,
                        fma(
                            14158.178124971111,
                            eta3,
                            fma(
                                -2707.9572896559484,
                                eta2,
                                fma(180.0378964169756, eta1, -0.8870770459576875),
                            ),
                        ),
                    )),
                s2,
                fma(
                    -7.1976397556450715,
                    fma(
                        1109.0598539312573,
                        eta4,
                        fma(
                            -635.1835857158857,
                            eta3,
                            fma(
                                120.87934060776237,
                                eta2,
                                fma(-6.587249958654174, eta1, 0.7603360145475428),
                            ),
                        ),
                    ),
                    -(0.0811847192323969
                        * fma(
                            -93272.17990295641,
                            eta4,
                            fma(
                                52586.63520999897,
                                eta3,
                                fma(
                                    -9548.970156895082,
                                    eta2,
                                    fma(517.4039644814231, eta1, 7.951454648295709),
                                ),
                            ),
                        )
                        * s1),
                ),
            ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    1034.4897198648962,
                    eta3,
                    fma(29.731007410186827, eta1, -(372.09609843131386 * eta2)),
                )
                * s1,
            sqroot,
            fma(
                delta
                    * fma(
                        -3914.3390484239667,
                        eta4,
                        fma(
                            2435.4098320959706,
                            eta3,
                            fma(
                                -568.0027716789259,
                                eta2,
                                fma(69.01761987509603, eta1, 5.390240326328237),
                            ),
                        ),
                    ),
                sqroot,
                fma(
                    chidiff1 * eta5,
                    fma(
                        16.379731479465033,
                        s2,
                        fma(157.94791186394787, s1, -2282.9983216879655),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            1350.476411541137,
                            eta3,
                            fma(21.935833431534224, eta1, -(460.7130131927895 * eta2)),
                        )
                        * sqroot,
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_Insp_Amp_33_iv1` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_33_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    fma(
        delta
            * s1
            * fma(
                -(0.26859293613553986
                    * fma(
                        639989.1197424605,
                        eta5,
                        fma(
                            -462437.2690007375,
                            eta4,
                            fma(
                                128973.72291098491,
                                eta3,
                                fma(
                                    -16664.276483466252,
                                    eta2,
                                    fma(888.7768796633982, eta1, -8.029069605349488),
                                ),
                            ),
                        ),
                    )),
                s2,
                fma(
                    -9.753610761811967,
                    fma(
                        -26407.40988450443,
                        eta5,
                        fma(
                            19794.855873435758,
                            eta4,
                            fma(
                                -5652.44754829634,
                                eta3,
                                fma(
                                    750.8933447725581,
                                    eta2,
                                    fma(-44.07982999150369, eta1, 1.7819678168496158),
                                ),
                            ),
                        ),
                    ),
                    0.014210376114848208
                        * fma(
                            7.911048294733888e6,
                            eta5,
                            fma(
                                -4.948175330328345e6,
                                eta4,
                                fma(
                                    1.1458022059130718e6,
                                    eta3,
                                    fma(
                                        -125763.47850622259,
                                        eta2,
                                        fma(7264.159472864562, eta1, -196.97328616330392),
                                    ),
                                ),
                            ),
                        )
                        * s1,
                ),
            ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    965.4670353331536,
                    eta3,
                    fma(29.221863857271703, eta1, -(348.1658322276406 * eta2)),
                )
                * s1,
            sqroot,
            fma(
                delta
                    * fma(
                        -8113.199260593833,
                        eta4,
                        fma(
                            4598.64854748635,
                            eta3,
                            fma(
                                -898.470298591732,
                                eta2,
                                fma(68.14577638017842, eta1, 20.94372147101354),
                            ),
                        ),
                    ),
                sqroot,
                fma(
                    chidiff1 * eta5,
                    fma(
                        1.4157230717300835,
                        s2,
                        fma(26.852777193715088, s1, 155.1434307076563),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            -141.85538063933927,
                            eta3,
                            fma(6.296698171560171, eta1, 15.81328761563562 * eta2),
                        )
                        * sqroot,
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_33_iv2` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_33_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    fma(
        delta
            * s1
            * fma(
                -(0.15528289788512326
                    * fma(
                        1.4516597366230418e6,
                        eta5,
                        fma(
                            -953010.5908118751,
                            eta4,
                            fma(
                                231384.03447562453,
                                eta3,
                                fma(
                                    -25091.14888164863,
                                    eta2,
                                    fma(1129.44485109116, eta1, -12.738184090548508),
                                ),
                            ),
                        ),
                    )),
                s2,
                fma(
                    -5.2580116732827085,
                    fma(
                        -36383.824902258915,
                        eta5,
                        fma(
                            25678.53964955809,
                            eta4,
                            fma(
                                -6879.681319382729,
                                eta3,
                                fma(
                                    861.1650630146937,
                                    eta2,
                                    fma(-48.20753331991333, eta1, 1.7794900975289085),
                                ),
                            ),
                        ),
                    ),
                    0.028627002336747746
                        * fma(
                            454655.37836367317,
                            eta5,
                            fma(
                                -168311.2370167227,
                                eta4,
                                fma(
                                    15062.821881048358,
                                    eta3,
                                    fma(
                                        -2287.0465658878725,
                                        eta2,
                                        fma(734.7581857539398, eta1, -50.57295946557892),
                                    ),
                                ),
                            ),
                        )
                        * s1,
                ),
            ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    327.2467193417936,
                    eta3,
                    fma(13.854211287141906, eta1, -(135.54916401086845 * eta2)),
                )
                * s1,
            sqroot,
            fma(
                delta
                    * fma(
                        -3100.2396041211605,
                        eta4,
                        fma(
                            1773.8883882162215,
                            eta3,
                            fma(
                                -357.6437296928763,
                                eta2,
                                fma(35.75290806646574, eta1, 12.501397517602173),
                            ),
                        ),
                    ),
                sqroot,
                fma(
                    chidiff1 * eta5,
                    fma(
                        -0.16889712161410445,
                        s2,
                        fma(37.141092711336846, s1, 161.62678370819597),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            -249.71072528701757,
                            eta3,
                            fma(3.4895829486899825, eta1, 51.07954458810889 * eta2),
                        )
                        * sqroot,
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_33_iv3` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_33_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let eta6 = eta1 * eta5;
    let eta7 = eta1 * eta6;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        delta
            * fma(
                -521279.50781871413,
                eta6,
                fma(
                    460912.76402476896,
                    eta5,
                    fma(
                        -162498.31493332976,
                        eta4,
                        fma(
                            29052.156454239772,
                            eta3,
                            fma(
                                -2774.5868652930294,
                                eta2,
                                fma(148.98231189649468, eta1, 5.148353856800232),
                            ),
                        ),
                    ),
                ),
            ),
        sqroot,
        fma(
            chidiff1 * eta7,
            fma(
                293.89657731755483,
                s2,
                fma(
                    793.5988192446642,
                    s1,
                    fma(47.94869769580204, chidiff2, 2667.3441342894776),
                ),
            ),
            fma(
                delta * s1,
                fma(
                    -1.2846335585108297,
                    fma(
                        -70.83896459998664,
                        eta4,
                        fma(
                            -34.842921754693386,
                            eta3,
                            fma(
                                23.26434219690476,
                                eta2,
                                fma(1.37856806162599, eta1, 0.09991079016763821),
                            ),
                        ),
                    ),
                    -(0.03496714763391888
                        * fma(
                            -43221.0453556626,
                            eta4,
                            fma(
                                22714.70643022915,
                                eta3,
                                fma(
                                    -3736.1574640444287,
                                    eta2,
                                    fma(188.38585449575902, eta1, -0.230558571912664),
                                ),
                            ),
                        )
                        * s1),
                ),
                fma(
                    chidiff1 * delta,
                    fma(
                        -110.10276573567405,
                        eta3,
                        fma(-0.5869777957488564, eta1, 32.65536124256588 * eta2),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            113.77466499598913,
                            eta3,
                            fma(3.524800489907584, eta1, -(40.26479860265549 * eta2)),
                        )
                        * s1,
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_32_iv1` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_32_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let eta6 = eta1 * eta5;
    let s1 = s;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1 * fma(
            0.6343034651912586,
            fma(
                30602.210950069973,
                eta5,
                fma(
                    -24776.529123239357,
                    eta4,
                    fma(
                        7616.234910399297,
                        eta3,
                        fma(
                            -1087.6241783616488,
                            eta2,
                            fma(78.98200041834092, eta1, -2.5844888818001737),
                        ),
                    ),
                ),
            ),
            -(0.062088720220899465
                * fma(
                    198366.10615196894,
                    eta5,
                    fma(
                        -138572.93700180828,
                        eta4,
                        fma(
                            33822.437731298516,
                            eta3,
                            fma(
                                -3124.4712274775407,
                                eta2,
                                fma(36.01386705325694, eta1, 6.5586380356588565),
                            ),
                        ),
                    ),
                )
                * s1),
        ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    -2069.5349737049073,
                    eta4,
                    fma(
                        1160.9833883184604,
                        eta3,
                        fma(12.488240781993923, eta1, -(209.32038774208385 * eta2)),
                    ),
                )
                * s1,
            sqroot,
            fma(
                fma(
                    chidiff1 * delta,
                    fma(
                        -566.6420939162068,
                        eta4,
                        fma(
                            278.9717709112207,
                            eta3,
                            fma(-0.739317114582042, eta1, -(47.473246070362634 * eta2)),
                        ),
                    ),
                    chidiff2
                        * fma(
                            23.783684827838247,
                            eta4,
                            fma(
                                -24.37776782232888,
                                eta3,
                                fma(-0.5873680378268906, eta1, 6.692187014925888 * eta2),
                            ),
                        ),
                ),
                sqroot,
                fma(
                    -55436.68857586653,
                    eta6,
                    fma(
                        51904.91007211022,
                        eta5,
                        fma(
                            -19279.186145681153,
                            eta4,
                            fma(
                                3585.9269057886418,
                                eta3,
                                fma(
                                    -343.3143244815765,
                                    eta2,
                                    fma(4.94285331708559, eta1, 3.2940434453819694),
                                ),
                            ),
                        ),
                    ),
                ) * sqroot,
            ),
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_32_iv2` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_32_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let eta6 = eta1 * eta5;
    let s1 = s;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1 * fma(
            0.7718201508695763,
            fma(
                -970.2126139442341,
                eta4,
                fma(
                    712.9101229418721,
                    eta3,
                    fma(
                        -186.5001124789369,
                        eta2,
                        fma(26.432880113146012, eta1, -1.3012906461000349),
                    ),
                ),
            ),
            0.04832734931068797
                * fma(
                    6858.003573909322,
                    eta4,
                    fma(
                        -2537.636347529708,
                        eta3,
                        fma(
                            1.8360177574514709,
                            eta2,
                            fma(78.98681284391004, eta1, -5.9999628512498315),
                        ),
                    ),
                )
                * s1,
        ),
        sqroot,
        fma(
            fma(
                chidiff2,
                fma(
                    -8.78587250202154,
                    eta3,
                    fma(-0.03940151060321499, eta1, 1.9034209537174116 * eta2),
                ),
                chidiff1
                    * delta
                    * fma(
                        0.36790005839460627,
                        eta3,
                        fma(-1.704299788495861, eta1, -(4.923510922214181 * eta2)),
                    ),
            ),
            sqroot,
            fma(
                30950.67053883646,
                eta6,
                fma(
                    -25472.73682092519,
                    eta5,
                    fma(
                        7970.929280907627,
                        eta4,
                        fma(
                            -1139.110227749627,
                            eta3,
                            fma(
                                60.10368251688146,
                                eta2,
                                fma(-5.1846950040514335, eta1, 2.2911849711339123),
                            ),
                        ),
                    ),
                ),
            ) * sqroot,
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_32_iv3` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_32_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let eta6 = eta1 * eta5;
    let s1 = s;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1 * fma(
            0.785171689871352,
            fma(
                -1115.008825153004,
                eta4,
                fma(
                    744.7699462372949,
                    eta3,
                    fma(
                        -167.9231114864569,
                        eta2,
                        fma(18.70856733065619, eta1, -0.4634745514643032),
                    ),
                ),
            ),
            0.13469300326662165
                * fma(
                    -472.02962341590774,
                    eta4,
                    fma(
                        1136.8367114738041,
                        eta3,
                        fma(
                            -483.7040402103785,
                            eta2,
                            fma(72.17373498208947, eta1, -2.7311391326835133),
                        ),
                    ),
                )
                * s1,
        ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    80.67120303906654,
                    eta3,
                    fma(1.904350804857431, eta1, -(25.565242391371093 * eta2)),
                )
                * s1,
            sqroot,
            fma(
                fma(
                    chidiff2,
                    fma(
                        -14.078156877577733,
                        eta3,
                        fma(-0.6358511175987503, eta1, 5.555088747533164 * eta2),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            36.20685853857613,
                            eta3,
                            fma(0.23205448591711159, eta1, -(19.46049432345157 * eta2)),
                        ),
                ),
                sqroot,
                fma(
                    38785.53683146884,
                    eta6,
                    fma(
                        -30849.228936891763,
                        eta5,
                        fma(
                            8880.536804741967,
                            eta4,
                            fma(
                                -976.9940395257111,
                                eta3,
                                fma(
                                    -17.51265776660515,
                                    eta2,
                                    fma(7.380126197972549, eta1, 1.1525594672495008),
                                ),
                            ),
                        ),
                    ),
                ) * sqroot,
            ),
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_44_iv1` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_44_iv1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let eta6 = eta1 * eta5;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1 * fma(
            -(0.023693246676754775
                * fma(
                    1.883460394974573e7,
                    eta6,
                    fma(
                        -1.7434160932177313e7,
                        eta5,
                        fma(
                            6.255666490084672e6,
                            eta4,
                            fma(
                                -1.0803340811860575e6,
                                eta3,
                                fma(
                                    89338.32554683842,
                                    eta2,
                                    fma(-2911.2154288136217, eta1, 16.437107575918503),
                                ),
                            ),
                        ),
                    ),
                )),
            s2,
            fma(
                -0.41839636169678796,
                fma(
                    -788783.0559069642,
                    eta6,
                    fma(
                        758333.2655404843,
                        eta5,
                        fma(
                            -283097.33358214336,
                            eta4,
                            fma(
                                50694.024801783446,
                                eta3,
                                fma(
                                    -4232.508174045782,
                                    eta2,
                                    fma(104.62892421207803, eta1, 5.143510231379954),
                                ),
                            ),
                        ),
                    ),
                ),
                -(0.05653522061311774
                    * fma(
                        -1.6010489769238676e6,
                        eta6,
                        fma(
                            1.7440315410044065e6,
                            eta5,
                            fma(
                                -761661.9645651339,
                                eta4,
                                fma(
                                    165236.6480734229,
                                    eta3,
                                    fma(
                                        -17551.398321516353,
                                        eta2,
                                        fma(694.00652410087, eta1, 5.605483124564013),
                                    ),
                                ),
                            ),
                        ),
                    )
                    * s1),
            ),
        ),
        sqroot,
        fma(
            fma(
                chidiff1 * delta,
                fma(
                    -34.98399465240273,
                    eta3,
                    fma(0.5697308729057493, eta1, 8.895576813118867 * eta2),
                ),
                chidiff2
                    * fma(
                        33.182723737396294,
                        eta3,
                        fma(1.6370346538130884, eta1, -(14.597095790380884 * eta2)),
                    ),
            ),
            sqroot,
            fma(
                -923.5495700703648,
                eta4,
                fma(
                    603.7453704122706,
                    eta3,
                    fma(
                        -138.9749850448088,
                        eta2,
                        fma(-3.557926105832778, eta1, 5.2601381002242595),
                    ),
                ),
            ) * sqroot,
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_44_iv2` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_44_iv2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let eta6 = eta1 * eta5;
    let s1 = s;
    let s2 = s1 * s1;
    let s3 = s1 * s2;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1 * fma(
            -(0.027167813927224657
                * fma(
                    1.8855858173287075e7,
                    eta6,
                    fma(
                        -1.6270212243584689e7,
                        eta5,
                        fma(
                            5.570748509263883e6,
                            eta4,
                            fma(
                                -956052.0170024392,
                                eta3,
                                fma(
                                    84126.14531363266,
                                    eta2,
                                    fma(-3312.7658728016568, eta1, 34.281932237450256),
                                ),
                            ),
                        ),
                    ),
                )),
            s3,
            fma(
                -(0.008830881730801855
                    * fma(
                        -7.756341721552802e6,
                        eta6,
                        fma(
                            5.858455766230802e6,
                            eta5,
                            fma(
                                -1.723865744472182e6,
                                eta4,
                                fma(
                                    260027.4121353132,
                                    eta3,
                                    fma(
                                        -23328.560285901796,
                                        eta2,
                                        fma(1359.7883958101172, eta1, -37.49992494976597),
                                    ),
                                ),
                            ),
                        ),
                    )),
                s2,
                fma(
                    -0.2627203100303006,
                    fma(
                        -234487.3875219032,
                        eta6,
                        fma(
                            227254.37897941095,
                            eta5,
                            fma(
                                -81525.50289542897,
                                eta4,
                                fma(
                                    12546.255587592654,
                                    eta3,
                                    fma(
                                        -552.1725902144143,
                                        eta2,
                                        fma(-52.82425783851536, eta1, 6.460396349297595),
                                    ),
                                ),
                            ),
                        ),
                    ),
                    -(0.008424003742397579
                        * fma(
                            -8.510662871580401e7,
                            eta6,
                            fma(
                                7.539643948937692e7,
                                eta5,
                                fma(
                                    -2.6446950627957724e7,
                                    eta4,
                                    fma(
                                        4.620165968920881e6,
                                        eta3,
                                        fma(
                                            -408022.6805482195,
                                            eta2,
                                            fma(15514.571912666677, eta1, -109.26773035716548),
                                        ),
                                    ),
                                ),
                            ),
                        )
                        * s1),
                ),
            ),
        ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    97.01712147349483,
                    eta3,
                    fma(4.749663394334708, eta1, -(42.62996105525792 * eta2)),
                )
                * s1,
            sqroot,
            fma(
                fma(
                    chidiff2,
                    fma(
                        -16.648660653220123,
                        eta3,
                        fma(-0.8318312659717388, eta1, 7.6541168007977864 * eta2),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            5.56587823143958,
                            eta3,
                            fma(2.214478316304753, eta1, -(7.028104574328955 * eta2)),
                        ),
                ),
                sqroot,
                fma(
                    -876.5048298377644,
                    eta4,
                    fma(
                        604.0067075996933,
                        eta3,
                        fma(
                            -155.22519772642607,
                            eta2,
                            fma(6.707695566702527, eta1, 3.173191054680422),
                        ),
                    ),
                ) * sqroot,
            ),
        ),
    )
}

/// `IMRPhenomXHM_Insp_Amp_44_iv3` (version 122022): inspiral amplitude
/// collocation point (strain units, f^(-7/6) included).
pub(crate) fn amp_44_iv3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let sqroot = eta.sqrt();
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let eta6 = eta1 * eta5;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1 * fma(
            -(0.026407497690308382
                * fma(
                    128167.89151782403,
                    eta6,
                    fma(
                        148102.4257785174,
                        eta5,
                        fma(
                            -184787.79377504133,
                            eta4,
                            fma(
                                66150.52694069289,
                                eta3,
                                fma(
                                    -10494.512487701073,
                                    eta2,
                                    fma(744.4743953122965, eta1, -17.184685557542196),
                                ),
                            ),
                        ),
                    ),
                )),
            s2,
            fma(
                -0.15782329022461472,
                fma(
                    307714.3954026392,
                    eta6,
                    fma(
                        -324891.4001578353,
                        eta5,
                        fma(
                            139319.8239886073,
                            eta4,
                            fma(
                                -31109.307191210424,
                                eta3,
                                fma(
                                    3777.477006415653,
                                    eta2,
                                    fma(-229.16361598098678, eta1, 6.0309399412954345),
                                ),
                            ),
                        ),
                    ),
                ),
                -(0.03050157254864058
                    * fma(
                        -1.138861230369246e7,
                        eta6,
                        fma(
                            1.0147654212772278e7,
                            eta5,
                            fma(
                                -3.5589766538499263e6,
                                eta4,
                                fma(
                                    612317.1751155312,
                                    eta3,
                                    fma(
                                        -51213.27604422822,
                                        eta2,
                                        fma(1609.4251694451375, eta1, 4.232861441291087),
                                    ),
                                ),
                            ),
                        ),
                    )
                    * s1),
            ),
        ),
        sqroot,
        fma(
            chidiff1
                * delta
                * fma(
                    53.35003268489743,
                    eta3,
                    fma(2.3532739003216254, eta1, -(21.37216554136868 * eta2)),
                )
                * s1,
            sqroot,
            fma(
                fma(
                    chidiff1 * delta,
                    fma(
                        -9.473290540936633,
                        eta3,
                        fma(1.4739380748149558, eta1, 0.06541707987699942 * eta2),
                    ),
                    chidiff2
                        * fma(
                            -8.709159662885131,
                            eta3,
                            fma(-0.3640838331639651, eta1, 3.7369795937033756 * eta2),
                        ),
                ),
                sqroot,
                fma(
                    -599.2709938848282,
                    eta4,
                    fma(
                        456.78649322753824,
                        eta3,
                        fma(
                            -139.6610487470118,
                            eta2,
                            fma(12.656614578053683, eta1, 1.7335503724888923),
                        ),
                    ),
                ) * sqroot,
            ),
        ),
    )
}

/// `IMRPhenomXHM_Insp_Phase_21_lambda` (version 122019): linear-in-f
/// correction to the rescaled inspiral phase for extreme mass ratios.
pub(crate) fn phase_21_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let s2 = s * s;
    let no_spin = fma(
        -133820.89317471132,
        eta5,
        fma(
            96463.68163125668,
            eta4,
            fma(
                -26847.690494515424,
                eta3,
                fma(
                    3535.657736681598,
                    eta2,
                    fma(-170.08866400251395, eta, 13.664473636545068),
                ),
            ),
        ),
    );
    let eq_spin = s * fma(
        eta2,
        fma(
            62149.00902591944,
            s2,
            fma(-7233.191207000735, s, -35031.361998480075),
        ),
        fma(
            33.8697137964237,
            s2,
            fma(
                eta,
                fma(
                    -6961.484805326852,
                    s2,
                    fma(1198.4163078715173, s, 3642.5891077598003),
                ),
                fma(
                    eta3,
                    fma(
                        -149700.4915210766,
                        s2,
                        fma(16501.749243703132, s, 83493.24265292779),
                    ),
                    fma(-41.55066592130464, s, 18.52571430563905),
                ),
            ),
        ),
    ) / fma(1.0, s, 6.880288191574696);
    let uneq_spin = -134.27742343186577 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Phase_33_lambda` (version 122019): linear-in-f
/// correction to the rescaled inspiral phase for extreme mass ratios.
pub(crate) fn phase_33_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let no_spin = fma(
        242.3428625556764,
        eta3,
        fma(
            -103.92956504520747,
            eta2,
            fma(9.772510519809892, eta, 4.1138398568400705),
        ),
    );
    let eq_spin = fma(
        -105.09339163109497,
        eta2,
        fma(26.644159828590055, eta, -0.13253553909611435),
    ) * s
        / fma(0.11322426762297967, s, 1.0);
    let uneq_spin = -19.705359163581168 * wf.dchi * eta2 * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Phase_32_lambda` (version 122019): linear-in-f
/// correction to the rescaled inspiral phase for extreme mass ratios.
pub(crate) fn phase_32_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -6244.001932443913,
        eta4,
        fma(
            2671.7813055097877,
            eta3,
            fma(
                -574.8672384388947,
                eta2,
                fma(18.424900617803107, eta, 9.913819875501506),
            ),
        ),
    ) / fma(-0.9103118343073325, eta, 1.0);
    let eq_spin = fma(
        fma(
            42379.205045791656,
            eta3,
            fma(
                -18636.271414571278,
                eta2,
                fma(2009.7361967331492, eta, 1.8403161863444328),
            ),
        ),
        s4,
        fma(
            fma(
                -273.23750933544176,
                eta3,
                fma(
                    708.1455154504333,
                    eta2,
                    fma(-138.36970336658558, eta, -1.28324755577382),
                ),
            ),
            s3,
            fma(
                fma(
                    5894.355429022858,
                    eta3,
                    fma(
                        -2233.9319708029775,
                        eta2,
                        fma(245.06757304950986, eta, -4.367632806613781),
                    ),
                ),
                s,
                fma(
                    -40928.07304790013,
                    eta3,
                    fma(
                        17608.172965575013,
                        eta2,
                        fma(-1876.760129419146, eta, -1.375112297530783),
                    ),
                ) * s2,
            ),
        ),
    );
    let uneq_spin = wf.dchi
        * fma(-4.0, eta, 1.0).sqrt()
        * eta2
        * fma(
            2155.472229664981 * eta,
            s,
            fma(
                1566.1242344157668 * wf.chi2l,
                eta,
                fma(-(1566.1242344157668 * wf.chi1l), eta, -105.34550407768225),
            ),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Insp_Phase_44_lambda` (version 122019): linear-in-f
/// correction to the rescaled inspiral phase for extreme mass ratios.
pub(crate) fn phase_44_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let s2 = s * s;
    let no_spin = fma(
        -1538.6661348106031,
        eta5,
        fma(
            1685.5912722190276,
            eta4,
            fma(
                -1162.954360723399,
                eta3,
                fma(
                    160.43721442910618,
                    eta2,
                    fma(-21.277760168559862, eta, 5.254484747463392),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            -810.7773268493444,
            eta3,
            fma(
                246.8787141453734,
                eta2,
                fma(-10.945895160727437, eta, 0.007067861615983771),
            ),
        ),
        s,
        fma(
            621.6920322846844,
            eta3,
            fma(
                -176.4987316167203,
                eta2,
                fma(4.530539154777984, eta, 0.17447830920234977),
            ),
        ) * s2,
    );
    let uneq_spin = -8.384066369867833 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/* End of the parameter-space fits. */

/// Inspiral amplitude collocation fits `iv1..iv3` of mode `mode_int`.
fn amp_fit(mode_int: usize, i: usize, wf: &WaveformStruct) -> f64 {
    const FITS: [[Fit; 3]; 4] = [
        [amp_21_iv1, amp_21_iv2, amp_21_iv3],
        [amp_33_iv1, amp_33_iv2, amp_33_iv3],
        [amp_32_iv1, amp_32_iv2, amp_32_iv3],
        [amp_44_iv1, amp_44_iv2, amp_44_iv3],
    ];
    FITS[mode_int][i](wf)
}

/// `pPhase->InspiralPhaseFits[modeInt]`: the linear-in-f phase term for
/// extreme mass ratios.
pub(crate) fn phase_lambda_fit(mode_int: usize, wf: &WaveformStruct) -> f64 {
    const FITS: [Fit; 4] = [
        phase_21_lambda,
        phase_33_lambda,
        phase_32_lambda,
        phase_44_lambda,
    ];
    FITS[mode_int](wf)
}

/// `IMRPhenomXHM_Insp_Phase_LambdaPN`: linear-in-f term from the complex PN
/// amplitude added to the rescaled orbital phase, Eq. 4.12 of
/// arXiv:2001.10914.
pub(crate) fn lambda_pn(eta: f64, mode_tag: u32) -> f64 {
    let pi = POWERS_OF_LALPI.itself;
    let output = match mode_tag {
        21 => 2. * pi * fma(-2., LN_2, -0.5),
        33 => 2. / 3. * pi * fma(6., 1.5f64.ln(), -21. / 5.),
        32 => -((2376. * pi * fma(22., eta, -5.)) / fma(11880., eta, -3960.)),
        44 => {
            45045. * pi * fma(320. * fma(3., eta, -1.), LN_2, fma(-1193., eta, 336.))
                / (2. * fma(5405400., eta, -1801800.))
        }
        _ => 0.,
    };
    -output
}

/// `IMRPhenomXHM_GetPNAmplitudeCoefficients`: coefficients of the
/// re-expanded Fourier-domain PN amplitude (power series in f^(1/3)) and the
/// leading PN factor of the mode relative to the 22.
///
/// Some coefficients are complex. LAL evaluates them in C complex
/// arithmetic, where sums of complex and real operands are not contracted
/// into fused multiply-adds; the real and imaginary parts below follow that
/// order of operations.
pub(crate) fn get_pn_amplitude_coefficients(
    amp: &mut HMAmpCoefficients,
    wfhm: &HMWaveformStruct,
    wf: &WaveformStruct,
) {
    let chi_a = wfhm.chi_a;
    let chi_s = wfhm.chi_s;
    let eta = wf.eta;
    let delta = wf.delta;
    let lalpi = &*POWERS_OF_LALPI;
    let pi = lalpi.itself;
    let sqrt_5_7 = (5. / 7f64).sqrt();
    let prefactors = [
        2f64.sqrt() / 3.,
        0.75 * sqrt_5_7,
        sqrt_5_7 / 3.,
        4. * 2f64.sqrt() / 9. * sqrt_5_7,
    ];
    amp.pn_global_factor = (2. / wfhm.emm as f64).powf(-7. / 6.) * prefactors[wfhm.mode_int];
    let real = |x: f64| Complex::new(x, 0.0);
    let zero = real(0.0);
    amp.pn = match wfhm.mode_tag {
        21 => {
            let p2 = UsefulPowers::new(2.);
            let four_thirds = {
                let mut re = 3427. * chi_a;
                re += 3427. * chi_s * delta;
                re -= 8404. * chi_a * eta;
                re -= 3860. * chi_s * delta * eta;
                re -= 1344. * delta * pi;
                let mut im = -(672. * delta);
                im -= 672. * delta * 16f64.ln();
                let scale = |x: f64| x / 1344. * lalpi.four_thirds * p2.four_thirds;
                Complex::new(scale(re), scale(im))
            };
            let six_thirds = {
                let mut re = 143063173. * chi_a;
                re += 143063173. * chi_s * delta;
                re -= 546199608. * chi_a * eta;
                re -= 169191096. * chi_s * delta * eta;
                re -= 9898560. * delta * pi;
                re += 20176128. * delta * eta * pi;
                re += 61725888. * chi_s * delta * chi_a * chi_a;
                re -= 81285120. * chi_s * delta * eta * chi_a * chi_a;
                re += 20575296. * chi_a.powf(3.);
                re -= 81285120. * eta * chi_a.powf(3.);
                re += 61725888. * chi_a * chi_s * chi_s;
                re -= 165618432. * chi_a * eta * chi_s * chi_s;
                re += 20575296. * delta * chi_s.powf(3.);
                re -= 1016064. * delta * eta * chi_s * chi_s * chi_s;
                re += 128873808. * chi_a * eta * eta;
                re -= 3859632. * chi_s * delta * eta * eta;
                let mut im = -(1350720. * delta);
                im -= 72043776. * delta * eta;
                im -= 5402880. * delta * LN_2;
                im -= 17224704. * delta * eta * LN_2;
                let scale = |x: f64| x / 5.419008e6 * lalpi.two * p2.two;
                Complex::new(scale(re), scale(im))
            };
            [
                zero,
                real(delta * lalpi.one_third * p2.one_third),
                real(-3.0 * fma(chi_s, delta, chi_a) / 2.0 * lalpi.two_thirds * p2.two_thirds),
                real(fma(335.0, delta, 1404.0 * delta * eta) / 672.0 * lalpi.itself * p2.itself),
                four_thirds,
                real(
                    fma(
                        21535920.0 * delta * eta,
                        eta,
                        fma(
                            39626496.0 * delta * eta * chi_s,
                            chi_s,
                            fma(
                                -(77982912.0 * delta * chi_s),
                                chi_s,
                                fma(
                                    81285120.0 * delta * eta * chi_a,
                                    chi_a,
                                    fma(
                                        -(77982912.0 * delta * chi_a),
                                        chi_a,
                                        fma(
                                            24385536.0 * chi_s * delta,
                                            pi,
                                            fma(
                                                24385536.0 * chi_a,
                                                pi,
                                                fma(
                                                    -(23670792.0 * delta),
                                                    eta,
                                                    fma(
                                                        432843264.0 * chi_a * chi_s,
                                                        eta,
                                                        fma(
                                                            -155965824.0 * chi_a,
                                                            chi_s,
                                                            -(964357.0 * delta),
                                                        ),
                                                    ),
                                                ),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ) / 8.128512e6
                        * lalpi.five_thirds
                        * p2.five_thirds,
                ),
                six_thirds,
            ]
        }
        33 => {
            let p2 = UsefulPowers::new(2. / 3.);
            let four_thirds = {
                let mut re = 325. * chi_a;
                re += 325. * chi_s * delta;
                re -= 1120. * chi_a * eta;
                re -= 80. * chi_s * delta * eta;
                re += 120. * delta * pi;
                let mut im = -(504. * delta);
                im += 720. * delta * 1.5f64.ln();
                let scale = |x: f64| x / 120. * lalpi.four_thirds * p2.four_thirds;
                Complex::new(scale(re), scale(im))
            };
            let six_thirds = {
                let mut re = 22007835. * chi_a;
                re += 22007835. * chi_s * delta;
                re -= 80190540. * chi_a * eta;
                re -= 31722300. * chi_s * delta * eta;
                re -= 9193500. * delta * pi;
                re += 17826480. * delta * eta * pi;
                re -= 12428640. * chi_a * eta * eta;
                re -= 6078240. * chi_s * delta * eta * eta;
                let mut im = 26467560. * delta;
                im -= 98774368. * delta * eta;
                im -= 37810800. * delta * 1.5f64.ln();
                im += 37558080. * delta * eta * 1.5f64.ln();
                let scale = |x: f64| x / 2.17728e6 * lalpi.two * p2.two;
                Complex::new(scale(re), scale(im))
            };
            [
                zero,
                real(delta * lalpi.one_third * p2.one_third),
                zero,
                real(fma(-1945.0, delta, 2268.0 * delta * eta) / 672.0 * lalpi.itself * p2.itself),
                four_thirds,
                real(
                    fma(
                        2966264784.0 * delta * eta,
                        eta,
                        fma(
                            55883520.0 * delta * eta * chi_s,
                            chi_s,
                            fma(
                                -(1131641280.0 * delta * chi_s),
                                chi_s,
                                fma(
                                    4470681600.0 * delta * eta * chi_a,
                                    chi_a,
                                    fma(
                                        -(1131641280.0 * delta * chi_a),
                                        chi_a,
                                        fma(
                                            -(5926068792.0 * delta),
                                            eta,
                                            fma(
                                                9053130240.0 * chi_a * chi_s,
                                                eta,
                                                fma(
                                                    -2263282560.0 * chi_a,
                                                    chi_s,
                                                    -(1077664867.0 * delta),
                                                ),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ) / 4.4706816e8
                        * lalpi.five_thirds
                        * p2.five_thirds,
                ),
                six_thirds,
            ]
        }
        32 => {
            let five_thirds = {
                let mut re = 0.0 - 3955. * chi_s;
                re -= 3955. * chi_a * delta;
                re += 10810. * chi_s * eta;
                re += 11865. * chi_a * delta * eta;
                re -= 12600. * chi_s * eta * eta;
                let mut im = 2520.;
                im -= 11088. * eta;
                let scale = |x: f64| x / 840. * lalpi.five_thirds;
                Complex::new(scale(re), scale(im))
            };
            [
                zero,
                zero,
                real(fma(3.0, eta, -1.0) * lalpi.two_thirds),
                real(-4.0 * chi_s * eta * lalpi.itself),
                real(
                    fma(82460.0 * eta, eta, fma(-61625.0, eta, 10471.0)) / 10080.0
                        * lalpi.four_thirds,
                ),
                five_thirds,
                real(
                    fma(
                        5907445488.0,
                        eta.powf(3.0),
                        fma(
                            5830513920.0 * chi_s * chi_s * eta,
                            eta,
                            fma(
                                13412044800.0 * chi_a * chi_a * eta,
                                eta,
                                fma(
                                    -(13202119896.0 * eta),
                                    eta,
                                    fma(
                                        -(11870591040.0 * eta * chi_s),
                                        chi_s,
                                        fma(
                                            1131641280.0 * chi_s,
                                            chi_s,
                                            fma(
                                                -(7865605440.0 * eta * chi_a),
                                                chi_a,
                                                fma(
                                                    1131641280.0 * chi_a,
                                                    chi_a,
                                                    fma(
                                                        3576545280.0 * chi_s * eta,
                                                        pi,
                                                        fma(
                                                            -(15209631360.0
                                                                * chi_a
                                                                * chi_s
                                                                * delta),
                                                            eta,
                                                            fma(
                                                                -26069649.0,
                                                                eta,
                                                                fma(
                                                                    2263282560.0 * chi_a * chi_s,
                                                                    delta,
                                                                    824173699.0,
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
                    ) / 4.4706816e8
                        * lalpi.two,
                ),
            ]
        }
        _ => {
            let p2 = UsefulPowers::new(0.5);
            let five_thirds = {
                let mut re = -0.0 + 565. * chi_s;
                re += 565. * chi_a * delta;
                re -= 2075. * chi_s * eta;
                re -= 1695. * chi_a * delta * eta;
                re += 240. * pi;
                re -= 720. * eta * pi;
                re += 1140. * chi_s * eta * eta;
                let mut im = -1008.;
                im += 3579. * eta;
                im += 960. * LN_2;
                im -= 2880. * eta * LN_2;
                let scale = |x: f64| x / 120. * lalpi.five_thirds * p2.five_thirds;
                Complex::new(scale(re), scale(im))
            };
            [
                zero,
                zero,
                real(fma(-3.0, eta, 1.0) * lalpi.two_thirds * p2.two_thirds),
                zero,
                real(
                    fma(-(446460.0 * eta), eta, fma(641105.0, eta, -158383.0)) / 36960.0
                        * lalpi.four_thirds
                        * p2.four_thirds,
                ),
                five_thirds,
                real(
                    fma(
                        -805075876080.0,
                        eta.powf(3.0),
                        fma(
                            -(10897286400.0 * chi_s * chi_s * eta),
                            eta,
                            fma(
                                -(871782912000.0 * chi_a * chi_a * eta),
                                eta,
                                fma(
                                    2271682065240.0 * eta,
                                    eta,
                                    fma(
                                        224302478400.0 * eta * chi_s,
                                        chi_s,
                                        fma(
                                            -(73556683200.0 * chi_s),
                                            chi_s,
                                            fma(
                                                511264353600.0 * eta * chi_a,
                                                chi_a,
                                                fma(
                                                    -(73556683200.0 * chi_a),
                                                    chi_a,
                                                    fma(
                                                        441340099200.0 * chi_a * chi_s * delta,
                                                        eta,
                                                        fma(
                                                            -745140957231.0,
                                                            eta,
                                                            fma(
                                                                -(147113366400.0 * chi_a * chi_s),
                                                                delta,
                                                                7888301437.0,
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
                    ) / 2.90594304e10
                        * lalpi.two
                        * p2.two,
                ),
            ]
        }
    };
}

/// `IMRPhenomXHM_Inspiral_PNAmp_Ansatz` (power series, rescale factor 0):
/// the PN amplitude in strain units without `amp0`.
#[inline]
pub(crate) fn pn_amp_ansatz(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    let pn = &amp.pn;
    // C complex arithmetic: componentwise products and sums, no contraction.
    let terms = [
        powers.one_third,
        powers.two_thirds,
        powers.itself,
        powers.four_thirds,
        powers.five_thirds,
        powers.two,
    ];
    let mut re = pn[0].re;
    let mut im = pn[0].im;
    for (p, c) in terms.iter().zip(&pn[1..]) {
        re += p * c.re;
        im += p * c.im;
    }
    re.hypot(im) * (amp.pn_global_factor * powers.m_seven_sixths * wfhm.amp_norm)
}

/// `IMRPhenomXHM_Inspiral_Amp_Ansatz` (122022 release): PN amplitude plus
/// pseudo-PN terms, in strain units without `amp0`.
#[inline]
pub(crate) fn amp_ansatz(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    let insp = pn_amp_ansatz(powers, wfhm, amp);
    let c = &amp.inspiral_coefficients;
    let mut pseudo = fma(
        powers.seven_thirds / amp.fcut_insp_seven_thirds,
        c[0],
        powers.eight_thirds / amp.fcut_insp_eight_thirds * c[1],
    );
    pseudo = fma(powers.three / amp.fcut_insp_three, c[2], pseudo);
    pseudo *= powers.m_seven_sixths * amp.pn_dominant;
    insp + pseudo
}

/// `|pn[0] + pn[1] f^(1/3) + ... + pn[6] f^2|` with its first and second
/// derivatives.
pub(crate) fn pn_modulus_jet(ff: f64, pn: &[Complex<f64>; 7]) -> Jet {
    let mut re = Jet::constant(pn[0].re);
    let mut im = Jet::constant(pn[0].im);
    for (k, c) in pn.iter().enumerate().skip(1) {
        let p = Jet::power(ff, k as f64 / 3.);
        re = re + p * c.re;
        im = im + p * c.im;
    }
    re.hypot(im)
}

/// [`amp_ansatz`] with its first and second derivatives.
fn amp_ansatz_jet(ff: f64, wfhm: &HMWaveformStruct, amp: &HMAmpCoefficients) -> Jet {
    let m_seven_sixths = Jet::power(ff, -7. / 6.);
    let insp =
        pn_modulus_jet(ff, &amp.pn) * m_seven_sixths * (amp.pn_global_factor * wfhm.amp_norm);
    let c = &amp.inspiral_coefficients;
    let pseudo = Jet::power(ff, 7. / 3.) * (c[0] / amp.fcut_insp_seven_thirds)
        + Jet::power(ff, 8. / 3.) * (c[1] / amp.fcut_insp_eight_thirds)
        + Jet::power(ff, 3.) * (c[2] / amp.fcut_insp_three);
    insp + pseudo * m_seven_sixths * amp.pn_dominant
}

/// `IMRPhenomXHM_Inspiral_Amp_NDAnsatz`: derivative of the inspiral
/// amplitude. LALSimulation takes it by fourth-order finite differences
/// with step 1e-9; this is the analytic derivative they approximate.
pub(crate) fn amp_nd_ansatz(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    amp_ansatz_jet(powers.itself, wfhm, amp).d1
}

/// `IMRPhenomXHM_Get_Inspiral_Amp_Coefficients` with
/// `IMRPhenomXHM_Inspiral_Amp_CollocationPoints` and
/// `IMRPhenomXHM_Inspiral_Amp_Coefficients` (122022 release, three
/// collocation points).
pub(crate) fn get_inspiral_amp_coefficients(
    amp: &mut HMAmpCoefficients,
    wfhm: &HMWaveformStruct,
    wf: &WaveformStruct,
) {
    let fcut = amp.f_amp_match_in;
    amp.colloc_freqs_insp = [0.5 * fcut, 0.75 * fcut, fcut];
    for i in 0..3 {
        amp.colloc_values_insp[i] = amp_fit(wfhm.mode_int, i, wf).abs();
    }
    let f1 = UsefulPowers::new(amp.colloc_freqs_insp[0]);
    let f2 = UsefulPowers::new(amp.colloc_freqs_insp[1]);
    let f3 = UsefulPowers::new(amp.colloc_freqs_insp[2]);
    let finsp = UsefulPowers::new(fcut);
    amp.fcut_insp_seven_thirds = finsp.seven_thirds;
    amp.fcut_insp_eight_thirds = finsp.eight_thirds;
    amp.fcut_insp_three = finsp.three;
    for (i, powers) in [&f1, &f2, &f3].into_iter().enumerate() {
        amp.pn_amplitude_insp[i] = pn_amp_ansatz(powers, wfhm, amp);
    }

    let v = |i: usize, powers: &UsefulPowers| {
        (amp.colloc_values_insp[i] - amp.pn_amplitude_insp[i])
            / amp.pn_dominant
            / powers.m_seven_sixths
    };
    let (v1, v2, v3) = (v(0, &f1), v(1, &f2), v(2, &f3));
    let denominator = f1.seven_thirds
        * (f1.one_third - f2.one_third)
        * f2.seven_thirds
        * (f1.one_third - f3.one_third)
        * (f2.one_third - f3.one_third)
        * f3.seven_thirds;
    let c0 = finsp.seven_thirds
        * fma(
            f2.eight_thirds,
            fma(-f3.three, v1, f1.three * v3),
            fma(
                f2.three,
                fma(f3.eight_thirds, v1, -(f1.eight_thirds * v3)),
                fma(
                    -(f1.three * f3.eight_thirds),
                    v2,
                    f1.eight_thirds * f3.three * v2,
                ),
            ),
        )
        / denominator;
    let c1 = finsp.eight_thirds
        * fma(
            f2.seven_thirds,
            fma(f3.three, v1, -(f1.three * v3)),
            fma(
                f2.three,
                fma(-f3.seven_thirds, v1, f1.seven_thirds * v3),
                fma(
                    f1.three * f3.seven_thirds,
                    v2,
                    -(f1.seven_thirds * f3.three * v2),
                ),
            ),
        )
        / denominator;
    let c2 = finsp.three
        * fma(
            f2.eight_thirds,
            fma(f3.seven_thirds, v1, -(f1.seven_thirds * v3)),
            fma(
                f1.seven_thirds * (-f1.one_third + f3.one_third) * f3.seven_thirds,
                v2,
                f2.seven_thirds * fma(-f3.eight_thirds, v1, f1.eight_thirds * v3),
            ),
        )
        / denominator;
    amp.inspiral_coefficients = [c0, c1, c2];
}

/// `IMRPhenomXHM_Inspiral_Phase_AnsatzInt`: the rescaled 22 inspiral phase.
#[inline]
pub(crate) fn phase_ansatz_int(mf: f64, powers: &UsefulPowers, phase: &HMPhaseCoefficients) -> f64 {
    let freqs = [
        powers.m_five_thirds,
        powers.m_four_thirds,
        powers.m_one,
        powers.m_two_thirds,
        powers.m_one_third,
        1.,
        powers.one_third,
        powers.two_thirds,
        mf,
        powers.four_thirds,
        powers.five_thirds,
        powers.two,
        powers.seven_thirds,
    ];
    let log_mf = powers.log;
    let mut philm = 0.;
    for ((phi, phi_l), f) in phase.phi.iter().zip(&phase.phi_l).zip(freqs) {
        philm = fma(fma(*phi_l, log_mf, *phi), f, philm);
    }
    philm
}

/// `IMRPhenomXHM_Inspiral_Phase_Ansatz`: derivative of
/// [`phase_ansatz_int`].
pub(crate) fn phase_ansatz(mf: f64, powers: &UsefulPowers, phase: &HMPhaseCoefficients) -> f64 {
    const COEFFS: [f64; 13] = [
        -5. / 3.,
        -4. / 3.,
        -1.,
        -2. / 3.,
        -1. / 3.,
        0.,
        1. / 3.,
        2. / 3.,
        1.,
        4. / 3.,
        5. / 3.,
        2.,
        7. / 3.,
    ];
    let freqs = [
        powers.m_eight_thirds,
        powers.m_seven_thirds,
        powers.m_two,
        powers.m_five_thirds,
        powers.m_four_thirds,
        powers.m_one,
        powers.m_two_thirds,
        powers.m_one_third,
        1.,
        powers.one_third,
        powers.two_thirds,
        mf,
        powers.four_thirds,
    ];
    let log_mf = powers.log;
    let mut dphilm = 0.;
    for i in 0..13 {
        let coefficient = fma(phase.phi_l[i], log_mf, phase.phi[i]);
        dphilm = fma(
            fma(coefficient, COEFFS[i], phase.phi_l[i]),
            freqs[i],
            dphilm,
        );
    }
    dphilm
}
