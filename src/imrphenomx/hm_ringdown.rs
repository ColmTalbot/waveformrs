// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Ringdown region of the IMRPhenomXHM higher modes,
//! `LALSimIMRPhenomXHM_ringdown.c`: parameter-space fits of the amplitude
//! (version 122022) and phase (version 122019), the ringdown amplitude
//! coefficients, the spheroidal (3,2) ringdown and the ansatzes.

use crate::imrphenomx::hm_internals::{
    Fit, HMAmpCoefficients, HMPhaseCoefficients, HMWaveformStruct,
};
use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct};
use crate::imrphenomx::jet::Jet;
use crate::imrphenomx::linalg::lu_solve;

/// `IMRPhenomXHM_RD_Amp_21_rdcp1` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_21_rdcp1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
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
    fma(
        chidiff1 * eta5,
        fma(
            10.13490934572329,
            s2,
            fma(175.59081294852444, s1, -1512.439342647443),
        ),
        fma(
            delta * eta1 * s1,
            fma(
                -10.090867481062709,
                fma(
                    -1992.9074348833092,
                    eta6,
                    fma(
                        3950.7894052628735,
                        eta5,
                        fma(
                            -2688.211889175019,
                            eta4,
                            fma(
                                801.3948756800821,
                                eta3,
                                fma(
                                    -107.22993216128548,
                                    eta2,
                                    fma(5.388149112485179, eta1, 0.9580746052260011),
                                ),
                            ),
                        ),
                    ),
                ),
                -(0.42972412296628143
                    * fma(
                        658866.0983367155,
                        eta6,
                        fma(
                            -479649.75188253267,
                            eta5,
                            fma(
                                107980.65459735804,
                                eta4,
                                fma(
                                    -3176.4950303461164,
                                    eta3,
                                    fma(
                                        -1616.9974609915555,
                                        eta2,
                                        fma(139.73149069609775, eta1, 1.9193131231064235),
                                    ),
                                ),
                            ),
                        ),
                    )
                    * s1),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        591.1097623278739,
                        eta3,
                        fma(21.402708785047853, eta1, -(232.07306353130417 * eta2)),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(
                        -175.16681482428694,
                        eta3,
                        fma(
                            92.6090002736012,
                            eta2,
                            fma(-23.5291063016996, eta1, 12.880905080761432),
                        ),
                    ),
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            2255.040486907459,
                            eta3,
                            fma(26.89427230731867, eta1, -(710.8871223808559 * eta2)),
                        ),
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_RD_Amp_21_rdcp2` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_21_rdcp2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
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
        delta * eta1 * s1,
        fma(
            0.40226131643223145
                * fma(
                    -29331.524552164105,
                    eta4,
                    fma(
                        16087.433824017446,
                        eta3,
                        fma(
                            -2885.5110686004946,
                            eta2,
                            fma(190.42214703482531, eta1, -3.834742668014861),
                        ),
                    ),
                ),
            s2,
            fma(
                -8.236485576091717,
                fma(
                    -997.6961311502365,
                    eta4,
                    fma(
                        519.0146702141192,
                        eta3,
                        fma(
                            -87.05136337926156,
                            eta2,
                            fma(5.1501231412520285, eta1, 0.8917610178208336),
                        ),
                    ),
                ),
                0.2836840678615208
                    * fma(
                        2330.3534917059906,
                        eta4,
                        fma(
                            -1882.2040277496196,
                            eta3,
                            fma(
                                586.7942442434971,
                                eta2,
                                fma(-57.65586769647737, eta1, -0.19281297100324718),
                            ),
                        ),
                    )
                    * s1,
            ),
        ),
        fma(
            chidiff1 * eta5,
            fma(
                21.159179960295855,
                s2,
                fma(245.73200219767514, s1, -1116.4039232324135),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        -20.274691328125606,
                        eta3,
                        fma(2.7866868976718226, eta1, -(0.4210629980868266 * eta2)),
                    ),
                s1,
                fma(
                    delta * fma(-7.5304766811877455, eta1, 9.112452928978168),
                    eta1,
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            1618.0818430353293,
                            eta3,
                            fma(16.236533863306132, eta1, -(500.11964987628926 * eta2)),
                        ),
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_RD_Amp_21_rdcp3` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_21_rdcp3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
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
        delta * eta1 * s1,
        fma(
            0.1930861073906724
                * fma(
                    485.67801863289293,
                    eta3,
                    fma(
                        -78.07256429516346,
                        eta2,
                        fma(-19.11543562444476, eta1, 3.7735099269174106),
                    ),
                ),
            s2,
            fma(
                -2.7294297839371824,
                fma(
                    -87.26454353763077,
                    eta3,
                    fma(
                        36.120093043420326,
                        eta2,
                        fma(-4.384077347340523, eta1, 1.148166706456899),
                    ),
                ),
                0.23949142867803436
                    * fma(
                        862.3123076782916,
                        eta3,
                        fma(
                            -307.3404155231787,
                            eta2,
                            fma(33.33372867559165, eta1, -0.6931516433988293),
                        ),
                    )
                    * s1,
            ),
        ),
        fma(
            chidiff1 * eta5,
            fma(
                -5.353399481074393,
                s2,
                fma(92.91301644484749, s1, -334.0969956136684),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        -0.7347258030219584,
                        eta3,
                        fma(1.9701247529688362, eta1, -(2.8616711550845575 * eta2)),
                    ),
                s1,
                fma(
                    delta * fma(-3.038523690239521, eta1, 2.920930733198033),
                    eta1,
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            544.1978232314333,
                            eta3,
                            fma(6.3472251472354975, eta1, -(171.23657247338042 * eta2)),
                        ),
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_RD_Amp_33_rdcp1` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_33_rdcp1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
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
        13.499999999999998
            * delta
            * eta1
            * s1
            * fma(
                -0.13654149379906394,
                fma(
                    -3591.1757995710486,
                    eta5,
                    fma(
                        -6167.08766058184,
                        eta4,
                        fma(
                            4142.974510926698,
                            eta3,
                            fma(
                                -742.1357702210267,
                                eta2,
                                fma(29.023992126142304, eta1, 2.719687834084113),
                            ),
                        ),
                    ),
                ),
                -(0.06248535354306988
                    * fma(
                        -10056.98422430965,
                        eta5,
                        fma(
                            6601.918552659412,
                            eta4,
                            fma(
                                -1907.008984765889,
                                eta3,
                                fma(
                                    444.79350113344543,
                                    eta2,
                                    fma(-78.23231700361792, eta1, 6.697567446351289),
                                ),
                            ),
                        ),
                    )
                    * s1),
            ),
        1.0 / (-3.9329308614837704 + s1),
        fma(
            chidiff1 * eta5,
            fma(
                -15.81624534976308,
                s2,
                fma(
                    -21.836905132600755,
                    s1,
                    fma(-34.5074027042393, chidiff2, 257.33227387984863),
                ),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        -1653.496174539196,
                        eta3,
                        fma(-52.85961155799673, eta1, 577.6347407795782 * eta2),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(
                        22.780673360839497,
                        eta2,
                        fma(-4.436329538596615, eta1, 12.439702602599235),
                    ),
                    delta
                        * eta1
                        * fma(
                            chidiff1,
                            fma(
                                -1524.2981907688634,
                                eta3,
                                fma(-41.04442169938298, eta1, 502.9246970179746 * eta2),
                            ),
                            chidiff2
                                * fma(
                                    1020.6734178547847,
                                    eta3,
                                    fma(32.23960072974939, eta1, -(365.1526474476759 * eta2)),
                                ),
                        ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_33_rdcp2` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_33_rdcp2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
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
        chidiff1 * eta5,
        fma(
            41.321480132899524,
            s2,
            fma(
                -1.873799747678187,
                s1,
                fma(-29.10412708633528, chidiff2, 159.72371180117415),
            ),
        ),
        fma(
            delta * eta1 * s1,
            fma(
                0.1206817303851239,
                fma(
                    -265581.70142471837,
                    eta5,
                    fma(
                        163644.2192178668,
                        eta4,
                        fma(
                            -35378.156133055556,
                            eta3,
                            fma(
                                3188.189172446398,
                                eta2,
                                fma(-144.08062755162752, eta1, 8.667503604073314),
                            ),
                        ),
                    ),
                ),
                0.08028332044013944
                    * fma(
                        -148782.33687815256,
                        eta5,
                        fma(
                            121293.97832549023,
                            eta4,
                            fma(
                                -35625.58409457366,
                                eta3,
                                fma(
                                    4777.45310151897,
                                    eta2,
                                    fma(-322.95832000179297, eta1, 12.632478544060636),
                                ),
                            ),
                        ),
                    )
                    * s1,
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        -636.4802385540647,
                        eta3,
                        fma(-17.18006888428382, eta1, 190.73514518113845 * eta2),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(4.543696144846763, eta1, 8.425057692276933),
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            -1293.422289802462,
                            eta3,
                            fma(-32.18860840414171, eta1, 412.07321398189293 * eta2),
                        ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_33_rdcp3` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_33_rdcp3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
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
        13.499999999999998
            * delta
            * eta1
            * s1
            * fma(
                0.0011484326782460882,
                fma(
                    -517034.11171277676,
                    eta5,
                    fma(
                        316864.0589150127,
                        eta4,
                        fma(
                            -67368.27397765424,
                            eta3,
                            fma(
                                5709.239330076732,
                                eta2,
                                fma(-172.58816646768219, eta1, 4.1815722950796035),
                            ),
                        ),
                    ),
                ),
                -(0.009496797093329243
                    * fma(
                        -140089.40725211444,
                        eta5,
                        fma(
                            94309.57566199072,
                            eta4,
                            fma(
                                -23464.64953722729,
                                eta3,
                                fma(
                                    2628.6024206791726,
                                    eta2,
                                    fma(-118.35865186626413, eta1, 0.9233282181397624),
                                ),
                            ),
                        ),
                    )
                    * s1),
            ),
        1.0 / (fma(-0.41099904730526465, s1, 0.09549360183532198) + s2),
        fma(
            chidiff1 * eta5,
            fma(
                67.23189816732847,
                s2,
                fma(
                    35.09131921815571,
                    s1,
                    fma(-3.071422453072518, chidiff2, 32.07819766300362),
                ),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        -473.30988537734174,
                        eta3,
                        fma(-10.27418232676514, eta1, 136.5150165348149 * eta2),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(2.321696430921996, eta1, 2.485784720088995),
                    delta
                        * eta1
                        * fma(
                            chidiff1,
                            fma(
                                -496.1564538739011,
                                eta3,
                                fma(-10.454376404653859, eta1, 147.10344302665484 * eta2),
                            ),
                            chidiff2
                                * fma(
                                    -197.51205149250532,
                                    eta3,
                                    fma(-5.9236399792925996, eta1, 65.86115501723127 * eta2),
                                ),
                        ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_44_rdcp1` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_44_rdcp1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
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
        eta1 * s1,
        fma(
            0.005856387555754387
                * fma(
                    -3.2369858387064277e6,
                    eta5,
                    fma(
                        3.054333578686424e6,
                        eta4,
                        fma(
                            -1.0461268207440731e6,
                            eta3,
                            fma(
                                156431.3764198244,
                                eta2,
                                fma(-9560.410655118145, eta1, 186.39698091707513),
                            ),
                        ),
                    ),
                ),
            s2,
            fma(
                -0.4360849737360132,
                fma(
                    -111910.71267453219,
                    eta5,
                    fma(
                        71236.86532610047,
                        eta4,
                        fma(
                            -16718.425586396803,
                            eta3,
                            fma(
                                1729.1839588870455,
                                eta2,
                                fma(-58.70494649755802, eta1, -0.9543114627170375),
                            ),
                        ),
                    ),
                ),
                -(0.024861802943501172
                    * fma(
                        -581776.5303770923,
                        eta5,
                        fma(
                            168937.32229060197,
                            eta4,
                            fma(
                                35332.328181283,
                                eta3,
                                fma(
                                    -15866.093368857853,
                                    eta2,
                                    fma(1585.462602954658, eta1, -52.25045490410733),
                                ),
                            ),
                        ),
                    )
                    * s1),
            ),
        ),
        fma(
            eta1,
            fma(
                chidiff1 * delta,
                fma(
                    -297.2592736781142,
                    eta3,
                    fma(-8.51952446214978, eta1, 117.76530248141987 * eta2),
                ),
                chidiff2
                    * fma(
                        -8.017569928870929,
                        eta3,
                        fma(-0.2750098647982238, eta1, 4.456900599347149 * eta2),
                    ),
            ),
            eta1 * fma(
                158292.58792484726,
                eta6,
                fma(
                    -98374.18361532023,
                    eta5,
                    fma(
                        25108.811524802128,
                        eta4,
                        fma(
                            -3514.3385364216438,
                            eta3,
                            fma(
                                287.9418482197136,
                                eta2,
                                fma(-33.67252878543393, eta1, 5.635069974807398),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_44_rdcp2` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_44_rdcp2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
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
        eta1 * s1,
        fma(
            0.019740635678180102
                * fma(
                    910729.0219043035,
                    eta5,
                    fma(
                        -696891.1349708183,
                        eta4,
                        fma(
                            201981.0743810629,
                            eta3,
                            fma(
                                -27212.67022616794,
                                eta2,
                                fma(1535.781942095697, eta1, -7.046295936301379),
                            ),
                        ),
                    ),
                ),
            s2,
            fma(
                -0.23311779185707152,
                fma(
                    -127148.27131388368,
                    eta5,
                    fma(
                        73397.72190288734,
                        eta4,
                        fma(
                            -14672.23221502075,
                            eta3,
                            fma(
                                1163.9107546486134,
                                eta2,
                                fma(-20.12558747513885, eta1, -1.0795711755430002),
                            ),
                        ),
                    ),
                ),
                0.025805905356653
                    * fma(
                        1.3238624538095295e6,
                        eta5,
                        fma(
                            -819148.9390278616,
                            eta4,
                            fma(
                                174164.91607515427,
                                eta3,
                                fma(
                                    -14580.02701600596,
                                    eta2,
                                    fma(350.93274421955806, eta1, 11.929946153728276),
                                ),
                            ),
                        ),
                    )
                    * s1,
            ),
        ),
        fma(
            eta1,
            fma(
                chidiff1 * delta,
                fma(
                    -123.94152825700999,
                    eta3,
                    fma(-2.861653255976984, eta1, 50.50227103211222 * eta2),
                ),
                chidiff2
                    * fma(
                        72.40230240887851,
                        eta3,
                        fma(2.9415751419018865, eta1, -(28.79779545444817 * eta2)),
                    ),
            ),
            eta1 * fma(
                49148.57065911245,
                eta5,
                fma(
                    -30595.4993786313,
                    eta4,
                    fma(
                        7168.843978909433,
                        eta3,
                        fma(
                            -792.0167314124681,
                            eta2,
                            fma(25.15310593958783, eta1, 3.2461722686239307),
                        ),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_44_rdcp3` (version 122022): ringdown amplitude
/// collocation point.
pub(crate) fn amp_44_rdcp3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        eta1 * s1,
        fma(
            -0.048465588779596405,
            fma(
                -133340.93723954144,
                eta5,
                fma(
                    76620.18243090731,
                    eta4,
                    fma(
                        -16522.343895064576,
                        eta3,
                        fma(
                            1744.329802302927,
                            eta2,
                            fma(-81.33156665674845, eta1, -0.34041762314288154),
                        ),
                    ),
                ),
            ),
            0.024804027856323612
                * fma(
                    584744.0406581408,
                    eta5,
                    fma(
                        -422282.0368440555,
                        eta4,
                        fma(
                            112832.04975245205,
                            eta3,
                            fma(
                                -13644.988225595187,
                                eta2,
                                fma(711.8727878341302, eta1, -8.666095805675418),
                            ),
                        ),
                    ),
                )
                * s1,
        ),
        fma(
            eta1,
            fma(
                chidiff1 * delta,
                fma(
                    66.58241012629095,
                    eta3,
                    fma(2.4286414692113816, eta1, -(23.213332913737403 * eta2)),
                ),
                chidiff2
                    * fma(
                        78.49621016381445,
                        eta3,
                        fma(3.085167288859442, eta1, -(31.60440418701438 * eta2)),
                    ),
            ),
            eta1 * fma(
                18536.937955014455,
                eta5,
                fma(
                    -12028.786386004691,
                    eta4,
                    fma(
                        2932.3415281149432,
                        eta3,
                        fma(
                            -337.70598252897696,
                            eta2,
                            fma(13.695204704208976, eta1, 0.861883217178703),
                        ),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_32_alambda` (version 122022): spheroidal ringdown
/// amplitude coefficient.
pub(crate) fn amp_32_alambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
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
        fma(
            14.163264397061505,
            eta5,
            fma(
                -15.707327897569396,
                eta4,
                fma(
                    6.5972460654253515,
                    eta3,
                    fma(0.08408469319155859, eta1, -(1.223794846617597 * eta2)),
                ),
            ),
        ),
        1.0 / fma(18.93655612952139, eta2, fma(-8.612447115134758, eta1, 1.0)),
        fma(
            s1,
            fma(
                -(0.00043088431510840695
                    * fma(
                        -6.094089675611567e6,
                        eta6,
                        fma(
                            3.9072674134789007e6,
                            eta5,
                            fma(
                                -870035.1394696251,
                                eta4,
                                fma(
                                    76517.21042363928,
                                    eta3,
                                    fma(-30.014415072587354, eta1, -(1900.5495690280086 * eta2)),
                                ),
                            ),
                        ),
                    )),
                s2,
                fma(
                    -0.0013330716557843666,
                    fma(
                        -150504.2927996263,
                        eta6,
                        fma(
                            252228.94931931415,
                            eta5,
                            fma(
                                -102969.88022112886,
                                eta4,
                                fma(
                                    16509.69605686402,
                                    eta3,
                                    fma(56.35538385647113, eta1, -(1218.1550992423377 * eta2)),
                                ),
                            ),
                        ),
                    ),
                    0.0010126460331462495
                        * fma(
                            -597155.7619486618,
                            eta6,
                            fma(
                                295469.28617550555,
                                eta5,
                                fma(
                                    -36980.079328277505,
                                    eta4,
                                    fma(
                                        -1304.9210590188136,
                                        eta3,
                                        fma(-33.87083889060834, eta1, 502.6221651850776 * eta2),
                                    ),
                                ),
                            ),
                        )
                        * s1,
                ),
            ),
            fma(
                chidiff1
                    * delta
                    * fma(
                        117.66195692191727,
                        eta5,
                        fma(5.786292428422767, eta3, -(53.60467819078566 * eta4)),
                    ),
                s1,
                fma(
                    chidiff2
                        * fma(
                            -42.74689906831975,
                            eta5,
                            fma(-0.4600401291210382, eta3, 12.23450117663151 * eta4),
                        ),
                    s1,
                    fma(
                        chidiff2,
                        fma(
                            -85.19723511005235,
                            eta5,
                            fma(-3.4614418482110163, eta3, 35.464117772624164 * eta4),
                        ),
                        chidiff1
                            * delta
                            * fma(
                                170.9266105597438,
                                eta5,
                                fma(2.0328561081997463, eta3, -(46.18751757691501 * eta4)),
                            ),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_32_lambda` (version 122022): spheroidal ringdown
/// amplitude coefficient.
pub(crate) fn amp_32_lambda(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1,
        fma(
            0.22358643406992756,
            fma(
                -915.2489655397206,
                eta3,
                fma(
                    365.2485484044132,
                    eta2,
                    fma(-32.06682257944444, eta1, 0.48943645614341924),
                ),
            ),
            0.0792473022309144
                * fma(
                    -3206.340850767219,
                    eta3,
                    fma(
                        1202.174780792418,
                        eta2,
                        fma(-103.65639889587327, eta1, 1.877251717679991),
                    ),
                )
                * s1,
        ),
        fma(
            chidiff2,
            fma(
                2657.6998831179444,
                eta5,
                fma(115.3689949926392, eta3, -(1159.8596972989067 * eta4)),
            ),
            fma(
                chidiff1 * delta,
                fma(
                    -6328.182823770599,
                    eta5,
                    fma(-188.37119473865533, eta3, 2151.8731700399308 * eta4),
                ),
                fma(
                    49.19719473681921,
                    eta3,
                    fma(
                        -12.259851752618998,
                        eta2,
                        fma(0.36457571743142897, eta1, 0.978510781593996),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_32_sigma` (version 122022): spheroidal ringdown
/// amplitude coefficient.
pub(crate) fn amp_32_sigma(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let eta5 = eta1 * eta4;
    let s1 = s;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1,
        fma(
            -0.012086025709597246,
            fma(
                17277.60594350428,
                eta3,
                fma(
                    -6606.1009717965735,
                    eta2,
                    fma(600.5968613752918, eta1, -6.230497473791485),
                ),
            ),
            -(0.06066548829900489
                * fma(
                    4119.373703246675,
                    eta3,
                    fma(
                        -1567.249168668069,
                        eta2,
                        fma(142.0346574366267, eta1, -0.9208054306316676),
                    ),
                )
                * s1),
        ),
        fma(
            chidiff2,
            fma(
                3592.287297392387,
                eta5,
                fma(209.09202210427972, eta3, -(1769.4658099037918 * eta4)),
            ),
            fma(
                chidiff1 * delta,
                fma(
                    123.86194078913776,
                    eta5,
                    fma(144.37065005786636, eta3, -(754.4085447486738 * eta4)),
                ),
                fma(0.13401718687342024, eta1, 1.3353917551819414),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_32_rdaux1` (version 122022): collocation point of
/// the auxiliary ringdown region.
pub(crate) fn amp_32_rdaux1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
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
        fma(
            2048.333892310575,
            eta5,
            fma(
                -1642.3057499049708,
                eta4,
                fma(
                    522.5070635563092,
                    eta3,
                    fma(4.45156488896258, eta1, -(77.39303992494544 * eta2)),
                ),
            ),
        ),
        1.0 / fma(24.249594730050312, eta2, fma(-9.611489164758915, eta1, 1.0)),
        fma(
            s1,
            fma(
                -0.022968727462555794,
                fma(
                    -3.572465179268999e6,
                    eta6,
                    fma(
                        2.6934917075009225e6,
                        eta5,
                        fma(
                            -710313.3016672594,
                            eta4,
                            fma(
                                77424.12614733395,
                                eta3,
                                fma(83.66854837403105, eta1, -(3330.6261333413177 * eta2)),
                            ),
                        ),
                    ),
                ),
                0.0014795114305436387
                    * fma(
                        7.06066766311127e7,
                        eta6,
                        fma(
                            -5.116110998398143e7,
                            eta5,
                            fma(
                                1.3705532554135624e7,
                                eta4,
                                fma(
                                    -1.6690169155105734e6,
                                    eta3,
                                    fma(-1672.7273629876313, eta1, 90877.38260964208 * eta2),
                                ),
                            ),
                        ),
                    )
                    * s1,
            ),
            fma(
                chidiff2,
                fma(
                    -131.19660856923554,
                    eta4,
                    fma(-4.188795724777721, eta2, 53.39200466700963 * eta3),
                ),
                chidiff1
                    * delta
                    * fma(
                        1242.865584938088,
                        eta4,
                        fma(14.284921364132623, eta2, -(321.26423637658746 * eta3)),
                    ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Amp_32_rdaux2` (version 122022): collocation point of
/// the auxiliary ringdown region.
pub(crate) fn amp_32_rdaux2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
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
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        fma(
            636.0141528226893,
            eta5,
            fma(
                -693.6803787318279,
                eta4,
                fma(
                    294.9496843041973,
                    eta3,
                    fma(3.902154247490771, eta1, -(55.77521071924907 * eta2)),
                ),
            ),
        ),
        1.0 / fma(19.119341007236955, eta2, fma(-8.56699762573719, eta1, 1.0)),
        fma(
            s1,
            fma(
                -0.16707403272774676,
                fma(
                    -1.1006903622406831e7,
                    eta7,
                    fma(
                        9.51553823212259e6,
                        eta6,
                        fma(
                            -3.1894624909922685e6,
                            eta5,
                            fma(
                                520786.0225074669,
                                eta4,
                                fma(
                                    -41908.45179140144,
                                    eta3,
                                    fma(6.678916447469937, eta1, 1331.480396625797 * eta2),
                                ),
                            ),
                        ),
                    ),
                ),
                0.015205286051218441
                    * fma(
                        1.1201307979786257e8,
                        eta7,
                        fma(
                            -9.865815275452062e7,
                            eta6,
                            fma(
                                3.379925277713386e7,
                                eta5,
                                fma(
                                    -5.635028227588545e6,
                                    eta4,
                                    fma(
                                        462957.5593513407,
                                        eta3,
                                        fma(108.10032279461095, eta1, -(16084.215590200103 * eta2)),
                                    ),
                                ),
                            ),
                        ),
                    )
                    * s1,
            ),
            fma(
                chidiff2,
                fma(
                    -440.26516625611,
                    eta4,
                    fma(-18.550171209458394, eta2, 188.99161055445936 * eta3),
                ),
                chidiff1
                    * delta
                    * fma(
                        1327.1224176812448,
                        eta4,
                        fma(13.132625215315063, eta2, -(340.5204040505528 * eta3)),
                    ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_RD_Phase_22_alpha2` (version 122019).
pub(crate) fn phase_22_alpha2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let no_spin = fma(
        55.45508989446867,
        eta4,
        fma(
            -31.330215053905395,
            eta3,
            fma(
                6.510807976353186,
                eta2,
                fma(-0.37138987533788487, eta, 0.2088669311744758),
            ),
        ),
    );
    let eq_spin = fma(
        38.61300158832203,
        eta3,
        fma(
            -16.874355161681766,
            eta2,
            fma(1.6966740823756759, eta, 0.2393965714370633),
        ),
    ) * s
        / fma(-0.633218538432246, s, 1.0);
    let uneq_spin = wf.dchi
        * fma(
            0.9088578269496244,
            eta.powf(2.5),
            15.619592332008951 * wf.dchi * eta.powf(3.5),
        )
        * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Phase_22_alphaL` (version 122019).
pub(crate) fn phase_22_alpha_l(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let s2 = s * s;
    let no_spin = eta
        * fma(
            27.623649807617376,
            eta3,
            fma(
                -16.504334734464244,
                eta2,
                fma(2.5400257699690143, eta, -1.1926122248825484),
            ),
        );
    let eq_spin = fma(
        eta4 * s,
        fma(
            180.8345521274853,
            s2,
            fma(-39.83460092417137, s, -72.50310678862684),
        ),
        fma(
            eta2 * s,
            fma(
                8.170616575493503,
                s2,
                fma(0.8743576195364632, s, -4.7282481007397825),
            ),
            fma(
                eta3 * s,
                fma(
                    -77.2346297158916,
                    s2,
                    fma(9.700178927988006, s, 35.803988443700824),
                ),
                eta * s
                    * fma(
                        -0.06417449517826644,
                        s2,
                        fma(-0.21477847929548569, s, 0.1034526554654983),
                    ),
            ),
        ),
    );
    let uneq_spin = fma(
        17.588573345324154 * (wf.chi2l * wf.chi2l),
        eta.powf(4.5),
        fma(
            -(35.17714669064831 * wf.chi1l * wf.chi2l),
            eta.powf(4.5),
            fma(
                17.588573345324154 * (wf.chi1l * wf.chi1l),
                eta.powf(4.5),
                fma(
                    -0.7428134042821221 * wf.chi1l,
                    eta.powf(3.5),
                    0.7428134042821221 * wf.chi2l * eta.powf(3.5),
                ),
            ),
        ),
    ) * delta;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift` (version 122019): fit of
/// `dphiS(fref) - dphi22(fref)`.
pub(crate) fn phase_32_spheroidal_time_shift(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        739270.8814129328,
        eta5,
        fma(
            -364129.24735853914,
            eta4,
            fma(
                61559.132976189896,
                eta3,
                fma(
                    -4565.033758777737,
                    eta2,
                    fma(167.95086712701223, eta, 11.851438981981772),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            -65993.84497923202,
            eta3,
            fma(
                24793.809334683185,
                eta2,
                fma(-2292.778112523392, eta, 15.585154698977842),
            ),
        ),
        s4,
        fma(
            fma(
                -60929.87465551446,
                eta3,
                fma(
                    24145.932943269272,
                    eta2,
                    fma(-2244.4613237912527, eta, 3.4881615575084797),
                ),
            ),
            s3,
            fma(
                fma(
                    26929.677144312944,
                    eta3,
                    fma(
                        -8046.364492927503,
                        eta2,
                        fma(434.31707030999445, eta, 9.506768471271634),
                    ),
                ),
                s,
                fma(
                    3575.347142399199,
                    eta3,
                    fma(
                        1334.1062451631644,
                        eta2,
                        fma(-307.67253970367034, eta, -5.949655484033632),
                    ),
                ) * s2,
            ),
        ),
    );
    let uneq_spin = 465.7904934097202 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift` (version 122019): fit of
/// `phiS(fref) - phi22(fref)` modulo 2 pi.
pub(crate) fn phase_32_spheroidal_phase_shift(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let eta7 = eta.powf(7.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -7.612233821752137e6,
        eta7,
        fma(
            5.614617173188322e6,
            eta6,
            fma(
                -1.6312968467540336e6,
                eta5,
                fma(
                    246313.12887984765,
                    eta4,
                    fma(
                        -21256.376324666326,
                        eta3,
                        fma(
                            1056.2426481245027,
                            eta2,
                            fma(-22.209549522908667, eta, -1.3328895897490733),
                        ),
                    ),
                ),
            ),
        ),
    );
    let eq_spin = s * fma(
        eta3,
        fma(
            805.0266029288334,
            s4,
            fma(
                -1493.177074873678,
                s3,
                fma(
                    -2683.947280170815,
                    s2,
                    fma(-22.59781899512552, s, 62.632493533037625),
                ),
            ),
        ),
        fma(
            eta,
            fma(
                24.02858084544326,
                s4,
                fma(
                    -34.498125910065156,
                    s3,
                    fma(
                        -97.16014978329918,
                        s2,
                        fma(-35.074072557604445, s, -47.475568056234245),
                    ),
                ),
            ),
            fma(
                0.7638287749489343,
                s4,
                fma(
                    eta2,
                    fma(
                        -312.1010566468068,
                        s4,
                        fma(
                            546.6246807461155,
                            s3,
                            fma(
                                1086.7504825459278,
                                s2,
                                fma(135.25216875394943, s, 196.7055281997748),
                            ),
                        ),
                    ),
                    fma(
                        -1.9586085340429995,
                        s3,
                        fma(
                            -1.1239505323267036,
                            s2,
                            fma(0.9960210841611344, s, -1.622727240110213),
                        ),
                    ),
                ),
            ),
        ),
    ) / fma(1.0, s, -2.950271397057221);
    let uneq_spin = fma(-4.0, eta, 1.0).sqrt()
        * fma(
            wf.chi1l * eta.powf(2.5),
            fma(30.01812659282717, s, -88.56162028006072),
            fma(
                wf.chi1l * eta2,
                fma(14.617728550838805, s, -43.126266433486435),
                fma(
                    wf.chi2l * eta.powf(2.5),
                    fma(-30.01812659282717, s, 88.56162028006072),
                    wf.chi2l * eta2 * fma(-14.617728550838805, s, 43.126266433486435),
                ),
            ),
        )
        / fma(1.0, s, -2.950271397057221);
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Phase_32_p1` (version 122019): spheroidal ringdown
/// phase-derivative collocation point.
pub(crate) fn phase_32_p1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let s5 = s.powf(5.0);
    let no_spin = fma(
        1.5674053633767858e6,
        eta5,
        fma(
            -817182.2896823225,
            eta4,
            fma(
                149846.7281073725,
                eta3,
                fma(
                    -12569.748101922158,
                    eta2,
                    fma(426.8372805022653, eta, 3169.372056189274),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            -384336.57477596396,
            eta3,
            fma(
                162941.8009556697,
                eta2,
                fma(-16556.65370439302, eta, -61.01931672442576),
            ),
        ),
        s5,
        fma(
            fma(
                -139322.02687945773,
                eta3,
                fma(
                    53893.9479532688,
                    eta2,
                    fma(-4901.068176970293, eta, -24.72321125342808),
                ),
            ),
            s4,
            fma(
                fma(
                    201204.5096556517,
                    eta3,
                    fma(
                        -91664.12621864353,
                        eta2,
                        fma(10145.53174542332, eta, 13.918237996338371),
                    ),
                ),
                s3,
                fma(
                    fma(
                        -3785.49764771212,
                        eta3,
                        fma(
                            7855.316419853637,
                            eta2,
                            fma(-1762.6573670619173, eta, 19.23408352151287),
                        ),
                    ),
                    s,
                    fma(
                        20497.5021807654,
                        eta3,
                        fma(
                            -5615.908682338113,
                            eta2,
                            fma(336.8340966473415, eta, -42.88446003698396),
                        ),
                    ) * s2,
                ),
            ),
        ),
    );
    let uneq_spin = wf.dchi
        * fma(-4.0, eta, 1.0).sqrt()
        * eta2
        * fma(
            13275.623692212472 * eta,
            s,
            fma(
                1600.240100295189 * wf.chi2l,
                eta,
                fma(-(1600.240100295189 * wf.chi1l), eta, 641.2473192044652),
            ),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Phase_32_p2` (version 122019): spheroidal ringdown
/// phase-derivative collocation point.
pub(crate) fn phase_32_p2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        7475.062269742079,
        eta3,
        fma(
            -2636.4344627081873,
            eta2,
            fma(206.09687819102305, eta, 3131.0260952676376),
        ),
    );
    let eq_spin = fma(
        fma(
            69666.01833764123,
            eta3,
            fma(
                -34444.2579678986,
                eta2,
                fma(5119.138772157134, eta, -227.8023564332453),
            ),
        ),
        s4,
        fma(
            fma(
                -27654.067482558712,
                eta3,
                fma(
                    7269.262546204149,
                    eta2,
                    fma(429.0787542102207, eta, -132.49105074500454),
                ),
            ),
            s3,
            fma(
                fma(
                    10514.68111669422,
                    eta3,
                    fma(
                        -434.60154548208334,
                        eta2,
                        fma(-691.9815135740145, eta, 49.90874152040307),
                    ),
                ),
                s,
                fma(
                    -56142.13736008524,
                    eta3,
                    fma(
                        26748.805404989867,
                        eta2,
                        fma(-3458.2579971189534, eta, 97.3078084654917),
                    ),
                ) * s2,
            ),
        ),
    );
    let uneq_spin = 477.51566939885424 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Phase_32_p3` (version 122019): spheroidal ringdown
/// phase-derivative collocation point.
pub(crate) fn phase_32_p3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        977.6115755788503,
        eta3,
        fma(
            -586.2469821978381,
            eta2,
            fma(76.94679795837645, eta, 3082.803556599222),
        ),
    );
    let eq_spin = fma(
        fma(
            69666.01833764123,
            eta3,
            fma(
                -34444.2579678986,
                eta2,
                fma(5119.138772157134, eta, -227.8023564332453),
            ),
        ),
        s4,
        fma(
            fma(
                16259.911572457446,
                eta3,
                fma(
                    -10971.864789147161,
                    eta2,
                    fma(2268.4693587493093, eta, -144.96370424517866),
                ),
            ),
            s3,
            fma(
                fma(
                    2472.6476419567534,
                    eta3,
                    fma(
                        1775.4343704616288,
                        eta2,
                        fma(-807.7353772747749, eta, 45.08944710349874),
                    ),
                ),
                s,
                fma(
                    -25583.314298758105,
                    eta3,
                    fma(
                        13821.251641893134,
                        eta2,
                        fma(-2224.9613131172046, eta, 95.57355060136699),
                    ),
                ) * s2,
            ),
        ),
    );
    let uneq_spin = 378.2359918274837 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_RD_Phase_32_p4` (version 122019): spheroidal ringdown
/// phase-derivative collocation point.
pub(crate) fn phase_32_p4(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -357.38692756785395,
        eta2,
        fma(64.99844502520415, eta, 3077.0657367004565),
    );
    let eq_spin = fma(
        fma(
            -41414.05689348732,
            eta2,
            fma(
                95314.99988114933,
                eta.powf(3.0),
                fma(4900.706039920717, eta, -169.30909412804587),
            ),
        ),
        s4,
        fma(
            fma(
                -26766.730897942085,
                eta2,
                fma(
                    62863.79877094988,
                    eta.powf(3.0),
                    fma(3169.3397351803583, eta, -126.02362949830213),
                ),
            ),
            s3,
            fma(
                fma(
                    5700.682624203565,
                    eta2,
                    fma(
                        -9490.641676924794,
                        eta.powf(3.0),
                        fma(-986.7751755509875, eta, 34.793450080444714),
                    ),
                ),
                s,
                fma(
                    11008.881935880598,
                    eta2,
                    fma(
                        -19906.416384606226,
                        eta.powf(3.0),
                        fma(-1644.6690499868596, eta, 57.38106384558743),
                    ),
                ) * s2,
            ),
        ),
    );
    let uneq_spin = 390.5443469721231 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/* End of the parameter-space fits. */

/// Ringdown amplitude collocation fits `rdcp1..rdcp3` (21, 33, 44).
fn amp_rdcp_fit(mode_int: usize, k: usize, wf: &WaveformStruct) -> f64 {
    const FITS: [[Fit; 3]; 4] = [
        [amp_21_rdcp1, amp_21_rdcp2, amp_21_rdcp3],
        [amp_33_rdcp1, amp_33_rdcp2, amp_33_rdcp3],
        [amp_44_rdcp1, amp_44_rdcp2, amp_44_rdcp3], // (3,2) uses the coefficient fits
        [amp_44_rdcp1, amp_44_rdcp2, amp_44_rdcp3],
    ];
    FITS[mode_int][k](wf)
}

/// `pPhase->RingdownPhaseFits`: the spheroidal phase-derivative collocation
/// fits `p1..p4` of the (3,2) mode.
pub(crate) fn phase_32_fit(i: usize, wf: &WaveformStruct) -> f64 {
    const FITS: [Fit; 4] = [phase_32_p1, phase_32_p2, phase_32_p3, phase_32_p4];
    FITS[i](wf)
}

/// `IMRPhenomXHM_RD_Amp_Coefficients` (ringdown amplitude versions 1 and 2)
/// with the exponential fall-off and, for the (3,2) mode,
/// `IMRPhenomXHM_RDAux_Amp_Coefficients`.
pub(crate) fn rd_amp_coefficients(
    amp: &mut HMAmpCoefficients,
    wfhm: &HMWaveformStruct,
    wf: &WaveformStruct,
) {
    let fdamp = wfhm.f_damp;
    if wfhm.ringdown_amp_version == 1 {
        amp.rd_coefficients[0] = amp_32_alambda(wf).abs();
        amp.rd_coefficients[1] = amp_32_lambda(wf);
        amp.rd_coefficients[2] = amp_32_sigma(wf);
    } else {
        let rdcp1 = amp_rdcp_fit(wfhm.mode_int, 0, wf).abs();
        let rdcp2 = amp_rdcp_fit(wfhm.mode_int, 1, wf).abs();
        let mut rdcp3 = amp_rdcp_fit(wfhm.mode_int, 2, wf).abs();
        amp.colloc_freqs_rd = [wfhm.f_ring - fdamp, wfhm.f_ring, wfhm.f_ring + fdamp];
        // Vetos
        if rdcp3 >= rdcp2 * rdcp2 / rdcp1 {
            rdcp3 = 0.5 * rdcp2 * rdcp2 / rdcp1;
        }
        if rdcp3 > rdcp2 {
            rdcp3 = 0.5 * rdcp2;
        }
        if rdcp1 < rdcp2 && rdcp3 > rdcp1 {
            rdcp3 = rdcp1;
        }
        amp.colloc_values_rd = [rdcp1, rdcp2, rdcp3];
        let mut deno = (rdcp1 / rdcp3).sqrt() - (rdcp1 / rdcp2);
        if deno <= 0. {
            deno = 1e-16;
        }
        amp.rd_coefficients[0] = rdcp1 * fdamp / deno;
        amp.rd_coefficients[2] = (amp.rd_coefficients[0] / (rdcp2 * fdamp)).sqrt();
        amp.rd_coefficients[1] = 0.5 * amp.rd_coefficients[2] * (rdcp1 / rdcp3).ln();
    }

    if wfhm.f_amp_rd_falloff > 0. {
        let powers = UsefulPowers::new(wfhm.f_amp_rd_falloff);
        let value = amp_ansatz_falloff(&powers, wfhm, amp, 0.);
        amp.rd_coefficients[3] = value;
        amp.rd_coefficients[4] = -1. * amp_dansatz(&powers, wfhm, amp) / value;
    }
    if amp.n_coefficients_rd_aux > 0 {
        rd_aux_amp_coefficients(amp, wfhm, wf);
    }
}

/// `IMRPhenomXHM_RDAux_Amp_Coefficients`: cubic polynomial joining the
/// intermediate region to the spheroidal ringdown of the (3,2) mode.
fn rd_aux_amp_coefficients(
    amp: &mut HMAmpCoefficients,
    wfhm: &HMWaveformStruct,
    wf: &WaveformStruct,
) {
    amp.colloc_values_rd_aux[0] = amp_32_rdaux1(wf).abs();
    amp.colloc_values_rd_aux[1] = amp_32_rdaux2(wf).abs();
    let powers = UsefulPowers::new(amp.f_rd_aux);
    amp.colloc_values_rd_aux[2] = amp_ansatz(&powers, wfhm, amp);
    amp.colloc_values_rd_aux[3] = amp_dansatz(&powers, wfhm, amp);
    amp.colloc_freqs_rd_aux = [
        amp.f_amp_match_im,
        0.5 * (amp.f_amp_match_im + amp.f_rd_aux),
        amp.f_rd_aux,
        amp.f_rd_aux,
    ];
    let mut matrix = [[0.; 4]; 4];
    for (i, row) in matrix.iter_mut().enumerate() {
        let fcoll = amp.colloc_freqs_rd_aux[i];
        let mut fpower = 1.;
        if i < 3 {
            for entry in row.iter_mut() {
                *entry = fpower;
                fpower *= fcoll;
            }
        } else {
            // Derivative at the right boundary.
            row[0] = 0.;
            for (j, entry) in row.iter_mut().enumerate().skip(1) {
                *entry = j as f64 * fpower;
                fpower *= fcoll;
            }
        }
    }
    amp.rd_aux_coefficients = lu_solve(matrix, amp.colloc_values_rd_aux);
}

/// `IMRPhenomXHM_RD_Amp_Ansatz` (versions 1 and 2): Lorentzian with an
/// exponential decay, the polynomial of the auxiliary region of the (3,2)
/// below `fRDAux`, and the exponential fall-off above `fAmpRDfalloff`. For
/// the (3,2) this is the spheroidal amplitude. Strain units without `amp0`.
#[inline]
pub(crate) fn amp_ansatz(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    amp_ansatz_falloff(powers, wfhm, amp, wfhm.f_amp_rd_falloff)
}

#[inline]
fn amp_ansatz_falloff(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
    f_amp_rd_falloff: f64,
) -> f64 {
    let ff = powers.itself;
    let rd = &amp.rd_coefficients;
    if amp.n_coefficients_rd_aux > 0 && ff < amp.f_rd_aux {
        let mut amp_rd = 0.;
        let mut fpower = 1.;
        for &coefficient in &amp.rd_aux_coefficients[..amp.n_coefficients_rd_aux] {
            amp_rd = fma(fpower, coefficient, amp_rd);
            fpower *= ff;
        }
        amp_rd
    } else if f_amp_rd_falloff > 0. && ff >= f_amp_rd_falloff {
        rd[3] * (-rd[4] * (ff - f_amp_rd_falloff)).exp()
    } else {
        let fda = wfhm.f_damp;
        let dfr = ff - wfhm.f_ring;
        let dfd = fda * rd[2];
        rd[0] * fda / ((rd[1] / dfd * dfr).exp() * fma(dfr, dfr, dfd * dfd))
    }
}

/// [`amp_ansatz`] with its first and second derivatives.
pub(crate) fn amp_ansatz_jet(ff: f64, wfhm: &HMWaveformStruct, amp: &HMAmpCoefficients) -> Jet {
    let x = Jet::variable(ff);
    let rd = &amp.rd_coefficients;
    let f_amp_rd_falloff = wfhm.f_amp_rd_falloff;
    if amp.n_coefficients_rd_aux > 0 && ff < amp.f_rd_aux {
        let mut amp_rd = Jet::constant(0.);
        let mut fpower = Jet::constant(1.);
        for &coefficient in &amp.rd_aux_coefficients[..amp.n_coefficients_rd_aux] {
            amp_rd = amp_rd + fpower * coefficient;
            fpower = fpower * x;
        }
        amp_rd
    } else if f_amp_rd_falloff > 0. && ff >= f_amp_rd_falloff {
        ((x - f_amp_rd_falloff) * -rd[4]).exp() * rd[3]
    } else {
        let fda = wfhm.f_damp;
        let dfr = x - wfhm.f_ring;
        let dfd = fda * rd[2];
        rd[0] * fda / ((dfr * (rd[1] / dfd)).exp() * (dfr * dfr + dfd * dfd))
    }
}

/// `IMRPhenomXHM_RD_Amp_DAnsatz` (versions 1 and 2): derivative of the
/// Lorentzian ringdown amplitude.
pub(crate) fn amp_dansatz(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    let ff = powers.itself;
    let fda = wfhm.f_damp;
    let rd = &amp.rd_coefficients;
    let dfr = ff - wfhm.f_ring;
    let numerator = rd[0]
        * fma(
            fda * fda * rd[1] * rd[2],
            rd[2],
            fma(dfr * dfr, rd[1], 2. * fda * dfr * rd[2]),
        );
    let mut denom = fma(dfr, dfr, fda * fda * rd[2] * rd[2]);
    denom = rd[2] * denom * denom * (dfr * rd[1] / (fda * rd[2])).exp();
    -numerator / denom
}

/// `IMRPhenomXHM_RD_Phase_Ansatz`: the ringdown phase derivative, the
/// rescaled 22 ringdown for (2,1), (3,3), (4,4) and the spheroidal ansatz
/// (version 122019) for the (3,2).
pub(crate) fn phase_ansatz(
    ff: f64,
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    phase: &HMPhaseCoefficients,
) -> f64 {
    let frd = wfhm.f_ring;
    let fda = wfhm.f_damp;
    let den = fma(fda, fda, (ff - frd) * (ff - frd));
    if !wfhm.mixing_on {
        fma(frd * frd * phase.alpha2, powers.m_two, phase.alpha0) + (phase.alpha_l * fda) / den
    } else {
        let d = fma(phase.alpha2_s, powers.m_two, phase.alpha0_s);
        fma(phase.alpha4_s, powers.m_four, d) + (phase.alpha_l_s * fda) / den
    }
}

/// [`phase_ansatz`] (the phase derivative) with its first and second
/// derivatives.
pub(crate) fn phase_ansatz_jet(
    ff: f64,
    wfhm: &HMWaveformStruct,
    phase: &HMPhaseCoefficients,
) -> Jet {
    let x = Jet::variable(ff);
    let frd = wfhm.f_ring;
    let fda = wfhm.f_damp;
    let lorentzian = fda / ((x - frd) * (x - frd) + fda * fda);
    if !wfhm.mixing_on {
        Jet::power(ff, -2.) * (frd * frd * phase.alpha2) + lorentzian * phase.alpha_l + phase.alpha0
    } else {
        Jet::power(ff, -2.) * phase.alpha2_s
            + Jet::power(ff, -4.) * phase.alpha4_s
            + lorentzian * phase.alpha_l_s
            + phase.alpha0_s
    }
}

/// `IMRPhenomXHM_RD_Phase_AnsatzInt`: the ringdown phase, see
/// [`phase_ansatz`].
#[inline]
pub(crate) fn phase_ansatz_int(
    ff: f64,
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    phase: &HMPhaseCoefficients,
) -> f64 {
    let invf = powers.m_one;
    let frd = wfhm.f_ring;
    let fda = wfhm.f_damp;
    let arctan = ((ff - frd) / fda).atan();
    if !wfhm.mixing_on {
        let phase_rd = fma(phase.alpha0, ff, -(frd * frd * phase.alpha2 * invf));
        fma(phase.alpha_l, arctan, phase_rd)
    } else {
        let mut phase_rd = fma(phase.alpha0_s, ff, phase.phi0_s);
        phase_rd = fma(-phase.alpha2_s, invf, phase_rd);
        phase_rd = fma(-(1. / 3. * phase.alpha4_s), powers.m_three, phase_rd);
        fma(phase.alpha_l_s, arctan, phase_rd)
    }
}
