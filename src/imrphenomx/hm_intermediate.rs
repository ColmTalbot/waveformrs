// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Intermediate region of the IMRPhenomXHM higher modes,
//! `LALSimIMRPhenomXHM_intermediate.c`: parameter-space fits of the
//! amplitude (version 122022) and phase-derivative (version 122019)
//! collocation points, the reconstruction of the intermediate amplitude from
//! equispaced collocation points, and the phase ansatzes.

use num_complex::Complex;

use crate::imrphenomx::hm_internals::{
    Fit, HMAmpCoefficients, HMPhaseCoefficients, HMWaveformStruct,
};
use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct};
use crate::imrphenomx::jet;
use crate::imrphenomx::linalg::lu_solve;
use crate::imrphenomx::{hm_inspiral, hm_ringdown};

/// `IMRPhenomXHM_Inter_Amp_21_int1` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_21_int1(wf: &WaveformStruct) -> f64 {
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
        delta * eta1 * s1,
        fma(
            -(0.46006975902558517
                * fma(
                    731179.456702448,
                    eta5,
                    fma(
                        -528521.6231015554,
                        eta4,
                        fma(
                            142904.08962903373,
                            eta3,
                            fma(
                                -17333.549622759732,
                                eta2,
                                fma(826.026625945615, eta1, -2.1663474937625975),
                            ),
                        ),
                    ),
                )),
            s2,
            fma(
                -8.276366844994188,
                fma(
                    -23896.392706809565,
                    eta5,
                    fma(
                        16770.447637953577,
                        eta4,
                        fma(
                            -4389.799658723288,
                            eta3,
                            fma(
                                516.7886322104276,
                                eta2,
                                fma(-24.12941323757896, eta1, 1.0677538075697492),
                            ),
                        ),
                    ),
                ),
                -(1.6908277400304084
                    * fma(
                        -2085.0821513466053,
                        eta5,
                        fma(
                            592.300353344717,
                            eta4,
                            fma(
                                -184.13091281984674,
                                eta3,
                                fma(
                                    114.8330693231833,
                                    eta2,
                                    fma(-29.00026389706585, eta1, 3.4799140066657928),
                                ),
                            ),
                        ),
                    )
                    * s1),
            ),
        ),
        fma(
            delta * eta1,
            fma(
                chidiff1
                    * fma(
                        1022.8553774274128,
                        eta3,
                        fma(34.17829404207186, eta1, -(386.34587928670015 * eta2)),
                    ),
                s1,
                chidiff1
                    * fma(
                        1016.6019654342113,
                        eta3,
                        fma(56.76554600963724, eta1, -(491.4593694689354 * eta2)),
                    )
                    * s2,
            ),
            fma(
                chidiff1 * eta5,
                fma(
                    -1.3098861736238572,
                    s2,
                    fma(-195.27885562092342, s1, -1698.841763891122),
                ),
                fma(
                    delta * eta1,
                    fma(
                        chidiff2,
                        fma(
                            63.70715919820867,
                            eta3,
                            fma(5.159755997682368, eta1, -(30.293198248154948 * eta2)),
                        ),
                        chidiff1
                            * fma(
                                1427.5951158851076,
                                eta3,
                                fma(8.262642080222694, eta1, -(415.88826990259116 * eta2)),
                            ),
                    ),
                    delta
                        * eta1
                        * fma(
                            2199.614895727586,
                            eta4,
                            fma(
                                -1614.178472020212,
                                eta3,
                                fma(
                                    447.2214642597892,
                                    eta2,
                                    fma(-66.46950491124205, eta1, 18.55363583212328),
                                ),
                            ),
                        ),
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_Inter_Amp_21_int2` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_21_int2(wf: &WaveformStruct) -> f64 {
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
    fma(
        delta * eta1 * s1,
        fma(
            -(0.5006664599166409
                * fma(
                    21974.8201355744,
                    eta4,
                    fma(
                        -11939.621138801092,
                        eta3,
                        fma(
                            2334.0029399672358,
                            eta2,
                            fma(-212.00728173197498, eta1, 10.289020582277626),
                        ),
                    ),
                )),
            s2,
            fma(
                -8.097876227116853,
                fma(
                    -1063.0692862578455,
                    eta4,
                    fma(
                        593.8035334117192,
                        eta3,
                        fma(
                            -116.88834714736281,
                            eta2,
                            fma(9.861355377849485, eta1, 0.6569459700232806),
                        ),
                    ),
                ),
                -(1.0546375154878165
                    * fma(
                        -7141.611333893155,
                        eta4,
                        fma(
                            4350.442990924205,
                            eta3,
                            fma(
                                -902.5751736558435,
                                eta2,
                                fma(65.25215540635162, eta1, 0.745557030602097),
                            ),
                        ),
                    )
                    * s1),
            ),
        ),
        fma(
            delta * eta1,
            fma(
                chidiff1
                    * fma(
                        984.6826660818275,
                        eta3,
                        fma(36.358142200869295, eta1, -(384.2123173145321 * eta2)),
                    ),
                s1,
                chidiff1
                    * fma(
                        -639.190132707115,
                        eta3,
                        fma(4.159271594881928, eta1, 105.10911749116399 * eta2),
                    )
                    * s2,
            ),
            fma(
                chidiff1 * eta5,
                fma(
                    -3.3637011939285015,
                    s2,
                    fma(-102.34716990474509, s1, -1549.6199518612063),
                ),
                fma(
                    delta * eta1,
                    fma(
                        12.048194546899204,
                        eta2,
                        fma(-12.783698329428516, eta1, 13.757856231617446),
                    ),
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            1333.6181181686939,
                            eta3,
                            fma(15.107530092096438, eta1, -(416.811753638553 * eta2)),
                        ),
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_Inter_Amp_21_int3` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_21_int3(wf: &WaveformStruct) -> f64 {
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
    fma(
        delta * eta1 * s1,
        fma(
            -(0.4308267743835775
                * fma(
                    29480.158198408903,
                    eta4,
                    fma(
                        -17712.883990278668,
                        eta3,
                        fma(
                            3662.099161055873,
                            eta2,
                            fma(-302.9708323417439, eta1, 9.970654092010587),
                        ),
                    ),
                )),
            s2,
            fma(
                -8.665646058245033,
                fma(
                    -1001.2370065269745,
                    eta4,
                    fma(
                        576.7172598056907,
                        eta3,
                        fma(
                            -111.70764910503321,
                            eta2,
                            fma(8.293609541933655, eta1, 0.7862132291286934),
                        ),
                    ),
                ),
                -(0.9459820574514348
                    * fma(
                        -7518.309520232795,
                        eta4,
                        fma(
                            4331.56002883546,
                            eta3,
                            fma(
                                -817.7854010574645,
                                eta2,
                                fma(48.94077040282239, eta1, 1.309016452198605),
                            ),
                        ),
                    )
                    * s1),
            ),
        ),
        fma(
            chidiff1 * eta5,
            fma(
                3.0053086651515843,
                s2,
                fma(-16.639896279650678, s1, -1522.8543551416456),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        343.79587860999874,
                        eta3,
                        fma(11.877985158418596, eta1, -(131.04937626836355 * eta2)),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(
                        -127.62161159798488,
                        eta3,
                        fma(
                            76.14884211156267,
                            eta2,
                            fma(-21.755549987331054, eta1, 13.318990196097973),
                        ),
                    ),
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            1366.2408490833282,
                            eta3,
                            fma(17.704321326939414, eta1, -(434.4390350012534 * eta2)),
                        ),
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_Inter_Amp_21_int4` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_21_int4(wf: &WaveformStruct) -> f64 {
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
    fma(
        delta * eta1 * s1,
        fma(
            -(0.21540359093306097
                * fma(
                    -9835.416414460055,
                    eta4,
                    fma(
                        1889.9613367802453,
                        eta3,
                        fma(
                            385.8345793110262,
                            eta2,
                            fma(-109.06544597367301, eta1, 9.42805409480658),
                        ),
                    ),
                )),
            s2,
            fma(
                -9.499741513411829,
                fma(
                    -248.5050377498615,
                    eta4,
                    fma(
                        166.04254881175257,
                        eta3,
                        fma(
                            -33.651192908287236,
                            eta2,
                            fma(2.400945118514037, eta1, 0.912120958549489),
                        ),
                    ),
                ),
                -(0.7850652143322492
                    * fma(
                        -9077.037917192794,
                        eta4,
                        fma(
                            5381.481380750608,
                            eta3,
                            fma(
                                -1032.1319480683567,
                                eta2,
                                fma(60.81773903539479, eta1, 1.534131218043425),
                            ),
                        ),
                    )
                    * s1),
            ),
        ),
        fma(
            chidiff1 * eta5,
            fma(
                -2.2456928156392197,
                s2,
                fma(74.44462583487345, s1, -1528.956015503355),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        280.6003655502327,
                        eta3,
                        fma(12.925363020014743, eta1, -(126.87194512915104 * eta2)),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(
                        -139.25546924151664,
                        eta3,
                        fma(
                            83.20619262213437,
                            eta2,
                            fma(-22.831152256559523, eta1, 13.094382343446163),
                        ),
                    ),
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            1430.3698681181,
                            eta3,
                            fma(20.120192352555357, eta1, -(458.2592421214168 * eta2)),
                        ),
                ),
            ),
        ),
    )
    .abs()
}

/// `IMRPhenomXHM_Inter_Amp_33_int1` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_33_int1(wf: &WaveformStruct) -> f64 {
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
    fma(
        chidiff1 * eta5,
        fma(
            7.650565415855383,
            s2,
            fma(53.36849263931051, s1, 199.6863414922219),
        ),
        fma(
            delta * eta1 * s1,
            fma(
                -0.5383304368673182,
                fma(
                    11959.83533107835,
                    eta5,
                    fma(
                        -9675.194644814854,
                        eta4,
                        fma(
                            3507.3655567272644,
                            eta3,
                            fma(
                                -843.7897535238325,
                                eta2,
                                fma(129.36947401891433, eta1, -7.456619067234563),
                            ),
                        ),
                    ),
                ),
                -(0.28042799223829407
                    * fma(
                        166056.27237509796,
                        eta5,
                        fma(
                            -119209.70783201039,
                            eta4,
                            fma(
                                32634.43965039936,
                                eta3,
                                fma(
                                    -4241.537539226717,
                                    eta2,
                                    fma(266.69059813274475, eta1, -6.212827413930676),
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
                        99.91399033058927,
                        eta3,
                        fma(5.575274516197629, eta1, -(44.592719238427094 * eta2)),
                    ),
                s1,
                fma(
                    chidiff1 * delta * eta1,
                    fma(
                        -148.3162618111991,
                        eta3,
                        fma(-0.3516244197696068, eta1, 40.425151307421416 * eta2),
                    ),
                    delta
                        * eta1
                        * fma(
                            3996.417635728702,
                            eta4,
                            fma(
                                -3047.513201789169,
                                eta3,
                                fma(
                                    914.5350366065115,
                                    eta2,
                                    fma(-146.29035440932105, eta1, 26.998512565991778),
                                ),
                            ),
                        ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_33_int2` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_33_int2(wf: &WaveformStruct) -> f64 {
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
    fma(
        delta * eta1 * s1,
        fma(
            0.30983324991828787
                * fma(
                    -504089.07690692553,
                    eta5,
                    fma(
                        376409.2335857397,
                        eta4,
                        fma(
                            -104863.50453518033,
                            eta3,
                            fma(
                                13278.047199998868,
                                eta2,
                                fma(-719.9854052004307, eta1, 11.302200127272357),
                            ),
                        ),
                    ),
                ),
            s2,
            fma(
                0.3104472390387834,
                fma(
                    205525.13910847943,
                    eta5,
                    fma(
                        -154255.3408672655,
                        eta4,
                        fma(
                            42566.59565666228,
                            eta3,
                            fma(
                                -4964.199967099143,
                                eta2,
                                fma(169.85423386969634, eta1, 6.073591341439855),
                            ),
                        ),
                    ),
                ),
                0.2295327944679772
                    * fma(
                        79210.33893514868,
                        eta5,
                        fma(
                            -41567.07396803811,
                            eta4,
                            fma(
                                2404.4151687877525,
                                eta3,
                                fma(
                                    1876.408148917458,
                                    eta2,
                                    fma(-354.7914372697625, eta1, 19.236275867648594),
                                ),
                            ),
                        ),
                    )
                    * s1,
            ),
        ),
        fma(
            chidiff1 * eta5,
            fma(
                -10.05114600643587,
                s2,
                fma(51.2153901931183, s1, 234.13476431269862),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        60.38649425562178,
                        eta3,
                        fma(-0.3265505633310564, eta1, -(9.861644053348053 * eta2)),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(
                        50.726220750178435,
                        eta2,
                        fma(-28.970875603981295, eta1, 17.42562079069636),
                    ),
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            -273.1170921735085,
                            eta3,
                            fma(-7.861956897615623, eta1, 93.45476935080045 * eta2),
                        ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_33_int3` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_33_int3(wf: &WaveformStruct) -> f64 {
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
    fma(
        delta * eta1 * s1,
        fma(
            0.3695850566805098
                * fma(
                    -534437.8832860565,
                    eta5,
                    fma(
                        360950.6686625352,
                        eta4,
                        fma(
                            -91096.64982858274,
                            eta3,
                            fma(
                                10426.348944657859,
                                eta2,
                                fma(-498.95526910874986, eta1, 6.077183107132255),
                            ),
                        ),
                    ),
                ),
            s2,
            fma(
                0.4590550434774332,
                fma(
                    282064.174913521,
                    eta5,
                    fma(
                        -195187.55422847517,
                        eta4,
                        fma(
                            50003.94075934942,
                            eta3,
                            fma(
                                -5545.0086137386525,
                                eta2,
                                fma(207.74562213604057, eta1, 3.0594364612798635),
                            ),
                        ),
                    ),
                ),
                0.657748992123043
                    * fma(
                        -20879.999628681435,
                        eta5,
                        fma(
                            19714.675715229736,
                            eta4,
                            fma(
                                -6999.7659193505915,
                                eta3,
                                fma(
                                    1276.6209573025596,
                                    eta2,
                                    fma(-124.06189543062042, eta1, 5.57939137343977),
                                ),
                            ),
                        ),
                    )
                    * s1,
            ),
        ),
        fma(
            chidiff1 * eta5,
            fma(
                10.457035444964653,
                s2,
                fma(44.408527093031566, s1, 263.5151703818307),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        -73.62033985436068,
                        eta3,
                        fma(-7.048898856045782, eta1, 49.941617405768135 * eta2),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(
                        16.79500349318081,
                        eta2,
                        fma(-12.799844096694798, eta1, 14.555522136327964),
                    ),
                    chidiff1
                        * delta
                        * eta1
                        * fma(
                            -562.0936797781519,
                            eta3,
                            fma(-16.292654447108134, eta1, 190.3516012682791 * eta2),
                        ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_33_int4` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_33_int4(wf: &WaveformStruct) -> f64 {
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
        delta * eta1 * s1,
        fma(
            0.2181399561677432
                * fma(
                    277256.38289831823,
                    eta5,
                    fma(
                        -216137.90827404748,
                        eta4,
                        fma(
                            59288.71069769708,
                            eta3,
                            fma(
                                -6857.936471898544,
                                eta2,
                                fma(303.9394443302189, eta1, -1.6754158383043574),
                            ),
                        ),
                    ),
                ),
            s2,
            fma(
                0.6893003654021495,
                fma(
                    -27679.226582207506,
                    eta5,
                    fma(
                        19514.812242035154,
                        eta4,
                        fma(
                            -4797.429130246123,
                            eta3,
                            fma(
                                565.3767256471909,
                                eta2,
                                fma(-44.83989278653052, eta1, 3.1014226377197027),
                            ),
                        ),
                    ),
                ),
                0.7068016563068026
                    * fma(
                        -65658.74746265226,
                        eta5,
                        fma(
                            48603.96661003743,
                            eta4,
                            fma(
                                -13485.270489656365,
                                eta3,
                                fma(
                                    1788.1730303291356,
                                    eta2,
                                    fma(-118.51094098279343, eta1, 4.071212304920691),
                                ),
                            ),
                        ),
                    )
                    * s1,
            ),
        ),
        fma(
            chidiff1 * eta5,
            fma(
                -29.471275170013012,
                s2,
                fma(
                    7.635582066682054,
                    s1,
                    fma(-21.102665129433042, chidiff2, 287.90444670305715),
                ),
            ),
            fma(
                chidiff1
                    * delta
                    * eta1
                    * fma(
                        -1485.5134187612182,
                        eta3,
                        fma(-46.85669289852532, eta1, 522.3965959942979 * eta2),
                    ),
                s1,
                fma(
                    delta * eta1,
                    fma(
                        17.098576301150125,
                        eta2,
                        fma(-7.449975618083432, eta1, 13.312095699772305),
                    ),
                    delta
                        * eta1
                        * fma(
                            chidiff1,
                            fma(
                                -1103.1917047361735,
                                eta3,
                                fma(-31.171150896110156, eta1, 371.1389274783572 * eta2),
                            ),
                            chidiff2
                                * fma(
                                    1164.9282236341376,
                                    eta3,
                                    fma(32.78644599730888, eta1, -(395.15713118955387 * eta2)),
                                ),
                        ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_32_int1` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_32_int1(wf: &WaveformStruct) -> f64 {
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
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        fma(
            106683.98302194536,
            eta4,
            fma(
                -99873.4707358776,
                eta3,
                fma(480.0434256230109, eta1, 25346.341240810478 * eta2),
            ),
        ) * sqroot,
        1.0 / fma(10083.297670051445, eta2, fma(1082.6574834474493, eta1, 1.0)),
        fma(
            s1 * fma(
                0.7067243321652764,
                fma(
                    -1929.9084724384807,
                    eta5,
                    fma(
                        175.4194342269804,
                        eta4,
                        fma(
                            325.71940024255775,
                            eta3,
                            fma(
                                -85.46581740280585,
                                eta2,
                                fma(9.608999847549535, eta1, 0.12885110296881636),
                            ),
                        ),
                    ),
                ),
                0.1540566313813899
                    * fma(
                        43309.69769878964,
                        eta5,
                        fma(
                            -28804.241518798244,
                            eta4,
                            fma(
                                7184.647314370326,
                                eta3,
                                fma(
                                    -827.591235943271,
                                    eta2,
                                    fma(45.55785402900492, eta1, -0.3261041495083288),
                                ),
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
                        -8.68296921440857,
                        eta3,
                        fma(-0.2341404256829785, eta1, 2.606326837996192 * eta2),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            73.40268975811729,
                            eta3,
                            fma(0.5454562486736877, eta1, -(25.19759222940851 * eta2)),
                        ),
                ),
                sqroot,
                chidiff1
                    * delta
                    * fma(
                        32.22238925527844,
                        eta3,
                        fma(0.4422257616009941, eta1, -(8.490112284851655 * eta2)),
                    )
                    * s1
                    * sqroot,
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_32_int2` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_32_int2(wf: &WaveformStruct) -> f64 {
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
        eta1 * s1,
        fma(
            1.5098291294292217,
            fma(
                -1633.3295758142904,
                eta4,
                fma(
                    856.9917885742416,
                    eta3,
                    fma(
                        -143.66427232396376,
                        eta2,
                        fma(9.848766999273414, eta1, 0.4844667556328104),
                    ),
                ),
            ),
            0.32413108737204144
                * fma(
                    6660.304740652403,
                    eta4,
                    fma(
                        -3811.5047139343505,
                        eta3,
                        fma(
                            761.6103793011912,
                            eta2,
                            fma(-62.37317183581803, eta1, 2.835358206961064),
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
                    68.72798751937934,
                    eta3,
                    fma(1.2198776533959694, eta1, -(26.816651899746475 * eta2)),
                ),
            s1,
            fma(
                eta1,
                fma(
                    chidiff2,
                    fma(
                        -128.88334273588077,
                        eta3,
                        fma(-4.175680729484314, eta1, 47.54281549129226 * eta2),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            208.07105580635888,
                            eta3,
                            fma(-0.18274358639599947, eta1, -(71.01128541687838 * eta2)),
                        ),
                ),
                eta1 * fma(
                    77839.74180254337,
                    eta6,
                    fma(
                        -69430.10468748478,
                        eta5,
                        fma(
                            24956.9592553473,
                            eta4,
                            fma(
                                -4544.076411013166,
                                eta3,
                                fma(
                                    456.2188780552874,
                                    eta2,
                                    fma(-38.57900689641654, eta1, 4.760999387359598),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_32_int3` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_32_int3(wf: &WaveformStruct) -> f64 {
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
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        s1,
        fma(
            -(0.03985156529181297
                * fma(
                    162398.81483375572,
                    eta6,
                    fma(
                        -46454.57664737511,
                        eta5,
                        fma(
                            -7859.4784979916085,
                            eta4,
                            fma(-304.92981902871617, eta2, 3614.518459296278 * eta3),
                        ),
                    ),
                )),
            s2,
            fma(
                0.3650871458400108,
                fma(
                    16001.318492586077,
                    eta6,
                    fma(
                        -18580.859291282686,
                        eta5,
                        fma(
                            6734.058809060536,
                            eta4,
                            fma(71.57390929624825, eta2, -(994.5272351916166 * eta3)),
                        ),
                    ),
                ),
                0.0960146077440495
                    * fma(
                        451546.88775684836,
                        eta6,
                        fma(
                            -318877.43061174755,
                            eta5,
                            fma(
                                83403.5743434538,
                                eta4,
                                fma(451.74917589707513, eta2, -(9719.470997418284 * eta3)),
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
                    -26.46699529970884,
                    eta4,
                    fma(-0.2566958540737833, eta2, 14.488550203412675 * eta3),
                ),
            s1,
            fma(
                403.5748987560612,
                eta5,
                fma(
                    chidiff1 * delta,
                    fma(
                        387.00162600306226,
                        eta4,
                        fma(1.924723094787216, eta2, -(90.92716917757797 * eta3)),
                    ),
                    fma(
                        -58.408542930248046,
                        eta4,
                        fma(
                            chidiff2,
                            fma(
                                -203.65975422378446,
                                eta4,
                                fma(-7.065118823041031, eta2, 77.97950589523865 * eta3),
                            ),
                            fma(
                                1.7262466525848588,
                                eta3,
                                fma(3.881450518842405, eta1, -(12.580316392558837 * eta2)),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_32_int4` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_32_int4(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.s_tot_r;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let s1 = s;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        eta1 * s1,
        fma(
            1.242750265695504,
            fma(
                -3797.3561904195085,
                eta4,
                fma(
                    1858.1518762309654,
                    eta3,
                    fma(
                        -303.5223688400309,
                        eta2,
                        fma(25.168480118950065, eta1, -0.522172424518215),
                    ),
                ),
            ),
            0.2927045241764365
                * fma(
                    6097.887891566872,
                    eta4,
                    fma(
                        -3131.5783196211587,
                        eta3,
                        fma(
                            471.64047356915603,
                            eta2,
                            fma(-15.488754837330958, eta1, 0.5056957789079993),
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
                    -144.21651114700492,
                    eta3,
                    fma(-3.582040878719185, eta1, 57.75888914133383 * eta2),
                ),
            s1,
            fma(
                eta1,
                fma(
                    chidiff2,
                    fma(
                        -236.2438921965621,
                        eta3,
                        fma(-8.572797326909152, eta1, 92.95723645687826 * eta2),
                    ),
                    chidiff1
                        * delta
                        * fma(
                            645.2760206304703,
                            eta3,
                            fma(6.674358856924571, eta1, -(171.4826985994883 * eta2)),
                        ),
                ),
                eta1 * fma(
                    692.926425981414,
                    eta4,
                    fma(
                        -143.41033155133266,
                        eta3,
                        fma(
                            25.254017911686333,
                            eta2,
                            fma(-16.57299637423352, eta1, 3.921660532875504),
                        ),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_44_int1` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_44_int1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        eta1 * s1,
        fma(
            -(0.07704777584463054
                * fma(
                    -858.9168637051405,
                    eta3,
                    fma(
                        344.9177692251726,
                        eta2,
                        fma(-50.35070009227704, eta1, 4.581767671445529),
                    ),
                )),
            s2,
            fma(
                -0.062333275821238224,
                fma(
                    -1013.1110639471394,
                    eta3,
                    fma(
                        64.92115530780129,
                        eta2,
                        fma(137.4395894877131, eta1, -21.630297087123807),
                    ),
                ),
                -(0.11014697070998722
                    * fma(
                        -1828.2527520190122,
                        eta3,
                        fma(
                            831.6073263887092,
                            eta2,
                            fma(-108.6912882442823, eta1, 4.149721483857751),
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
                    16.879953490422782,
                    eta3,
                    fma(1.5378890240544967, eta1, -(3.4499418893734903 * eta2)),
                ),
                chidiff2
                    * fma(
                        23.259283336239545,
                        eta3,
                        fma(1.720226708214248, eta1, -(11.87925165364241 * eta2)),
                    ),
            ),
            eta1 * fma(
                1575.602443847111,
                eta4,
                fma(
                    -1111.9864921907126,
                    eta3,
                    fma(
                        324.1998823562892,
                        eta2,
                        fma(-64.95499142822892, eta1, 8.790173464969538),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_44_int2` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_44_int2(wf: &WaveformStruct) -> f64 {
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
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        eta1 * s1,
        fma(
            -0.07801583359561987,
            fma(
                -702172.2934567826,
                eta6,
                fma(
                    760889.6919776166,
                    eta5,
                    fma(
                        -350025.2112501252,
                        eta4,
                        fma(
                            83447.0461509547,
                            eta3,
                            fma(
                                -10671.072606753183,
                                eta2,
                                fma(752.1603553640072, eta1, -28.29972282146242),
                            ),
                        ),
                    ),
                ),
            ),
            0.013159545629626014
                * fma(
                    268756.85659532435,
                    eta6,
                    fma(
                        -1.1505186449682908e6,
                        eta5,
                        fma(
                            1.0141877915334814e6,
                            eta4,
                            fma(
                                -344254.9973814295,
                                eta3,
                                fma(
                                    52391.684517955284,
                                    eta2,
                                    fma(-3557.5003799977294, eta1, 91.1469833190294),
                                ),
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
                    44.78225529547671,
                    eta3,
                    fma(2.3123974306694057, eta1, -(12.237594841284904 * eta2)),
                ),
                chidiff2
                    * fma(
                        61.05270871360041,
                        eta3,
                        fma(2.9282931698944292, eta1, -(25.624210264341933 * eta2)),
                    ),
            ),
            eta1 * fma(
                1478.050456337336,
                eta4,
                fma(
                    -920.358408667518,
                    eta3,
                    fma(
                        236.76146303619544,
                        eta2,
                        fma(-46.81443520117986, eta1, 6.98072197826729),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_44_int3` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_44_int3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let delta = wf.delta;
    let s = wf.chi_pn_hat;
    let chidiff = wf.dchi_half;
    let eta1 = eta;
    let eta2 = eta1 * eta1;
    let eta3 = eta1 * eta2;
    let eta4 = eta1 * eta3;
    let s1 = s;
    let s2 = s1 * s1;
    let chidiff1 = chidiff;
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        eta1 * s1,
        fma(
            -(0.060477799540741804
                * fma(
                    426.0861915646515,
                    eta3,
                    fma(
                        -569.4799973930861,
                        eta2,
                        fma(175.6116682068523, eta1, -14.578189130145661),
                    ),
                )),
            s2,
            fma(
                -0.23316184683282615,
                fma(
                    -1115.4147472977265,
                    eta3,
                    fma(
                        409.3387618483284,
                        eta2,
                        fma(-23.606399143993716, eta1, -1.7279023138971559),
                    ),
                ),
                -(0.09653777612560172
                    * fma(
                        -1815.7529908827103,
                        eta3,
                        fma(
                            541.1248219096527,
                            eta2,
                            fma(-2.1852511802701264, eta1, -5.310643306559746),
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
                    -43.675503209991184,
                    eta3,
                    fma(-0.8765502142143329, eta1, 22.806632458441996 * eta2),
                ),
                chidiff2
                    * fma(
                        16.18571810759235,
                        eta3,
                        fma(0.48698617426180074, eta1, -(4.302527065360426 * eta2)),
                    ),
            ),
            eta1 * fma(
                2379.538739132234,
                eta4,
                fma(
                    -1285.7635006711453,
                    eta3,
                    fma(
                        269.44092930942793,
                        eta2,
                        fma(-44.10631039734796, eta1, 6.379772583015967),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Amp_44_int4` (version 122022): intermediate
/// amplitude collocation point.
pub(crate) fn amp_44_int4(wf: &WaveformStruct) -> f64 {
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
    let chidiff2 = chidiff1 * chidiff1;
    fma(
        eta1 * s1,
        fma(
            -(0.06565232123453196
                * fma(
                    1.0442282862506858e6,
                    eta5,
                    fma(
                        -845341.6001856657,
                        eta4,
                        fma(
                            253694.6032002248,
                            eta3,
                            fma(
                                -33951.35186039629,
                                eta2,
                                fma(1869.0168486099005, eta1, -26.15227471380236),
                            ),
                        ),
                    ),
                )),
            s2,
            fma(
                -0.36068516166901304,
                fma(
                    -91269.3745059139,
                    eta5,
                    fma(
                        60287.45399959349,
                        eta4,
                        fma(
                            -14925.514654896673,
                            eta3,
                            fma(
                                1618.4222330016048,
                                eta2,
                                fma(-47.56175350408845, eta1, -2.120354236840677),
                            ),
                        ),
                    ),
                ),
                -(0.09635801207669747
                    * fma(
                        -67024.98285179552,
                        eta5,
                        fma(
                            -4102.218189945819,
                            eta4,
                            fma(
                                16655.87939259747,
                                eta3,
                                fma(
                                    -4176.398139238679,
                                    eta2,
                                    fma(371.7551657959369, eta1, -11.824692837267394),
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
                    -112.2714974622516,
                    eta3,
                    fma(-2.461738962276138, eta1, 45.3240543970684 * eta2),
                ),
                chidiff2
                    * fma(
                        28.44633544874233,
                        eta3,
                        fma(0.9158352037567031, eta1, -(8.724582331021695 * eta2)),
                    ),
            ),
            eta1 * fma(
                4067.1685640401033,
                eta4,
                fma(
                    -2002.2013283876834,
                    eta3,
                    fma(
                        350.97192927929433,
                        eta2,
                        fma(-45.42463610529546, eta1, 6.098676337298138),
                    ),
                ),
            ),
        ),
    )
}

/// `IMRPhenomXHM_Inter_Phase_21_p1` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_21_p1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -7.86254e6,
        eta6,
        fma(
            6.82533e6,
            eta5,
            fma(
                -2.28832e6,
                eta4,
                fma(
                    369153.0,
                    eta3,
                    fma(-23428.1, eta2, fma(-1956.93, eta, 4045.84 + 7.63226 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        eta2,
        fma(
            2689.18,
            s4,
            fma(2392.73, s3, fma(-619.837, s, -(722.787 * s2))),
        ),
        fma(
            92.6654,
            s4,
            fma(
                eta,
                fma(
                    -1077.75,
                    s4,
                    fma(1114.51, s3, fma(-79.761, s, -(355.299 * s2))),
                ),
                fma(-355.67, s3, fma(-347.273, s, 83.5428 * s2))
                    + fma(13.6971, s3, fma(4.44457, s, 16.5548 * s2)) / eta,
            ),
        ),
    );
    let uneq_spin = fma(
        fma(
            918.976 * wf.chi1l,
            fma(-4.0, eta, 1.0).sqrt(),
            -(918.976 * wf.chi2l * fma(-4.0, eta, 1.0).sqrt()),
        ),
        eta,
        fma(
            91.7679 * wf.chi1l,
            fma(-4.0, eta, 1.0).sqrt(),
            -(91.7679 * wf.chi2l * fma(-4.0, eta, 1.0).sqrt()),
        ) * eta2,
    );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_21_p2` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_21_p2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta3 * eta2;
    let eta6 = eta4 * eta2;
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s3 * s;
    let no_spin = fma(
        -7.86254e6,
        eta6,
        fma(
            6.82533e6,
            eta5,
            fma(
                -2.28832e6,
                eta4,
                fma(
                    369153.0,
                    eta3,
                    fma(-27556.2, eta2, fma(194.72, eta, 3509.09 + 0.91868 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        eta * fma(2689.18, eta2, fma(-1077.75, eta, 92.6654)),
        s4,
        fma(
            fma(2392.73, eta3, fma(-1180.4, eta2, fma(117.913, eta, 5.7791))),
            s3,
            fma(
                fma(
                    -619.837,
                    eta3,
                    fma(131.815, eta2, fma(-60.1611, eta, 0.7083999999999999)),
                ),
                s,
                fma(
                    -722.787,
                    eta3,
                    fma(278.588, eta2, fma(-59.2068, eta, 6.104720000000001)),
                ) * s2,
            ),
        ),
    ) / eta;
    let uneq_spin = -91.7679
        * fma(-4.0, eta, 1.0).sqrt()
        * eta
        * fma(
            wf.chi1l,
            fma(-1.0, eta, -1.6012352903357276),
            wf.chi2l * fma(1.0, eta, 1.6012352903357276),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_21_p3` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_21_p3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let delta = fma(-4.0, eta, 1.0).sqrt();
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta2 * eta2;
    let eta5 = eta2 * eta3;
    let eta6 = eta3 * eta3;
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -7.86254e6,
        eta6,
        fma(
            6.82533e6,
            eta5,
            fma(
                -2.28832e6,
                eta4,
                fma(
                    369153.0,
                    eta3,
                    fma(-28651.9, eta2, fma(890.016, eta, 3241.68)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(2689.18, eta2, fma(-1077.75, eta, 92.6654)),
        s4,
        fma(
            fma(2392.73, eta2, fma(-1094.59, eta, 117.913)),
            s3,
            fma(
                fma(-619.837, eta2, fma(187.641, eta, -2.2484)),
                s,
                fma(-722.787, eta2, fma(166.323, eta, 3.22603)) * s2,
            ),
        ),
    );
    let uneq_spin = 91.7679 * wf.dchi * delta * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_21_p4` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_21_p4(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -7.86254e6,
        eta6,
        fma(
            6.82533e6,
            eta5,
            fma(
                -2.28832e6,
                eta4,
                fma(
                    369780.0,
                    eta3,
                    fma(-28932.5, eta2, fma(974.355, eta, 3160.88)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(-1096.81, eta2, fma(577.057, eta, -61.7446)),
        s4,
        fma(
            fma(-1613.0, eta2, fma(467.057, eta, -20.0261)),
            s3,
            fma(
                fma(438.401, eta2, fma(-196.851, eta, 26.3355)),
                s,
                fma(117.563, eta2, fma(-256.248, eta, 45.9957)) * s2,
            ),
        ),
    );
    let uneq_spin = 65.3326 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_21_p5` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_21_p5(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        3635.76,
        eta3,
        fma(-1688.26, eta2, fma(315.911, eta, 3102.36)),
    );
    let eq_spin = fma(
        fma(4689.22, eta2, fma(-1243.26, eta, 81.8379)),
        s4,
        fma(
            fma(1094.11, eta2, fma(-365.842, eta, 40.7054)),
            s3,
            fma(
                fma(-1029.76, eta2, fma(320.93, eta, -23.0959)),
                s,
                fma(-3079.39, eta2, fma(826.816, eta, -49.5435)) * s2,
            ),
        ),
    );
    let uneq_spin = 119.014 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_21_p6` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_21_p6(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -255.245,
        eta3,
        fma(190.008, eta2, fma(4.89194, eta, 3089.18)),
    );
    let eq_spin = fma(
        fma(3967.43, eta2, fma(-722.674, eta, -12.8798)),
        s4,
        fma(
            fma(549.085, eta2, -24.6193),
            s3,
            fma(
                fma(-432.223, eta2, fma(57.1612, eta, 2.96997)),
                s,
                fma(-2804.66, eta2, fma(630.516, eta, -18.8929)) * s2,
            ),
        ),
    );
    let uneq_spin = 74.0984 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_33_p1` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_33_p1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let no_spin = fma(
        -1.15659e7,
        eta6,
        fma(
            1.01799e7,
            eta5,
            fma(
                -3.22792e6,
                eta4,
                fma(
                    371303.00000000006,
                    eta3,
                    fma(18485.9, eta2, fma(-8727.4, eta, 4360.19 + 4.27128 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            -34860.2,
            eta4,
            fma(
                19614.6,
                eta3,
                fma(-3255.6400000000003, eta2, fma(-251.579, eta, 11.6635)),
            ),
        ),
        s,
        fma(
            -74299.5,
            eta4,
            fma(
                36587.3,
                eta3,
                fma(-5421.92, eta2, fma(204.025, eta, 14.8017)),
            ),
        ) * s2,
    ) / eta;
    let uneq_spin = eta
        * fma(
            223.65100000000004 * wf.chi1l * fma(-4.0, eta, 1.0).sqrt(),
            fma(1.0, eta, 3.9201300240106223),
            -(223.65100000000004
                * wf.chi2l
                * fma(-4.0, eta, 1.0).sqrt()
                * fma(1.0, eta, 3.9201300240106223)),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_33_p2` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_33_p2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let no_spin = fma(
        -1.15659e7,
        eta6,
        fma(
            1.01799e7,
            eta5,
            fma(
                -3.41708e6,
                eta4,
                fma(
                    518314.99999999994,
                    eta3,
                    fma(-25514.0, eta2, fma(-2397.09, eta, 3797.06 + 0.786684 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            -34860.2,
            eta4,
            fma(
                19614.6,
                eta3,
                fma(-3520.37, eta2, fma(39.4668, eta, 6.7812399999999995)),
            ),
        ),
        s,
        fma(
            -74299.5,
            eta4,
            fma(
                36587.3,
                eta3,
                fma(-5914.61, eta2, fma(293.215, eta, 4.80384)),
            ),
        ) * s2,
    ) / eta;
    let uneq_spin = -223.65100000000004
        * fma(-4.0, eta, 1.0).sqrt()
        * eta
        * fma(
            wf.chi1l,
            fma(-1.0, eta, -1.3095134830606614),
            wf.chi2l * fma(1.0, eta, 1.3095134830606614),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_33_p3` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_33_p3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let no_spin = fma(
        -1.15659e7,
        eta6,
        fma(
            1.01799e7,
            eta5,
            fma(
                -3.52532e6,
                eta4,
                fma(
                    605028.0,
                    eta3,
                    fma(-52406.1, eta2, fma(1796.03, eta, 3321.83)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            -34860.2,
            eta3,
            fma(19614.6, eta2, fma(-3714.77, eta, 223.601)),
        ),
        s,
        fma(
            -74299.5,
            eta3,
            fma(36587.3, eta2, fma(-5906.46, eta, 314.317)),
        ) * s2,
    );
    let uneq_spin = 223.651 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_33_p4` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_33_p4(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -1.70363e6,
        eta6,
        fma(
            1.17255e6,
            eta5,
            fma(
                -248477.0,
                eta4,
                fma(3456.2, eta3, fma(5139.79, eta2, fma(-661.15, eta, 3239.44))),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            67102.4,
            eta3,
            fma(-35337.1, eta2, fma(5946.15, eta, -310.489)),
        ),
        s3,
        fma(
            fma(
                -46537.3,
                eta3,
                fma(24364.0, eta2, fma(-4150.09, eta, 225.859)),
            ),
            s,
            fma(
                -23603.5,
                eta3,
                fma(8953.98, eta2, fma(-994.971, eta, 35.2439)),
            ) * s2,
        ),
    );
    let uneq_spin = 30.484 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_33_p5` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_33_p5(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -1.17123e7,
        eta6,
        fma(
            9.99745e6,
            eta5,
            fma(
                -3.35991e6,
                eta4,
                fma(
                    563997.0,
                    eta3,
                    fma(-49428.3, eta2, fma(2143.06, eta, 3114.3)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            67102.4,
            eta3,
            fma(-33001.4, eta2, fma(5179.97, eta, -257.524)),
        ),
        s3,
        fma(
            fma(
                -46537.3,
                eta3,
                fma(23046.2, eta2, fma(-3705.08, eta, 190.051)),
            ),
            s,
            fma(
                -23603.5,
                eta3,
                fma(10166.1, eta2, fma(-1414.2, eta, 63.6615)),
            ) * s2,
        ),
    );
    let uneq_spin = 54.9833 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_33_p6` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_33_p6(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -4.64858e6,
        eta6,
        fma(
            3.79886e6,
            eta5,
            fma(
                -1.19313e6,
                eta4,
                fma(
                    179537.0,
                    eta3,
                    fma(-13003.6, eta2, fma(384.121, eta, 3111.46)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            75501.8,
            eta3,
            fma(-37382.4, eta2, fma(5886.62, eta, -292.855)),
        ),
        s3,
        fma(
            fma(
                -50165.9,
                eta3,
                fma(24532.9, eta2, fma(-3834.22, eta, 182.864)),
            ),
            s,
            fma(
                -17842.3,
                eta3,
                fma(6701.33, eta2, fma(-746.957, eta, 21.0158)),
            ) * s2,
        ),
    );
    let uneq_spin = 75.5162 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_32_p1` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_32_p1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -3.26362e6,
        eta6,
        fma(
            2.86393e6,
            eta5,
            fma(
                -704442.0,
                eta4,
                fma(
                    -64068.40000000001,
                    eta3,
                    fma(58234.6, eta2, fma(-10687.8, eta, 4414.11 + 4.21564 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        eta * fma(4810.61, eta2, fma(-1569.0299999999997, eta, 106.047)),
        s4,
        fma(
            fma(
                294.289,
                eta3,
                fma(-541.261, eta2, fma(43.8812, eta, 10.4333)),
            ),
            s3,
            fma(
                fma(
                    -3970.89,
                    eta3,
                    fma(2095.72, eta2, fma(-610.267, eta, 6.39833)),
                ),
                s,
                fma(
                    -794.79,
                    eta3,
                    fma(331.593, eta2, fma(-99.1551, eta, 22.956700000000005)),
                ) * s2,
            ),
        ),
    ) / eta;
    let uneq_spin = 132.244
        * fma(-4.0, eta, 1.0).sqrt()
        * eta
        * fma(
            wf.chi1l,
            fma(-1.0, eta, 6.227738120444028),
            wf.chi2l * fma(1.0, eta, -6.227738120444028),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_32_p2` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_32_p2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -3.26362e6,
        eta6,
        fma(
            2.86393e6,
            eta5,
            fma(
                -807629.0,
                eta4,
                fma(
                    24484.2,
                    eta3,
                    fma(29218.1, eta2, fma(-6202.38, eta, 3980.7 + 0.956703 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        eta * fma(4810.61, eta2, fma(-1569.0299999999997, eta, 106.047)),
        s4,
        fma(
            fma(
                294.289,
                eta3,
                fma(-541.261, eta2, fma(43.8812, eta, 10.4333)),
            ),
            s3,
            fma(
                fma(
                    1291.56,
                    eta3,
                    fma(75.246, eta2, fma(-226.825, eta, 1.92692)),
                ),
                s,
                fma(
                    -2402.94,
                    eta3,
                    fma(608.328, eta2, fma(-99.1551, eta, 15.328700000000001)),
                ) * s2,
            ),
        ),
    ) / eta;
    let uneq_spin = 132.244
        * fma(-4.0, eta, 1.0).sqrt()
        * eta
        * fma(
            wf.chi1l,
            fma(-1.0, eta, 2.5769789177580837),
            wf.chi2l * fma(1.0, eta, -2.5769789177580837),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_32_p3` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_32_p3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -2.0065e7,
        eta6,
        fma(
            1.76399e7,
            eta5,
            fma(
                -6.0644e6,
                eta4,
                fma(
                    1.01936e6,
                    eta3,
                    fma(-84042.9, eta2, fma(2308.63, eta, 3416.57)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(4810.61, eta2, fma(-1569.03, eta, 106.047)),
        s4,
        fma(
            fma(
                -114029.0,
                eta3,
                fma(61581.0, eta2, fma(-10627.4, eta, 559.705)),
            ),
            s3,
            fma(
                fma(
                    12750.0,
                    eta3,
                    fma(-2582.55, eta2, fma(-282.354, eta, 24.6295)),
                ),
                s,
                fma(
                    -114798.0,
                    eta3,
                    fma(56407.8, eta2, fma(-8775.86, eta, 433.675)),
                ) * s2,
            ),
        ),
    );
    let uneq_spin = 63.9466 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_32_p4` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_32_p4(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = fma(
        -3.26362e6,
        eta6,
        fma(
            2.86393e6,
            eta5,
            fma(
                -919108.0,
                eta4,
                fma(
                    127610.0,
                    eta3,
                    fma(-5980.37, eta2, fma(-476.909, eta, 3307.49)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(4810.61, eta2, fma(-1569.03, eta, 106.047)),
        s4,
        fma(
            fma(294.289, eta2, fma(-370.362, eta, 43.8812)),
            s3,
            fma(
                fma(1291.56, eta2, fma(-282.354, eta, -5.02553)),
                s,
                fma(-2402.94, eta2, fma(740.123, eta, -43.8823)) * s2,
            ),
        ),
    );
    let uneq_spin = -132.244 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_32_p5` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_32_p5(wf: &WaveformStruct) -> f64 {
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
        -3.15866e8,
        eta7,
        fma(
            2.75056e8,
            eta6,
            fma(
                -9.56988e7,
                eta5,
                fma(
                    1.73811e7,
                    eta4,
                    fma(
                        -1.81883e6,
                        eta3,
                        fma(111203.0, eta2, fma(-3967.58, eta, 3259.03)),
                    ),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(4810.61, eta2, fma(-1569.03, eta, 106.047)),
        s4,
        fma(
            fma(-6233.99, eta2, fma(2183.43, eta, -201.633)),
            s3,
            fma(
                fma(3810.18, eta2, fma(-1104.53, eta, 19.7509)),
                s,
                fma(-5944.49, eta2, fma(2314.51, eta, -230.07)) * s2,
            ),
        ),
    );
    let uneq_spin = 112.714 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_32_p6` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_32_p6(wf: &WaveformStruct) -> f64 {
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
        -3.15866e8,
        eta7,
        fma(
            2.75056e8,
            eta6,
            fma(
                -9.56988e7,
                eta5,
                fma(
                    1.73811e7,
                    eta4,
                    fma(
                        -1.81883e6,
                        eta3,
                        fma(111203.0, eta2, fma(-3967.58, eta, 3259.03)),
                    ),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(4810.61, eta2, fma(-1569.03, eta, 106.047)),
        s4,
        fma(
            fma(-6233.99, eta2, fma(2183.43, eta, -201.633)),
            s3,
            fma(
                fma(3810.18, eta2, fma(-1104.53, eta, 19.7509)),
                s,
                fma(-5944.49, eta2, fma(2314.51, eta, -230.07)) * s2,
            ),
        ),
    );
    let uneq_spin = 112.714 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_44_p1` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_44_p1(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -1.60802e7,
        eta6,
        fma(
            1.37792e7,
            eta5,
            fma(
                -4.33197e6,
                eta4,
                fma(
                    536500.0,
                    eta3,
                    fma(5534.1, eta2, fma(-8202.33, eta, 4349.66 + 4.34125 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        eta * fma(
            -109418.0,
            eta3,
            fma(57303.4, eta2, fma(-9508.5, eta, 492.60300000000007)),
        ),
        s3,
        fma(
            fma(
                17900.9,
                eta4,
                fma(
                    -9349.73,
                    eta3,
                    fma(1822.9100000000003, eta2, fma(-528.098, eta, 12.0704)),
                ),
            ),
            s,
            fma(
                -71523.0,
                eta4,
                fma(
                    35416.6,
                    eta3,
                    fma(-5452.04, eta2, fma(253.334, eta, 10.4092)),
                ),
            ) * s2,
        ),
    ) / eta;
    let uneq_spin = -262.143
        * fma(-4.0, eta, 1.0).sqrt()
        * eta
        * fma(
            wf.chi1l,
            fma(-1.0, eta, -3.0782778864970646),
            wf.chi2l * fma(1.0, eta, 3.0782778864970646),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_44_p2` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_44_p2(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -1.60802e7,
        eta6,
        fma(
            1.37792e7,
            eta5,
            fma(
                -4.50145e6,
                eta4,
                fma(
                    665951.0,
                    eta3,
                    fma(-33475.8, eta2, fma(-2421.77, eta, 3804.19 + 0.66144 / eta)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        eta * fma(
            -109418.0,
            eta3,
            fma(57303.4, eta2, fma(-9508.5, eta, 492.60300000000007)),
        ),
        s3,
        fma(
            fma(
                17900.9,
                eta4,
                fma(
                    -7676.87,
                    eta3,
                    fma(926.576, eta2, fma(-172.047, eta, 5.83038)),
                ),
            ),
            s,
            fma(
                -71523.0,
                eta4,
                fma(
                    35722.1,
                    eta3,
                    fma(-5672.02, eta2, fma(253.334, eta, 6.17601)),
                ),
            ) * s2,
        ),
    ) / eta;
    let uneq_spin = -262.143
        * fma(-4.0, eta, 1.0).sqrt()
        * eta
        * fma(
            wf.chi1l,
            fma(-1.0, eta, -1.0543062374352932),
            wf.chi2l * fma(1.0, eta, 1.0543062374352932),
        );
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_44_p3` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_44_p3(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -1.60802e7,
        eta6,
        fma(
            1.37792e7,
            eta5,
            fma(
                -4.64438e6,
                eta4,
                fma(
                    777272.0,
                    eta3,
                    fma(-66340.1, eta2, fma(2353.58, eta, 3308.97)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            -109418.0,
            eta3,
            fma(57303.4, eta2, fma(-9508.5, eta, 492.603)),
        ),
        s3,
        fma(
            fma(
                17900.9,
                eta3,
                fma(-7989.26, eta2, fma(926.576, eta, -21.5697)),
            ),
            s,
            fma(
                -71523.0,
                eta3,
                fma(37599.5, eta2, fma(-6403.24, eta, 353.539)),
            ) * s2,
        ),
    );
    let uneq_spin = 262.143 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_44_p4` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_44_p4(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -1.31138e6,
        eta6,
        fma(
            908279.0,
            eta5,
            fma(
                -165169.0,
                eta4,
                fma(
                    -17422.6,
                    eta3,
                    fma(8463.89, eta2, fma(-928.56, eta, 3245.63)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            9582.58,
            eta3,
            fma(-4679.24, eta2, fma(620.695, eta, -15.7439)),
        ),
        s3,
        fma(
            fma(
                -6758.52,
                eta3,
                fma(3536.61, eta2, fma(-590.293, eta, 32.506)),
            ),
            s,
            fma(
                9129.45,
                eta3,
                fma(-4867.87, eta2, fma(738.141, eta, -25.7716)),
            ) * s2,
        ),
    );
    let uneq_spin = 87.0832 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_44_p5` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_44_p5(wf: &WaveformStruct) -> f64 {
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
    let no_spin = fma(
        2.2829700000000003e8,
        eta7,
        fma(
            -2.17143e8,
            eta6,
            fma(
                8.39194e7,
                eta5,
                fma(
                    -1.69796e7,
                    eta4,
                    fma(
                        1.92148e6,
                        eta3,
                        fma(-119588.0, eta2, fma(3722.46, eta, 3108.38)),
                    ),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(999.417, eta2, fma(-356.241, eta, 20.3384)),
        s3,
        fma(
            fma(-529.854, eta, 118.319) * eta,
            s,
            fma(516.333, eta2, fma(-240.648, eta, 21.0314)) * s2,
        ),
    );
    let uneq_spin = 97.1364 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `IMRPhenomXHM_Inter_Phase_44_p6` (version 122019): intermediate
/// phase-derivative collocation point.
pub(crate) fn phase_44_p6(wf: &WaveformStruct) -> f64 {
    let eta = wf.eta;
    let s = wf.s_tot_r;
    let eta2 = eta * eta;
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let eta5 = eta.powf(5.0);
    let eta6 = eta.powf(6.0);
    let s2 = s * s;
    let s3 = s.powf(3.0);
    let no_spin = fma(
        -6.01179e6,
        eta6,
        fma(
            4.29193e6,
            eta5,
            fma(
                -1.31523e6,
                eta4,
                fma(
                    220332.0,
                    eta3,
                    fma(-20371.1, eta2, fma(986.752, eta, 3096.03)),
                ),
            ),
        ),
    );
    let eq_spin = fma(
        fma(
            9582.58,
            eta3,
            fma(-3070.08, eta2, fma(70.1931, eta, 17.0209)),
        ),
        s3,
        fma(
            fma(
                -6758.52,
                eta3,
                fma(2255.76, eta2, fma(-118.526, eta, -9.96292)),
            ),
            s,
            fma(
                9129.45,
                eta3,
                fma(-3605.8, eta2, fma(370.039, eta, -14.4869)),
            ) * s2,
        ),
    );
    let uneq_spin = 23.0759 * wf.dchi * fma(-4.0, eta, 1.0).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

/* End of the parameter-space fits. */

/// Intermediate amplitude collocation fits `int1..int4` of mode `mode_int`.
fn amp_fit(mode_int: usize, k: usize, wf: &WaveformStruct) -> f64 {
    const FITS: [[Fit; 4]; 4] = [
        [amp_21_int1, amp_21_int2, amp_21_int3, amp_21_int4],
        [amp_33_int1, amp_33_int2, amp_33_int3, amp_33_int4],
        [amp_32_int1, amp_32_int2, amp_32_int3, amp_32_int4],
        [amp_44_int1, amp_44_int2, amp_44_int3, amp_44_int4],
    ];
    FITS[mode_int][k](wf)
}

/// `pPhase->IntermediatePhaseFits`: the phase-derivative collocation fits
/// `p1..p6` of mode `mode_int`.
pub(crate) fn phase_fit(mode_int: usize, i: usize, wf: &WaveformStruct) -> f64 {
    const FITS: [[Fit; 6]; 4] = [
        [
            phase_21_p1,
            phase_21_p2,
            phase_21_p3,
            phase_21_p4,
            phase_21_p5,
            phase_21_p6,
        ],
        [
            phase_33_p1,
            phase_33_p2,
            phase_33_p3,
            phase_33_p4,
            phase_33_p5,
            phase_33_p6,
        ],
        [
            phase_32_p1,
            phase_32_p2,
            phase_32_p3,
            phase_32_p4,
            phase_32_p5,
            phase_32_p6,
        ],
        [
            phase_44_p1,
            phase_44_p2,
            phase_44_p3,
            phase_44_p4,
            phase_44_p5,
            phase_44_p6,
        ],
    ];
    FITS[mode_int][i](wf)
}

/// `IMRPhenomXHM_Intermediate_Amp_Ansatz` (122022 release): a polynomial
/// times f^(-7/6), in strain units without `amp0`.
#[inline]
pub(crate) fn amp_ansatz(powers: &UsefulPowers, amp: &HMAmpCoefficients) -> f64 {
    let mut result = 0.;
    let mut fpower = 1.;
    for &coefficient in &amp.inter_coefficients[..amp.n_coefficients_inter] {
        result = fma(coefficient, fpower, result);
        fpower *= powers.itself;
    }
    result * powers.m_seven_sixths
}

/// `IMRPhenomXHM_Intermediate_Amp_CollocationPoints` and
/// `IMRPhenomXHM_Intermediate_Amp_Coefficients` (122022 release):
/// equispaced collocation points between the inspiral and ringdown
/// boundaries, with values (and derivatives) from the inspiral and ringdown
/// ansatzes at the boundaries.
///
/// `rotated_32` evaluates the spheroidal-to-spherical rotation of the (3,2)
/// ringdown for the coefficients computed so far, with its first and second
/// derivatives.
pub(crate) fn intermediate_amp_coefficients<F>(
    amp: &mut HMAmpCoefficients,
    wfhm: &HMWaveformStruct,
    wf: &WaveformStruct,
    rotated_32: F,
) where
    F: Fn(f64, &HMAmpCoefficients) -> [Complex<f64>; 3],
{
    let version = wfhm.intermediate_amp_version;
    let n_colloc = version.len();
    let deltaf = (amp.f_amp_match_im - amp.f_amp_match_in) / (n_colloc - 1) as f64;
    let mut idx = 0;
    for (i, &digit) in version.iter().enumerate() {
        let f = fma(deltaf, i as f64, amp.f_amp_match_in);
        if digit >= 1 {
            amp.colloc_freqs_inter[idx] = f;
        }
        if digit == 2 {
            amp.colloc_freqs_inter[idx + 1] = f;
        }
        idx += digit as usize;
    }

    // Values at the inspiral boundary.
    let powers_of_finsp = UsefulPowers::new(amp.f_amp_match_in);
    let mut n = 0;
    match version[0] {
        1 => {
            amp.colloc_values_inter[0] = hm_inspiral::amp_ansatz(&powers_of_finsp, wfhm, amp);
            n += 1;
        }
        2 => {
            amp.colloc_values_inter[0] = hm_inspiral::amp_ansatz(&powers_of_finsp, wfhm, amp);
            amp.colloc_values_inter[1] = hm_inspiral::amp_nd_ansatz(&powers_of_finsp, wfhm, amp);
            n += 2;
        }
        _ => {}
    }
    // Parameter-space fits of the inner points.
    for (i, &digit) in version.iter().enumerate().take(n_colloc - 1).skip(1) {
        if digit == 1 {
            amp.colloc_values_inter[n] = amp_fit(wfhm.mode_int, i - 1, wf).abs();
            n += 1;
        }
    }
    // Value and derivative at the ringdown boundary.
    let powers_of_frd = UsefulPowers::new(amp.f_amp_match_im);
    let (rd_value, rd_derivative) = if wfhm.mixing_on {
        // LALSimulation (IMRPhenomXHM_RD_Amp_NDAnsatz) takes the derivative
        // by fourth-order finite differences with step 1e-9.
        let h = rotated_32(amp.f_amp_match_im, amp);
        (h[0].re.hypot(h[0].im), jet::abs_derivative(h))
    } else {
        (
            hm_ringdown::amp_ansatz(&powers_of_frd, wfhm, amp),
            hm_ringdown::amp_dansatz(&powers_of_frd, wfhm, amp),
        )
    };
    amp.colloc_values_inter[n] = rd_value;
    n += 1;
    if version[n_colloc - 1] == 2 {
        amp.colloc_values_inter[n] = rd_derivative;
    }

    // Linear system for the coefficients of the polynomial times f^(-7/6).
    let n_coefficients = amp.n_coefficients_inter;
    let mut matrix = [[0.; 8]; 8];
    let mut b = [0.; 8];
    let seven_sixths = 7. / 6.;
    let mut row = 0;
    for &digit in &version {
        if digit == 0 {
            continue;
        }
        let fcoll = amp.colloc_freqs_inter[row];
        let fcoll_m_seven_sixths = fcoll.powf(-seven_sixths);
        b[row] = amp.colloc_values_inter[row];
        let mut fpower = 1.;
        for j in 0..n_coefficients {
            matrix[row][j] = fpower * fcoll_m_seven_sixths;
            fpower *= fcoll;
        }
        row += 1;
        if digit == 2 {
            b[row] = amp.colloc_values_inter[row];
            let mut fpower = 1. / fcoll;
            for j in 0..n_coefficients {
                matrix[row][j] = (j as f64 - seven_sixths) * fpower * fcoll_m_seven_sixths;
                fpower *= fcoll;
            }
            row += 1;
        }
    }
    let x = match n_coefficients {
        8 => lu_solve(matrix, b),
        5 => {
            let mut small = [[0.; 5]; 5];
            for (s, m) in small.iter_mut().zip(&matrix) {
                s.copy_from_slice(&m[..5]);
            }
            let x5 = lu_solve(small, [b[0], b[1], b[2], b[3], b[4]]);
            [x5[0], x5[1], x5[2], x5[3], x5[4], 0., 0., 0.]
        }
        _ => unreachable!("IMRPhenomXHM: unexpected number of intermediate coefficients"),
    };
    amp.inter_coefficients = x;
}

/// `IMRPhenomXHM_Inter_Phase_Ansatz`: the intermediate phase derivative.
pub(crate) fn phase_ansatz(
    ff: f64,
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    phase: &HMPhaseCoefficients,
) -> f64 {
    let fda = wfhm.f_damp;
    let frd = wfhm.f_ring;
    let lorentzian = (phase.c_l * fda) / fma(fda, fda, (ff - frd) * (ff - frd));
    let mut dphase = fma(phase.c1, powers.m_one, phase.c0);
    dphase = fma(phase.c2, powers.m_two, dphase);
    dphase = fma(phase.c4, powers.m_four, dphase);
    if wfhm.mode_tag == 32 {
        dphase = fma(phase.c3, powers.m_three, dphase);
    }
    dphase + lorentzian
}

/// `IMRPhenomXHM_Inter_Phase_AnsatzInt`: the intermediate phase.
#[inline]
pub(crate) fn phase_ansatz_int(
    ff: f64,
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    phase: &HMPhaseCoefficients,
) -> f64 {
    let fda = wfhm.f_damp;
    let frd = wfhm.f_ring;
    let mut phase_ir = fma(phase.c0, ff, phase.c1 * powers.log);
    phase_ir = fma(-phase.c2, powers.m_one, phase_ir);
    phase_ir = fma(-(1. / 3. * phase.c4), powers.m_three, phase_ir);
    if wfhm.mode_tag == 32 {
        phase_ir = fma(-(0.5 * phase.c3), powers.m_two, phase_ir);
    }
    fma(phase.c_l, ((ff - frd) / fda).atan(), phase_ir)
}
