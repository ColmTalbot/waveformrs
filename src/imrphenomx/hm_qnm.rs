//! Quasi-normal mode fits of the higher modes, `LALSimIMRPhenomXHM_qnm.c`.
//!
//! The ringdown and damping frequencies of the (2,1), (3,3), (3,2) and (4,4)
//! modes are in units of the remnant mass (divide by the final mass in units
//! of the total mass), and the spheroidal-spherical mixing coefficients of
//! the (3,2) mode follow Berti & Klein (2014). The fits are evaluated with
//! LAL's fused multiply-adds: their denominators cancel strongly near
//! extremal spins.

use crate::imrphenomx::internals::fma;

/// `evaluate_QNMfit_fring21`: ringdown frequency of the (2,1) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fring21(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    fma(
        -0.0003387374454044957,
        x5,
        fma(
            0.001164879406179587,
            x4,
            fma(
                -0.0018964744613388146,
                x3,
                fma(
                    0.021967909664591865,
                    x2,
                    fma(-0.07585416297991414, af, 0.059471695665734674),
                ),
            ),
        ),
    ) / fma(0.49246920313191234, x2, fma(-1.4437415542456158, af, 1.0))
}

/// `evaluate_QNMfit_fring33`: ringdown frequency of the (3,3) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fring33(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    fma(
        0.0027301732074159835,
        x6,
        fma(
            -0.006756504382964637,
            x5,
            fma(
                -0.030848060170259615,
                x4,
                fma(
                    0.03343753057911253,
                    x3,
                    fma(
                        0.13402916709362475,
                        x2,
                        fma(-0.22799517865876945, af, 0.09540436245212061),
                    ),
                ),
            ),
        ),
    ) / fma(
        0.05321818246993958,
        x6,
        fma(
            -0.4706873667569393,
            x4,
            fma(2.144070539525238, x2, fma(-2.7265947806178334, af, 1.0)),
        ),
    )
}

/// `evaluate_QNMfit_fring32`: ringdown frequency of the (3,2) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fring32(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    fma(
        -0.0015206141948469621,
        x6,
        fma(
            0.0020533816327907334,
            x5,
            fma(
                0.0118277880067919,
                x4,
                fma(
                    -0.000673589757007597,
                    x3,
                    fma(
                        0.030099881830507727,
                        x2,
                        fma(-0.13628306966373951, af, 0.09540436245212061),
                    ),
                ),
            ),
        ),
    ) / fma(
        -0.027297817699401976,
        x6,
        fma(
            0.12256204148002939,
            x4,
            fma(0.5634705514193629, x2, fma(-1.6531854335715193, af, 1.0)),
        ),
    )
}

/// `evaluate_QNMfit_fring44`: ringdown frequency of the (4,4) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fring44(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    fma(
        0.00034775290179000503,
        x6,
        fma(
            0.000045483523029172406,
            x5,
            fma(
                -0.002795401084713644,
                x4,
                fma(
                    0.015487322972031054,
                    x3,
                    fma(
                        0.0710926778043916,
                        x2,
                        fma(-0.21224284094693793, af, 0.1287821193485683),
                    ),
                ),
            ),
        ),
    ) / fma(
        -0.06378640753152783,
        x4,
        fma(1.0593147376898773, x2, fma(-1.9931645124693607, af, 1.0)),
    )
}

/// `evaluate_QNMfit_fdamp21`: damping frequency of the (2,1) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fdamp21(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    fma(
        0.002943057145913646,
        x5,
        fma(
            -0.2947372598589144,
            x4,
            fma(
                1.2194717985037946,
                x3,
                fma(
                    0.14456081596393977,
                    x2,
                    fma(-3.1358071947583093, af, 2.0696914454467294),
                ),
            ),
        ),
    ) / fma(
        -18.975287709794745,
        x4,
        fma(
            75.90115083917898,
            x3,
            fma(
                17.7141194900164,
                x2,
                fma(-219.81790388304876, af, 146.1779212636481),
            ),
        ),
    )
}

/// `evaluate_QNMfit_fdamp33`: damping frequency of the (3,3) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fdamp33(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    fma(
        -0.0001414350555699256,
        x5,
        fma(
            -0.0015925124814622795,
            x4,
            fma(
                0.001034954865629645,
                x3,
                fma(
                    0.017192623913708124,
                    x2,
                    fma(-0.03124423610028678, af, 0.014754148319335946),
                ),
            ),
        ),
    ) / fma(
        -0.09874113387889819,
        x4,
        fma(1.196809702382645, x2, fma(-2.0963684630756894, af, 1.0)),
    )
}

/// `evaluate_QNMfit_fdamp32`: damping frequency of the (3,2) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fdamp32(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    fma(
        -0.0034761714223258693,
        x4,
        fma(
            0.0014945908223317514,
            x3,
            fma(
                0.02168855041940869,
                x2,
                fma(-0.03445752346074498, af, 0.014754148319335946),
            ),
        ),
    ) / fma(
        -0.18791187563554512,
        x4,
        fma(1.5096146036915865, x2, fma(-2.320722660848874, af, 1.0)),
    )
}

/// `evaluate_QNMfit_fdamp44`: damping frequency of the (4,4) mode
/// for dimensionless final spin `af`, in units of the remnant mass.
pub(crate) fn fdamp44(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    fma(
        0.000812185411753066,
        x6,
        fma(
            0.00031983835498725354,
            x5,
            fma(
                0.002510528746148588,
                x4,
                fma(
                    0.0002837322846047305,
                    x3,
                    fma(
                        -0.0016734788189065538,
                        x2,
                        fma(-0.01722587715950451, af, 0.014986847152355699),
                    ),
                ),
            ),
        ),
    ) / fma(
        0.051876225199833995,
        x6,
        fma(
            0.13983808071522857,
            x4,
            fma(-0.0500827971270845, x2, fma(-1.1350205970682399, af, 1.0)),
        ),
    )
}

/// `evaluate_QNMfit_re_l2m2lp2`: real part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn re_l2m2lp2(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    fma(
        -0.0006267849034466508,
        x5,
        fma(
            -0.1695667458204109,
            x4,
            fma(
                0.0043296365593147035,
                x3,
                fma(1.461988775298876, x2, fma(-2.2956993576253635, af, 1.0)),
            ),
        ),
    ) / fma(
        -0.00007150540890128118,
        x6,
        fma(
            -0.16843226886562457,
            x4,
            fma(1.4646339137818438, x2, fma(-2.2956977727459043, af, 1.0)),
        ),
    )
}

/// `evaluate_QNMfit_re_l2m2lp3`: real part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn re_l2m2lp3(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    af * fma(
        -0.37117112862247975,
        x5,
        fma(
            -2.210028290847915,
            x4,
            fma(
                3.0355668591803386,
                x3,
                fma(
                    19.030390708998894,
                    x2,
                    fma(-37.978140351289014, af, 18.522563276099167),
                ),
            ),
        ),
    ) / fma(
        1.0,
        x6,
        fma(
            -30.79738566181734,
            x4,
            fma(
                243.3353695550844,
                x2,
                fma(-377.9093045285145, af, 164.52480238697507),
            ),
        ),
    )
}

/// `evaluate_QNMfit_re_l3m2lp2`: real part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn re_l3m2lp2(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    af * fma(
        -0.006557249133752502,
        x5,
        fma(
            -0.03620067426672167,
            x4,
            fma(
                0.0855565148647794,
                x3,
                fma(
                    0.3844811236426182,
                    x2,
                    fma(-0.9016636384605536, af, 0.47513455283841244),
                ),
            ),
        ),
    ) / fma(
        1.0,
        x4,
        fma(
            -9.406169787571082,
            x2,
            fma(15.170831931186493, af, -6.76894063440646),
        ),
    )
}

/// `evaluate_QNMfit_re_l3m2lp3`: real part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn re_l3m2lp3(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    fma(
        -0.005319515989331159,
        x6,
        fma(
            -0.0021343381708933025,
            x5,
            fma(
                -0.09593447799423722,
                x4,
                fma(
                    0.02244848864087732,
                    x3,
                    fma(1.1906393634562715, x2, fma(-2.107852425643677, af, 1.0)),
                ),
            ),
        ),
    ) / fma(
        -0.005471749827809503,
        x6,
        fma(
            -0.08910191596778137,
            x4,
            fma(1.2043484690080966, x2, fma(-2.1078515887706324, af, 1.0)),
        ),
    )
}

/// `evaluate_QNMfit_im_l2m2lp2`: imaginary part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn im_l2m2lp2(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    af * fma(
        0.024544862718252784,
        x5,
        fma(
            0.054714637311702986,
            x4,
            fma(
                0.0724525431346487,
                x3,
                fma(
                    -0.05898102880105067,
                    x2,
                    fma(-0.47531267226013896, af, 0.3826673013161342),
                ),
            ),
        ),
    ) / fma(
        1.0,
        x6,
        fma(
            -4.152310472191899,
            x4,
            fma(
                -27.99036444363243,
                x2,
                fma(69.82140084545878, af, -38.70835035062785),
            ),
        ),
    )
}

/// `evaluate_QNMfit_im_l2m2lp3`: imaginary part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn im_l2m2lp3(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    af * fma(
        -0.9819458637589314,
        x5,
        fma(
            11.540237244397877,
            x4,
            fma(
                1.721453011852496,
                x3,
                fma(
                    -82.95323455645332,
                    x2,
                    fma(120.43773704442333, af, -49.7688437256778),
                ),
            ),
        ),
    ) / fma(
        1.0,
        x6,
        fma(
            -377.7822297815406,
            x4,
            fma(
                3825.6742092829054,
                x2,
                fma(-6305.619505422591, af, 2858.5790831181725),
            ),
        ),
    )
}

/// `evaluate_QNMfit_im_l3m2lp2`: imaginary part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn im_l3m2lp2(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    af * fma(
        -0.041787576033810676,
        x5,
        fma(
            -0.4236810894599512,
            x4,
            fma(
                -0.08965314412106505,
                x3,
                fma(
                    -1.0115343326360486,
                    x2,
                    fma(4.436434016918535, af, -2.8704762147145533),
                ),
            ),
        ),
    ) / fma(
        1.0,
        x6,
        fma(
            -25.14197656531123,
            x4,
            fma(
                -76.68544453077854,
                x2,
                fma(272.362882450877, af, -171.80908957903395),
            ),
        ),
    )
}

/// `evaluate_QNMfit_im_l3m2lp3`: imaginary part of the spheroidal-spherical
/// mixing coefficient.
pub(crate) fn im_l3m2lp3(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let x6 = x3 * x3;
    af * fma(
        -0.0910871567367674,
        x5,
        fma(
            -3.159763242921214,
            x4,
            fma(
                1.9308599142669403,
                x3,
                fma(
                    18.26221675782779,
                    x2,
                    fma(-29.398484595717147, af, 12.45701482868677),
                ),
            ),
        ),
    ) / fma(
        1.0,
        x6,
        fma(
            -69.3840921447381,
            x4,
            fma(
                538.3888932415709,
                x2,
                fma(-815.4349339779621, af, 345.52914639836257),
            ),
        ),
    )
}
