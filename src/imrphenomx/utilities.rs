//! Utility functions for the IMRPhenomX family, `LALSimIMRPhenomXUtilities.c`.
//!
//! These are the fits for characteristic frequencies, the final state and
//! the time alignment, together with the spin parameterisations. All take
//! the symmetric mass ratio and the aligned spins with the convention
//! `m1 >= m2`.
//!
//! The spin parameterisations and the final mass and spin fits are written
//! with the fused multiply-adds that clang emits for the arm64 builds of
//! LALSimulation (see [`crate::spintaylor::fp`]): they set the ringdown
//! frequencies, which amplify last-bit differences, and the time step of
//! the SpinTaylor evolution used by the precessing models.

use std::f64::consts::PI;

use crate::spintaylor::fp::{fma, fms, fnms};

/// `sqrt(1 - 4 eta)` as contracted by clang.
#[inline]
fn delta_of_eta(eta: f64) -> f64 {
    fnms(4.0, eta, 1.0).sqrt()
}

/// `XLALSimIMRPhenomXsign`: sign of `x`, zero for zero.
pub(crate) fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// `XLALSimIMRPhenomXchiEff`: effective aligned spin.
pub(crate) fn chi_eff(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    let delta = delta_of_eta(eta);
    let mm1 = 0.5 * (1.0 + delta);
    let mm2 = 0.5 * (1.0 - delta);
    fma(mm1, chi1l, mm2 * chi2l)
}

/// `XLALSimIMRPhenomXchiPNHat`: normalised PN reduced spin parameter.
pub(crate) fn chi_pn_hat(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    let delta = delta_of_eta(eta);
    let mm1 = 0.5 * (1.0 + delta);
    let mm2 = 0.5 * (1.0 - delta);
    let chi_eff = fma(mm1, chi1l, mm2 * chi2l);
    fnms((38.0 / 113.0) * eta, chi1l + chi2l, chi_eff) / (1.0 - (76.0 * eta / 113.0))
}

/// `XLALSimIMRPhenomXSTotR`: total spin normalised to [-1, 1].
pub(crate) fn s_tot_r(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    let delta = delta_of_eta(eta);
    let m1 = 0.5 * (1.0 + delta);
    let m2 = 0.5 * (1.0 - delta);
    let m1s = m1 * m1;
    let m2s = m2 * m2;
    fma(m1s, chi1l, m2s * chi2l) / (m1s + m2s)
}

/// `XLALSimIMRPhenomXdchi`: spin difference.
pub(crate) fn dchi(chi1l: f64, chi2l: f64) -> f64 {
    chi1l - chi2l
}

/// `XLALSimIMRPhenomXfMECO`: hybrid minimum energy circular orbit frequency
/// (Cabero et al. 2017), in units of the total mass.
pub(crate) fn f_meco(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let delta = delta_of_eta(eta);
    let s = chi_pn_hat(eta, chi1l, chi2l);
    let s2 = s * s;
    let s3 = s2 * s;
    let dchi = chi1l - chi2l;
    let dchi2 = dchi * dchi;

    let no_spin = fnms(
        0.00006693930988501673,
        eta3,
        fma(
            0.003940354686136861,
            eta2,
            fma(0.0077903147004616865, eta, 0.018744340279608845),
        ),
    ) / fnms(0.10423384680638834, eta, 1.);

    let p4 = fnms(
        0.0230825153005985,
        s2,
        fma(0.022386313074011715, s, -0.0006807631931297156),
    );
    let p2 = fnms(
        0.00038216081981505285,
        s2,
        fnms(0.000010021140796150737, s, 0.00036556167661117023),
    );
    let p1 = fnms(
        0.00035182990586857726,
        s2,
        fnms(0.00001049013062611254, s, 0.00024422562796266645),
    );
    let p3 = fma(
        4.038390455349854e-6,
        s2,
        fma(0.000030679548774047616, s, -0.0005418851224505745),
    );
    let mut e = fnms(0.00002585252361022052, s, 0.00027180386951683135);
    e = fma(eta4, p4, e);
    e = fma(eta2, p2, e);
    e = fma(eta, p1, e);
    e = fma(eta3, p3, e);
    e = fnms(0.00007547517256664526, s2, e);
    let d1 = fma(
        0.005146809717492324,
        eta3,
        fma(
            0.025066696514373803,
            eta2,
            fma(
                1.4861197211952053,
                eta4,
                fnms(0.012429476486138982, eta, -0.014590539285641243),
            ),
        ),
    );
    let d2 = fnms(
        0.001956646166089053,
        eta3,
        fnms(
            0.019641378027236502,
            eta2,
            fnms(
                2.551566872093786,
                eta4,
                fnms(0.02876774751921441, eta, -0.0058684526275074025),
            ),
        ),
    );
    let d3 = fnms(
        0.00767768214056772,
        eta3,
        fma(
            0.012622225233586283,
            eta2,
            fma(
                1.,
                eta4,
                fma(0.014176504653145768, eta, 0.003507640638496499),
            ),
        ),
    );
    let eq_spin = (s * e) / fma(d3, s3, fma(d2, s2, fma(d1, s, 0.026666543809890402)));

    let inner = fnms(
        0.0011713735642446144,
        s,
        fma(
            0.08064665214195679,
            eta2,
            eta * fnms(0.005746537021035632, s, -0.028476219509487793),
        ),
    );
    let uneq_spin = fma(
        dchi2 * fma(0.000016343732281057392, eta, 0.00034375176678815234),
        eta2,
        dchi * delta * eta * inner,
    );

    no_spin + eq_spin + uneq_spin
}

/// `XLALSimIMRPhenomXfISCO`: Kerr innermost stable circular orbit
/// frequency for final spin `chif`, in units of the total mass.
pub(crate) fn f_isco(chif: f64) -> f64 {
    let mut z1 = 1.0 + (1.0 - chif * chif).cbrt() * ((1.0 + chif).cbrt() + (1.0 - chif).cbrt());
    if z1 > 3.0 {
        // Finite precision may give Z1 > 3, but this can not happen.
        z1 = 3.0;
    }
    let z2 = (3.0 * chif * chif + z1 * z1).sqrt();
    let r_isco = 3.0 + z2 - sign(chif) * ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).sqrt();
    let r_isco_sq = r_isco.sqrt();
    let r_isco_3o2 = r_isco_sq * r_isco_sq * r_isco_sq;
    let omega_isco = 1.0 / (r_isco_3o2 + chif);
    omega_isco / PI
}

/// Shared body of `XLALSimIMRPhenomXErad2017` and
/// `XLALSimIMRPhenomXFinalMass2017` (Jimenez-Forteza et al. 2017).
fn erad_2017_terms(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    let delta = delta_of_eta(eta);
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let s = s_tot_r(eta, chi1l, chi2l);
    let s2 = s * s;
    let s3 = s2 * s;
    let dchi = chi1l - chi2l;
    let dchi2 = dchi * dchi;

    let no_spin_terms = || {
        fma(
            3.145145224278187,
            eta4,
            fnms(
                0.84667563764404,
                eta3,
                fma(0.057190958417936644, eta, 0.5609904135313374 * eta2),
            ),
        )
    };
    let no_spin = no_spin_terms();
    // Because of the way this is written, we need to subtract the noSpin term.
    let p1 = fma(
        5.49074464410971,
        eta2,
        fnms(1.1387311580238488, eta, -0.13084389181783257),
    );
    let p2 = fma(2.176667900182948, eta2, -0.17762802148331427);
    let p3 = fnms(
        10.023747993978121,
        eta2,
        fma(4.952698546796005, eta, -0.6320191645391563),
    );
    let q = fma(
        4.274567337924067,
        eta2,
        fma(0.367620218664352, eta, -0.9919475346968611),
    );
    let mut eq_spin =
        (no_spin_terms() * fma(p3, s3, fma(p2, s2, fma(p1, s, 1.0)))) / fma(q, s, 1.0);
    eq_spin -= no_spin;
    let uneq_spin = fnms(
        ((0.01978238971523653 * dchi) * delta) * fnms(4.91667749015812, eta, 1.0) * eta,
        s,
        fma(
            ((-0.09803730445895877 * dchi) * delta) * fnms(3.2283713377939134, eta, 1.0),
            eta2,
            (0.01118530335431078 * dchi2) * eta3,
        ),
    );
    no_spin + eq_spin + uneq_spin
}

/// `XLALSimIMRPhenomXErad2017`: radiated energy in units of the total mass.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn erad_2017(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    erad_2017_terms(eta, chi1l, chi2l)
}

/// `XLALSimIMRPhenomXFinalMass2017`: remnant mass in units of the total mass.
pub(crate) fn final_mass_2017(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    1.0 - erad_2017_terms(eta, chi1l, chi2l)
}

/// `XLALSimIMRPhenomXFinalSpin2017`: dimensionless remnant spin
/// (Jimenez-Forteza et al. 2017).
pub(crate) fn final_spin_2017(eta: f64, chi1l: f64, chi2l: f64) -> f64 {
    let delta = delta_of_eta(eta);
    let m1 = 0.5 * (1.0 + delta);
    let m2 = 0.5 * (1.0 - delta);
    let m1_sq = m1 * m1;
    let m2_sq = m2 * m2;
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let s = s_tot_r(eta, chi1l, chi2l);
    let s2 = s * s;
    let s3 = s2 * s;
    let dchi = chi1l - chi2l;
    let dchi2 = dchi * dchi;

    let no_spin = fnms(
        12.333573402277912,
        eta3,
        fma(3.4641016151377544, eta, 20.0830030082033 * eta2),
    ) / fma(7.2388440419467335, eta, 1.0);
    let p1 = fma(
        1.668810429851045,
        eta3,
        fms(-0.8561951310209386, eta, 0.09939065676370885 * eta2),
    );
    let p2 = fma(
        3.4768263932898678,
        eta3,
        fms(0.5881660363307388, eta, 2.149269067519131 * eta2),
    );
    let p3 = fma(
        1.9595643107593743,
        eta3,
        fms(0.142443244743048, eta, 0.9598353840147513 * eta2),
    );
    let q = fnms(
        9.710576749140989,
        eta3,
        fma(2.3191363426522633, eta, -0.9142232693081653),
    );
    let eq_spin = fma(
        m1_sq + m2_sq,
        s,
        fma(p3, s3, fma(p1, s, p2 * s2)) / fma(q, s, 1.0),
    );
    let uneq_spin = fma(
        ((2.3170397514509933 * dchi) * delta) * fnms(3.2624649875884852, eta, 1.0) * eta3,
        s,
        fms(
            ((0.3223660562764661 * dchi) * delta) * fma(9.332575956437443, eta, 1.0),
            eta2,
            (0.059808322561702126 * dchi2) * eta3,
        ),
    );
    no_spin + eq_spin + uneq_spin
}

/// `XLALSimIMRPhenomXLinb`: fit of the linear-in-frequency phase term
/// (time shift) used to align the model to the calibration hybrids.
///
/// `s` is `STotR` and `delta` the PN symmetry parameter of the waveform.
pub(crate) fn linb(eta: f64, s: f64, dchi: f64, delta: f64) -> f64 {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta3 * eta2;
    let eta6 = eta4 * eta2;
    let s2 = s * s;
    let s3 = s2 * s;
    let s4 = s3 * s;
    let no_spin = 3155.1635543201924 + 1257.9949740608242 * eta - 32243.28428870599 * eta2
        + 347213.65466875216 * eta3
        - 1.9223851649491738e6 * eta4
        + 5.3035911346921865e6 * eta5
        - 5.789128656876938e6 * eta6;
    let eq_spin = (-24.181508118588667 + 115.49264174560281 * eta - 380.19778216022763 * eta2) * s
        + (24.72585609641552 - 328.3762360751952 * eta + 725.6024119989094 * eta2) * s2
        + (23.404604124552 - 646.3410199799737 * eta + 1941.8836639529036 * eta2) * s3
        + (-12.814828278938885 - 325.92980012408367 * eta + 1320.102640190539 * eta2) * s4;
    let uneq_spin = -148.17317525117338 * dchi * delta * eta2;
    no_spin + eq_spin + uneq_spin
}

/// `XLALSimIMRPhenomXPsi4ToStrain`: fit of the time difference between the
/// peaks of the strain and of psi4.
///
/// LAL computes the powers with `pow`, reproduced here with `powf`.
pub(crate) fn psi4_to_strain(eta: f64, s: f64, dchi: f64) -> f64 {
    let eta2 = eta.powf(2.0);
    let eta3 = eta.powf(3.0);
    let eta4 = eta.powf(4.0);
    let s2 = s.powf(2.0);
    let s3 = s.powf(3.0);
    let s4 = s.powf(4.0);
    let no_spin = 13.39320482758057 - 175.42481512989315 * eta + 2097.425116152503 * eta2
        - 9862.84178637907 * eta3
        + 16026.897939722587 * eta4;
    let eq_spin = (4.7895602776763 - 163.04871764530466 * eta + 609.5575850476959 * eta2) * s
        + (1.3934428041390161 - 97.51812681228478 * eta + 376.9200932531847 * eta2) * s2
        + (15.649521097877374 + 137.33317057388916 * eta - 755.9566456906406 * eta2) * s3
        + (13.097315867845788 + 149.30405703643288 * eta - 764.5242164872267 * eta2) * s4;
    let uneq_spin = 105.37711654943146 * dchi * (1. - 4. * eta).sqrt() * eta2;
    no_spin + eq_spin + uneq_spin
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imrphenomx::qnm::{fdamp22, fring22};

    /// (eta, chi1, chi2, [final spin, final mass, erad, fring22, fdamp22, fMECO,
    /// fISCO, chiPNHat, STotR, chiEff, linb, psi4ToStrain]) from the
    /// lalsimulation 6.2.0 python bindings (`SimIMRPhenomXFinalSpin2017`, ...).
    #[rustfmt::skip]
    const LAL_VALUES: [(f64, f64, f64, [f64; 12]); 5] = [
        (0.25, 0.0, 0.0, [0.6863698793893069, 0.9515839428556587, 0.04841605714434129, 0.08382415886824093, 0.012938159666496444, 0.021497332614512636, 0.04478020279587562, 0.0, 0.0, 0.0, 3136.2800623853086, 9.124237969507114]),
        (0.2222222222222222, 0.3, -0.4, [0.694625703151086, 0.9574380880216553, 0.04256191197834466, 0.08439148246730528, 0.012890844381762334, 0.022190105694075733, 0.04539055868812298, 0.08716763005780352, 0.1600000000000001, 0.06666666666666671, 3136.7914141719752, 9.835000601651613]),
        (0.16, 0.8, 0.5, [0.8825833580508831, 0.9540333487125494, 0.045966651287450616, 0.10386758091612643, 0.010743850486461898, 0.03508567717116395, 0.06798706542073071, 0.7508528361761206, 0.7823529411764706, 0.7400000000000001, 3107.7445127366636, 17.637368968733767]),
        (0.0826446280991736, -0.7, 0.2, [-0.2762495494338788, 0.9943840766641855, 0.0056159233358144675, 0.0546782483486677, 0.01419467195355803, 0.014267076993033338, 0.017928474290187783, -0.6398513126306822, -0.6910891089108911, -0.618181818181818, 3176.2067274551796, 6.4461725554068]),
        (0.04535147392290249, 0.99, -0.99, [0.9872121260217379, 0.9820155475349732, 0.01798445246502673, 0.13627477149156209, 0.005159383168537771, 0.03593434480002259, 0.11234776383668363, 0.9238948109204561, 0.9850623441396509, 0.8957142857142857, 3136.8770722836057, 42.42214053271711]),
    ];

    /// LALSimulation is compiled with floating point contraction (fused
    /// multiply-adds), which this port does not reproduce in general, so
    /// fits agree to a few units in the last place rather than bitwise.
    fn assert_close(name: &str, value: f64, expected: f64, ulps: f64) {
        assert!(
            (value - expected).abs() <= ulps * f64::EPSILON * expected.abs(),
            "{name}: {value:e} != {expected:e} (lalsimulation)"
        );
    }

    #[test]
    fn fits_match_lalsimulation() {
        for (eta, chi1, chi2, expected) in LAL_VALUES {
            let s = s_tot_r(eta, chi1, chi2);
            let delta = (1.0 - 4.0 * eta).sqrt();
            let lal_af = expected[0];
            let values = [
                ("final spin", final_spin_2017(eta, chi1, chi2)),
                ("final mass", final_mass_2017(eta, chi1, chi2)),
                ("erad", erad_2017(eta, chi1, chi2)),
                // The QNM fits are evaluated at LAL's final spin: their
                // denominators amplify differences in the input by up to ~1e4.
                ("fring22", fring22(lal_af)),
                ("fdamp22", fdamp22(lal_af)),
                ("fMECO", f_meco(eta, chi1, chi2)),
                ("fISCO", f_isco(lal_af)),
                ("chiPNHat", chi_pn_hat(eta, chi1, chi2)),
                ("STotR", s),
                ("chiEff", chi_eff(eta, chi1, chi2)),
                ("linb", linb(eta, s, dchi(chi1, chi2), delta)),
                ("psi4ToStrain", psi4_to_strain(eta, s, dchi(chi1, chi2))),
            ];
            for ((name, value), expected) in values.iter().zip(expected) {
                // The QNM fits use fused multiply-adds like LAL's build and
                // agree exactly.
                // The QNM fits use fused multiply-adds like LAL's build and
                // agree exactly; so do the final-state fits and the spin
                // parameterisations where the target fuses like LAL's arm64
                // build (see `crate::spintaylor::fp`).
                let qnm = name.starts_with("fring") || name.starts_with("fdamp");
                let final_state = [
                    "final spin",
                    "final mass",
                    "erad",
                    "fMECO",
                    "chiPNHat",
                    "STotR",
                    "chiEff",
                ];
                let exact = qnm
                    || (crate::spintaylor::fp::FUSED_MULTIPLY_ADD && final_state.contains(name));
                let ulps = if exact { 0.0 } else { 8.0 };
                assert_close(name, *value, expected, ulps);
            }
        }
    }
}
