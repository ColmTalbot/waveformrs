// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The amplitude of the IMRPhenomXHM higher modes in the original release,
//! `PhenomXHMReleaseVersion = 122019` (arXiv:2001.10914), as used by
//! IMRPhenomXPNR: the parts of `IMRPhenomXHM_SetHMWaveformVariables`,
//! `IMRPhenomXHM_GetAmplitudeCoefficients` and the amplitude ansatzes of
//! `LALSimIMRPhenomXHM_internals.c` (and the files it includes) that differ
//! from the 122022 release.
//!
//! The release uses the 122018 parameter-space fits
//! ([`crate::imrphenomx::hm_fits_122018`]), three inspiral collocation points
//! on top of the re-expanded PN amplitude (or the TaylorT4 based amplitude of
//! the (2,1) mode), an intermediate amplitude that is the inverse of a
//! polynomial ([`crate::imrphenomx::hm_poly_122019`]) whose order is lowered
//! when collocation points are vetoed or when it crosses zero, two
//! intermediate regions for extreme mass ratios, and the original ringdown
//! ansatz (`IMRPhenomXHMRingdownAmpVersion = 0`). The phase is the same as in
//! the 122022 release, except that the (3,2) ringdown phase follows from the
//! spheroidal-to-spherical rotation of the amplitude.
//!
//! As in the rest of the crate the operations follow LAL (with fused
//! multiply-adds where the hardware has them, see [`fma`]); the complex
//! square and cube roots are those of [`crate::complex`].

use num_complex::Complex;

use crate::complex::{csqrt, principal_cbrt};
use crate::imrphenomx::hm_fits_122018 as fits;
use crate::imrphenomx::hm_internals::{HMAmpCoefficients, HMWaveformStruct};
use crate::imrphenomx::hm_poly_122019::{self as poly, PolyInputs};
use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct, POWERS_OF_LALPI};
use crate::imrphenomx::poly_roots;

/// `PhenomXHMReleaseVersion` of the original IMRPhenomXHM.
pub(crate) const RELEASE_122019: u32 = 122019;

/// The default `PhenomXHMReleaseVersion` (recalibrated amplitudes).
pub(crate) const RELEASE_122022: u32 = 122022;

/// Threshold of the imaginary part below which `CrossZeroP3`, `CrossZeroP4`
/// and `CrossZeroP5` consider a root real, `pow(10., -15)`.
const ROOT_THRESHOLD: f64 = 1e-15;

/// Deviations of the (3,3) mode in the PNR co-precessing frame,
/// `PNR_DEV_PARAMETER` and `MU1` ... `MU4` of
/// `IMRPhenomXHMWaveformStruct`. They are zero unless the tuned
/// co-precessing model is used (`PNRUseTunedCoprec33`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Pnr33Deviations {
    /// `PNR_DEV_PARAMETER` (already multiplied by `delta`).
    pub(crate) dev_parameter: f64,
    /// `MU1`, `MU2` (intermediate collocation points), `MU3` (`alambda`)
    /// and `MU4` (`sigma`).
    pub(crate) mu: [f64; 4],
}

/// The flags of `IMRPhenomXHMWaveformStruct` that only the 122019 release
/// uses, set in `IMRPhenomXHM_SetHMWaveformVariables`, and the
/// reconstruction versions as modified by the vetoes in
/// `IMRPhenomXHM_GetAmplitudeCoefficients`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Flags {
    /// `AmpEMR`: two intermediate regions (q > 70 and chi1 <= 0.9).
    pub(crate) amp_emr: bool,
    /// `useFAmpPN`: TaylorT4 based inspiral amplitude of the (2,1) (q <~ 40).
    pub(crate) use_famp_pn: bool,
    /// `InspiralAmpVeto`, `IntermediateAmpVeto`, `RingdownAmpVeto`.
    pub(crate) inspiral_amp_veto: bool,
    pub(crate) intermediate_amp_veto: bool,
    pub(crate) ringdown_amp_veto: bool,
    /// `IMRPhenomXHMInspiralAmpVersion`: number of inspiral collocation
    /// points used (3 before the vetoes).
    pub(crate) inspiral_amp_version: i32,
    /// `IMRPhenomXHMIntermediateAmpVersion`: reconstruction of the
    /// intermediate amplitude (2 before the coefficients are computed).
    pub(crate) intermediate_amp_version: i32,
    /// PNR deviations of the (3,3) mode.
    pub(crate) pnr: Pnr33Deviations,
    /// Not in LAL: evaluate the ringdown amplitude with one `exp`, `f^(-lc)
    /// exp(-df lambda / dfd) = exp(-df lambda / dfd - lc ln f)`, with the
    /// logarithm of the powers (from a grid), instead of `pow` and `exp`
    /// (the grid paths; equal to rounding).
    pub(crate) grid_ringdown: bool,
}

/// The release-specific part of `IMRPhenomXHM_SetHMWaveformVariables` for
/// `PhenomXHMReleaseVersion = 122019`.
pub(crate) fn set_waveform_variables(wfhm: &mut HMWaveformStruct, wf: &WaveformStruct) {
    let mut flags = Flags {
        inspiral_amp_version: 3,
        intermediate_amp_version: 2,
        ..Flags::default()
    };
    // For q > 70 and chi1 <= 0.9 the amplitude drops more steeply at the end
    // of the inspiral: use two intermediate regions.
    if wf.eta < 0.013886133703630232 && wf.chi1l <= 0.9 {
        flags.amp_emr = true;
    }
    match wfhm.mode_tag {
        21 => {
            if wf.q < 8. {
                flags.inspiral_amp_veto = true;
                flags.intermediate_amp_veto = true;
            }
            flags.ringdown_amp_veto = true;
            // TaylorT4 based inspiral for q <~ 40 (LAL's comment says
            // "for EMR (q>40.)", the condition is the opposite).
            if wf.eta >= 0.0237954 {
                flags.use_famp_pn = true;
            }
        }
        32 => flags.ringdown_amp_veto = true,
        // The (3,3) of equal masses and spins vanishes as in both releases.
        _ => {}
    }
    wfhm.ringdown_amp_version = 0;
    wfhm.v122019 = flags;
}

/// Coefficients of the TaylorT4 based inspiral amplitude of the (2,1),
/// `Get21PNAmplitudeCoefficients`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Pn21Coefficients {
    /// `PNTDfactor`
    pub(crate) pntd_factor: f64,
    /// The complex time-domain PN coefficients `x05`, `x1`, `x15`, `x2`,
    /// `x25` and `x3`.
    pub(crate) x: [Complex<f64>; 6],
    pub(crate) xdot5: f64,
    pub(crate) xdot6: f64,
    pub(crate) xdot65: f64,
    pub(crate) xdot7: f64,
    pub(crate) xdot75: f64,
    pub(crate) xdot8: f64,
    pub(crate) xdot8_log: f64,
    pub(crate) xdot85: f64,
    pub(crate) log2pi_two_thirds: f64,
}

/// The amplitude coefficients of `IMRPhenomXHMAmpCoefficients` that only
/// the 122019 release uses.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct AmpCoefficients {
    pub(crate) pn21: Pn21Coefficients,
    /// Pseudo-PN coefficients `rho1`, `rho2`, `rho3` of the inspiral.
    pub(crate) rho: [f64; 3],
    /// `RDCoefficient[0..4]`: `alambda`, `lambda`, `sigma` and `lc`.
    pub(crate) rd: [f64; 4],
    /// `useInspAnsatzRingdown`
    pub(crate) use_insp_ansatz_ringdown: bool,
    /// The inputs of the intermediate reconstruction (`v1` ... `d4`).
    pub(crate) poly_inputs: PolyInputs,
    /// `delta0` ... `delta5` of the intermediate polynomial.
    pub(crate) delta: [f64; 6],
    /// `alpha0` ... `alpha4` of the first intermediate region (`AmpEMR`).
    pub(crate) alpha: [f64; 5],
    /// End of the first intermediate region, `fAmpMatchInt12`.
    pub(crate) f_amp_match_int12: f64,
    /// `InterAmpPolOrder`: order of the intermediate polynomial.
    pub(crate) inter_amp_pol_order: i32,
}

/* ---------------------------------------------------------------------- */
/*                          Cutting frequencies                            */
/* ---------------------------------------------------------------------- */

/// `IMRPhenomXHM_Amplitude_fcutInsp` (version 122018).
fn amplitude_fcut_insp(wfhm: &HMWaveformStruct, wf: &WaveformStruct) -> f64 {
    let f_meco = wfhm.f_meco_lm;
    let emm = wfhm.emm as f64;
    let eta = wf.eta;
    let chi1 = wf.chi1l;
    let fcut_emr = 1.25
        * emm
        * (fma(
            0.000316707064291237,
            chi1 * chi1,
            fma(-0.0000858396080377194, chi1, 0.011671068725758493),
        ) * fma(6.2873167352395125, eta, 0.8447212540381764))
        / fma(-0.9977728883419751, chi1, 1.2857082764038923);
    let chieff = wf.chi_eff;
    let f_isco = wf.f_isco * emm * 0.5;
    let df = (f_isco - f_meco).abs();
    match wfhm.mode_tag {
        21 => {
            if eta < 0.023795359904818562 {
                fcut_emr
            } else {
                let factor = fma(-(5. / 6. * chieff), chieff, fma(-0.235, chieff, 0.75));
                fma(factor, df, f_meco)
            }
        }
        33 => {
            if eta < 0.04535147392290249 {
                fcut_emr
            } else {
                let factor = fma(-(5. / 6.), chieff, fma(-0.235, chieff, 0.75));
                fma(factor, df, f_meco)
            }
        }
        32 => {
            if eta < 0.04535147392290249 {
                fcut_emr
            } else {
                let fcut = fma(fma(-0.235, chieff.abs(), 0.75), df, f_meco);
                fcut * wfhm.f_ring / wf.f_ring
            }
        }
        _ => {
            if eta < 0.04535147392290249 {
                fcut_emr
            } else {
                fma(fma(-0.235, chieff, 0.75), df, f_meco)
            }
        }
    }
}

/// `IMRPhenomXHM_Amplitude_fcutRD` (version 122018).
fn amplitude_fcut_rd(wfhm: &HMWaveformStruct, wf: &WaveformStruct) -> f64 {
    let fring = wfhm.f_ring;
    let fdamp = wfhm.f_damp;
    let eta = wf.eta;
    let chi1 = wf.chi1l;
    match wfhm.mode_tag {
        21 => 0.75 * fring,
        33 => 0.95 * fring,
        32 => {
            let f_rd22 = wf.f_ring;
            let (c, r): (f64, f64) = (0.5, 5.);
            let mut fcut = if eta < 0.0453515 {
                // A smooth step between fRING (negative spins) and fRD22
                // (positive spins).
                fma(fring, (c * r).exp(), f_rd22 * (r * chi1).exp())
                    / ((c * r).exp() + (r * chi1).exp())
                    - fdamp
            } else {
                f_rd22
            };
            if 0.02126654064272212 < eta && eta < 0.12244897959183673 && chi1 > 0.95 {
                fcut = fma(-2., fdamp, fring);
            }
            fcut
        }
        _ => 0.9 * fring,
    }
}

/* ---------------------------------------------------------------------- */
/*                              Inspiral                                   */
/* ---------------------------------------------------------------------- */

/// `Get21PNAmplitudeCoefficients`. The complex coefficients follow C's
/// complex arithmetic: `I * x` is `(0 * x, x)`, products and sums are
/// evaluated componentwise without contraction.
fn pn21_coefficients(wf: &WaveformStruct) -> Pn21Coefficients {
    let lalpi = &*POWERS_OF_LALPI;
    let p2 = UsefulPowers::new(2.);
    let logof2 = p2.log;
    let delta = wf.delta;
    let eta = wf.eta;
    let m1 = wf.m1;
    let m2 = wf.m2;
    let chi1 = wf.chi1l;
    let chi2 = wf.chi2l;
    let chi_s = (chi1 + chi2) * 0.5;
    let chi_a = (chi1 - chi2) * 0.5;
    let sc = fma(m1 * m1, chi1, m2 * m2 * chi2);
    let sigmac = fma(m2, chi2, -(m1 * chi1));

    // `I` is the single-precision `_Complex_I`, promoted to double when
    // combined with a double; `I * delta` is `(0 * delta, delta)`.
    let i_delta = Complex::new(0.0 * delta, 1.0 * delta);
    let scale = |z: Complex<f64>, a: f64| Complex::new(z.re * a, z.im * a);
    let divide = |z: Complex<f64>, a: f64| Complex::new(z.re / a, z.im / a);
    let mul = |a: Complex<f64>, b: Complex<f64>| {
        Complex::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
    };

    // I*delta/3*powers_of_2d1.one_third*powers_of_lalpiHM.one_third
    let x05 = scale(scale(divide(i_delta, 3.), p2.one_third), lalpi.one_third);
    // -I*0.5*(chiA + chiS*delta)*powers_of_2d1.two_thirds*powers_of_lalpiHM.two_thirds
    let minus_i_half = Complex::new(-0.0 * 0.5, -1.0 * 0.5);
    let x1 = scale(
        scale(scale(minus_i_half, fma(chi_s, delta, chi_a)), p2.two_thirds),
        lalpi.two_thirds,
    );
    // I*delta*(-17./28 + 5*eta/7.)/3.*powers_of_2d1.itself*powers_of_lalpiHM.itself
    let x15 = scale(
        scale(
            divide(scale(i_delta, -17. / 28. + 5. * eta / 7.), 3.),
            p2.itself,
        ),
        lalpi.itself,
    );
    // ( I*(-43./21*delta*Sc + (-79.+139*eta)/42.*Sigmac)
    //   + I/3*delta*(powers_of_lalpiHM.itself + I*(-0.5-2*logof2)) )
    //   *powers_of_2d1.four_thirds*powers_of_lalpiHM.four_thirds
    let a = fma(-43. / 21. * delta, sc, fma(139., eta, -79.) / 42. * sigmac);
    let i_a = Complex::new(0.0 * a, 1.0 * a);
    let w = fma(-2., logof2, -0.5);
    let pi_plus_iw = Complex::new(lalpi.itself + 0.0 * w, 1.0 * w);
    // (LAL's `I/3` is a single-precision `complex float` divided by an int,
    // so 1/3 to ~3e-8; it is 1/3 in double precision here, see
    // docs/lal-quirks.md.)
    let i_third = 1.0 / 3.0;
    let i_third_delta = Complex::new(0.0 * delta, i_third * delta);
    let product = mul(i_third_delta, pi_plus_iw);
    let x2 = scale(
        scale(
            Complex::new(i_a.re + product.re, i_a.im + product.im),
            p2.four_thirds,
        ),
        lalpi.four_thirds,
    );
    // I*delta*((-43-509*eta)/126. + 79*eta*eta/168.)/3.*powers_of_2d1.five_thirds*powers_of_lalpiHM.five_thirds
    let t25 = fma(-509., eta, -43.) / 126. + 79. * eta * eta / 168.;
    let x25 = scale(
        scale(divide(scale(i_delta, t25), 3.), p2.five_thirds),
        lalpi.five_thirds,
    );
    // I*delta*( (-17. + 6.*eta)/28.*powers_of_lalpiHM.itself
    //   + I*(17/56. + eta*(-353/28. - 3.*logof2/7.) + 17.*logof2/14.) )/3.*powers_of_2d1.two*powers_of_lalpiHM.two
    let u = fma(6., eta, -17.) / 28. * lalpi.itself;
    let v = fma(eta, -353. / 28. - 3. * logof2 / 7., 17. / 56.) + 17. * logof2 / 14.;
    let u_plus_iv = Complex::new(u + 0.0 * v, 1.0 * v);
    let x3 = scale(
        scale(divide(mul(i_delta, u_plus_iv), 3.), p2.two),
        lalpi.two,
    );

    Pn21Coefficients {
        pntd_factor: 8. * eta * lalpi.two_thirds * p2.two_thirds * lalpi.sqrt / 5f64.sqrt(),
        x: [x05, x1, x15, x2, x25, x3],
        xdot5: xdot5(wf, &p2, lalpi),
        xdot6: xdot6(wf, &p2, lalpi),
        xdot65: xdot65(wf, &p2, lalpi),
        xdot7: xdot7(wf, &p2, lalpi),
        xdot75: xdot75(wf, &p2, lalpi),
        xdot8: xdot8(wf, &p2, lalpi),
        xdot8_log: xdot8_log(wf, &p2, lalpi),
        xdot85: xdot85(wf, &p2, lalpi),
        log2pi_two_thirds: (lalpi.two_thirds * p2.two_thirds).ln(),
    }
}

/// `IMRPhenomXHM_Inspiral_PNAmp_21Ansatz` (strain amplitude without `amp0`):
/// the time-domain PN amplitude times the stationary-phase factor of
/// TaylorT4.
fn pn21_ansatz(powers: &UsefulPowers, c: &Pn21Coefficients) -> f64 {
    let lalpi = &*POWERS_OF_LALPI;
    let two_to_m_one_sixths = 0.8908987181403393;
    let three_to_m_one_second = 0.5773502691896257;
    let x_to_m_one_four = two_to_m_one_sixths * lalpi.m_one_sixth * powers.m_one_sixth;
    let p = [
        powers.one_third,
        powers.two_thirds,
        powers.itself,
        powers.four_thirds,
        powers.five_thirds,
        powers.two,
    ];
    // Complex sum of real * complex products, no contraction.
    let mut re = p[0] * c.x[0].re;
    let mut im = p[0] * c.x[0].im;
    for (pk, xk) in p.iter().zip(&c.x).skip(1) {
        re += pk * xk.re;
        im += pk * xk.im;
    }
    re = re * powers.two_thirds * c.pntd_factor;
    im = im * powers.two_thirds * c.pntd_factor;
    let xdot = xdot_t4(powers, c);
    2. * lalpi.sqrt * three_to_m_one_second * re.hypot(im) * x_to_m_one_four / xdot.sqrt()
}

/// `IMRPhenomXHM_Inspiral_PNAmp_Ansatz` with `InspRescaleFactor = 1`: the PN
/// amplitude divided by the leading order of the 22, `ampNorm f^(-7/6)`.
fn pn_amp_rescaled(powers: &UsefulPowers, wfhm: &HMWaveformStruct, amp: &HMAmpCoefficients) -> f64 {
    if wfhm.v122019.use_famp_pn {
        // RescaleFactor 1
        return pn21_ansatz(powers, &amp.v122019.pn21) / (amp.amp_norm * powers.m_seven_sixths);
    }
    let pn = &amp.pn;
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
    re.hypot(im) * amp.pn_global_factor
}

/// `IMRPhenomXHM_Inspiral_Amp_Ansatz` (frequencies version 122018) during
/// the computation of the coefficients (`InspRescaleFactor = 1`): PN plus
/// pseudo-PN terms, rescaled by the leading order of the 22.
fn inspiral_ansatz_rescaled(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    let rho = &amp.v122019.rho;
    let pn = pn_amp_rescaled(powers, wfhm, amp);
    let mut insp = fma(powers.seven_thirds / amp.fcut_insp_seven_thirds, rho[0], pn);
    insp = fma(
        powers.eight_thirds / amp.fcut_insp_eight_thirds,
        rho[1],
        insp,
    );
    fma(powers.three / amp.fcut_insp_three, rho[2], insp)
}

/// `IMRPhenomXHM_Inspiral_Amp_Ansatz` (`InspRescaleFactor = 0`): the
/// inspiral amplitude in strain units without `amp0`.
#[inline]
pub(crate) fn inspiral_ansatz(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    inspiral_ansatz_rescaled(powers, wfhm, amp) * (powers.m_seven_sixths * wfhm.amp_norm)
}

/// The fourth-order finite-difference derivative with step 1e-9 of
/// `IMRPhenomXHM_Inspiral_Amp_NDAnsatz` and `IMRPhenomXHM_RD_Amp_NDAnsatz`.
fn nd_derivative(f: f64, at: impl Fn(f64) -> f64) -> f64 {
    let df = 10e-10;
    let fun2r = at(fma(2., df, f));
    let funr = at(f + df);
    let funl = at(f - df);
    let fun2l = at(fma(-2., df, f));
    (fma(-8., funl, fma(8., funr, -fun2r)) + fun2l) / (12. * df)
}

/// `WavyPoints`: whether three consecutive values are not monotonic.
fn wavy_points(p1: f64, p2: f64, p3: f64) -> bool {
    (p1 > p2 && p2 < p3) || (p1 < p2 && p2 > p3)
}

/// The pseudo-PN coefficients `rho1`, `rho2`, `rho3` for
/// `IMRPhenomXHMInspiralAmpVersion` collocation points.
fn rho_coefficients(
    version: i32,
    v: [f64; 3],
    fcut_insp: &UsefulPowers,
    f: [&UsefulPowers; 3],
) -> [f64; 3] {
    let [v1, v2, v3] = v;
    let [f1, f2, f3] = f;
    match version {
        3 => [
            rho1_3(v1, v2, v3, fcut_insp, f1, f2, f3),
            rho2_3(v1, v2, v3, fcut_insp, f1, f2, f3),
            rho3_3(v1, v2, v3, fcut_insp, f1, f2, f3),
        ],
        2 => [
            rho1_2(v1, v2, fcut_insp, f1, f2),
            rho2_2(v1, v2, fcut_insp, f1, f2),
            0.,
        ],
        1 => [rho1_1(v1, fcut_insp, f1), 0., 0.],
        0 => [0.; 3],
        _ => panic!("IMRPhenomXHM: inspiral amplitude version {version} is not valid"),
    }
}

/* ---------------------------------------------------------------------- */
/*                              Ringdown                                   */
/* ---------------------------------------------------------------------- */

/// `IMRPhenomXHM_RD_Amp_Ansatz` (version 0): the ringdown amplitude (the
/// spheroidal one for the (3,2)) rescaled by the leading order of the 22.
#[inline]
pub(crate) fn rd_amp_ansatz(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
) -> f64 {
    let rd = &amp.v122019.rd;
    if wfhm.v122019.grid_ringdown {
        let [a0, lambda, sigma, lc] = [rd[0], rd[1], rd[2], rd[3]];
        let (ff, frd, fda) = (powers.itself, wfhm.f_ring, wfhm.f_damp);
        let dfr = ff - frd;
        let dfd = fda * sigma;
        return fda * a0.abs() * sigma * fma(-lc, powers.log, -dfr * lambda / dfd).exp()
            / fma(dfr, dfr, dfd * dfd);
    }
    rd_amp_ansatz0(
        powers.itself,
        wfhm.f_ring,
        wfhm.f_damp,
        rd[0],
        rd[1],
        rd[2],
        rd[3],
    )
}

/// `IMRPhenomXHM_RD_Amp_DAnsatz` (version 0).
fn rd_amp_dansatz(powers: &UsefulPowers, wfhm: &HMWaveformStruct, amp: &HMAmpCoefficients) -> f64 {
    let rd = &amp.v122019.rd;
    rd_amp_dansatz0(
        powers.itself,
        wfhm.f_ring,
        wfhm.f_damp,
        rd[0],
        rd[1],
        rd[2],
        rd[3],
    )
}

/// `IMRPhenomXHM_Ringdown_Amplitude_Veto`: drop both intermediate
/// collocation points if the ringdown amplitude is very low.
fn ringdown_veto(v2: &mut f64, v3: &mut f64, v4: f64, version: &mut i32, wfhm: &HMWaveformStruct) {
    if wfhm.mode_tag == 32 {
        let threshold = 0.1 / wfhm.amp_norm;
        if 1. / v4 < threshold {
            *v2 = 1.;
            *v3 = 1.;
            *version = 101;
        }
    } else {
        let threshold = 0.01 / wfhm.amp_norm;
        if 1. / v4 < threshold {
            *v2 = 1.;
            *v3 = 1.;
            *version = 1032;
        }
    }
}

/* ---------------------------------------------------------------------- */
/*                             Intermediate                                */
/* ---------------------------------------------------------------------- */

/// `IMRPhenomXHM_Intermediate_Amplitude_Veto`: drop too low intermediate
/// collocation points (the values are inverse amplitudes).
fn intermediate_veto(int1: &mut f64, int2: &mut f64, version: &mut i32, amp_norm: f64) {
    let threshold = 0.2 / amp_norm;
    if 1. / *int1 < threshold {
        *int1 = 1.;
        *version = 1042;
        if 1. / *int2 < threshold {
            *int2 = 1.;
            *version = 1032;
        }
    } else if 1. / *int2 < threshold {
        *int2 = 1.;
        *version = 1042;
    }
}

/// `InsideInterval`
fn inside_interval(ftest: f64, fmin: f64, fmax: f64) -> bool {
    ftest >= fmin && ftest <= fmax
}

/// Complex division `(a + ib) / (c + id)`, scaled by a power of two to
/// avoid overflow and underflow (the algorithm of compiler-rt's `__divdc3`).
fn divdc3(a: f64, b: f64, c: f64, d: f64) -> Complex<f64> {
    let mut c = c;
    let mut d = d;
    let mut ilogbw = 0;
    let logbw = c.abs().max(d.abs()).log2().floor();
    if logbw.is_finite() {
        ilogbw = logbw as i32;
        c = scalbn(c, -ilogbw);
        d = scalbn(d, -ilogbw);
    }
    let denom = c * c + d * d;
    let mut x = scalbn((a * c + b * d) / denom, -ilogbw);
    let mut y = scalbn((b * c - a * d) / denom, -ilogbw);
    if x.is_nan() && y.is_nan() {
        // compiler-rt's recovery of infinities and zeros
        if denom == 0.0 && (!a.is_nan() || !b.is_nan()) {
            x = f64::INFINITY.copysign(c) * a;
            y = f64::INFINITY.copysign(c) * b;
        } else if (a.is_infinite() || b.is_infinite()) && c.is_finite() && d.is_finite() {
            let a = (if a.is_infinite() { 1.0f64 } else { 0.0 }).copysign(a);
            let b = (if b.is_infinite() { 1.0f64 } else { 0.0 }).copysign(b);
            x = f64::INFINITY * (a * c + b * d);
            y = f64::INFINITY * (b * c - a * d);
        } else if logbw == f64::INFINITY && logbw > 0.0 && a.is_finite() && b.is_finite() {
            let c = (if c.is_infinite() { 1.0f64 } else { 0.0 }).copysign(c);
            let d = (if d.is_infinite() { 1.0f64 } else { 0.0 }).copysign(d);
            x = 0.0 * (a * c + b * d);
            y = 0.0 * (b * c - a * d);
        }
    }
    Complex::new(x, y)
}

/// `scalbn(x, n)`: `x 2^n`, exact unless the result under- or overflows.
fn scalbn(x: f64, n: i32) -> f64 {
    // Two steps keep the factors representable for |n| <= 2046.
    let half = n / 2;
    x * 2f64.powi(half) * 2f64.powi(n - half)
}

/// Division of a complex number by a real one, `(re / x, im / x)`.
fn cdiv_real(z: Complex<f64>, x: f64) -> Complex<f64> {
    Complex::new(z.re / x, z.im / x)
}

/// Product of a real number and a complex one, `(x re, x im)`.
fn cmul_real(x: f64, z: Complex<f64>) -> Complex<f64> {
    Complex::new(x * z.re, x * z.im)
}

/// Product of two complex numbers.
fn cmul(a: Complex<f64>, b: Complex<f64>) -> Complex<f64> {
    Complex::new(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
}

/// Complex division ([`divdc3`]).
fn cdiv(a: Complex<f64>, b: Complex<f64>) -> Complex<f64> {
    divdc3(a.re, a.im, b.re, b.im)
}

/// Whether a root `f` is real (imaginary part below 1e-15) and inside
/// `[fstart, fend]`.
fn real_root_inside(f: Complex<f64>, fstart: f64, fend: f64) -> bool {
    f.im.abs() < ROOT_THRESHOLD && inside_interval(f.re, fstart, fend)
}

/// `CrossZeroP2`: whether the quadratic polynomial crosses zero in
/// `[fstart, fend]`.
fn cross_zero_p2(a: &[f64; 6], fstart: f64, fend: f64) -> bool {
    let [a0, a1, a2, ..] = *a;
    let discriminant = p2_discriminant(a0, a1, a2);
    if discriminant < 0. {
        return false;
    }
    let f1 = (-a1 - discriminant.sqrt()) / (2. * a2);
    let f2 = (-a1 + discriminant.sqrt()) / (2. * a2);
    inside_interval(f1, fstart, fend) || inside_interval(f2, fstart, fend)
}

/// `CrossZeroP3`: whether the cubic polynomial crosses zero in
/// `[fstart, fend]`, from Cardano's formula in complex arithmetic.
fn cross_zero_p3(a: &[f64; 6], fstart: f64, fend: f64) -> bool {
    let [a0, a1, a2, a3, ..] = *a;
    let q1 = -a2 / (3. * a3);
    let q2 = p3_q2(a1, a2, a3);
    let q3 = p3_q3(a0, a1, a2, a3);
    let discri = p3_discri(q2, q3);
    let root = csqrt(Complex::new(discri, 0.));
    let onethirdpower = principal_cbrt(Complex::new(q3 + root.re, root.im));
    let twothird = 2f64.cbrt();
    let sqrt3 = 3f64.sqrt();
    let z1 = Complex::new(1., sqrt3);
    let z2 = Complex::new(1., -sqrt3);
    let real = |x: f64| Complex::new(x, 0.);
    // q1 - q2*twothird/(3*a3*onethirdpower) + onethirdpower/(3.*a3*twothird)
    let t = cdiv(real(q2 * twothird), cmul_real(3. * a3, onethirdpower));
    let u = cdiv_real(onethirdpower, 3. * a3 * twothird);
    let f1 = Complex::new(q1 - t.re + u.re, -t.im + u.im);
    // q1 + q2*z1/(3*twothird*twothird*a3*onethirdpower) - onethirdpower*z2/(6*twothird*a3)
    let den = cmul_real(3. * twothird * twothird * a3, onethirdpower);
    let root_k = |za: Complex<f64>, zb: Complex<f64>| {
        let t = cdiv(cmul_real(q2, za), den);
        let u = cdiv_real(cmul(onethirdpower, zb), 6. * twothird * a3);
        Complex::new(q1 + t.re - u.re, t.im - u.im)
    };
    let f2 = root_k(z1, z2);
    let f3 = root_k(z2, z1);
    real_root_inside(f1, fstart, fend)
        || real_root_inside(f2, fstart, fend)
        || real_root_inside(f3, fstart, fend)
}

/// `CrossZeroP4`: whether the quartic polynomial crosses zero in
/// `[fstart, fend]`, from Ferrari's formula in complex arithmetic.
///
/// Like LAL, `sqrt(q1 + frac1 + frac2)` in `bigdenom` is the real square
/// root of the real part (the C `sqrt` of a complex argument).
fn cross_zero_p4(a: &[f64; 6], fstart: f64, fend: f64) -> bool {
    let [a0, a1, a2, a3, a4, _] = *a;
    let q0 = -a3 / (4. * a4);
    let q1 = p4_q1(a2, a3, a4);
    let q1b = 2. * q1;
    let q2 = p4_q2(a0, a1, a2, a3, a4);
    let q3 = p4_q3(a0, a1, a2, a3, a4);
    let squareroot = csqrt(Complex::new(p4_sqrt_arg(q2, q3), 0.));
    let onethird = principal_cbrt(Complex::new(q3 + squareroot.re, squareroot.im));
    let twothird = 2f64.cbrt();
    // twothird*q2/(3*a4*onethird)
    let frac1 = cdiv(
        Complex::new(twothird * q2, 0.),
        cmul_real(3. * a4, onethird),
    );
    // onethird/(3*twothird*a4)
    let frac2 = cdiv_real(onethird, 3. * twothird * a4);
    // q1 + frac1 + frac2
    let s = Complex::new(q1 + frac1.re + frac2.re, frac1.im + frac2.im);
    let bigdenom = Complex::new(4. * s.re.sqrt(), 0.);
    let bignum = p4_bignum(a1, a2, a3, a4);
    let sqrt_s = csqrt(s);
    let ratio = cdiv(Complex::new(bignum, 0.), bigdenom);
    // q1b - frac1 - frac2 -+ bignum/bigdenom
    let base = Complex::new(q1b - frac1.re - frac2.re, -frac1.im - frac2.im);
    let minus = csqrt(Complex::new(base.re - ratio.re, base.im - ratio.im));
    let plus = csqrt(Complex::new(base.re + ratio.re, base.im + ratio.im));
    let half_s = cmul_real(0.5, sqrt_s);
    let half_minus = cmul_real(0.5, minus);
    let half_plus = cmul_real(0.5, plus);
    // q0 - 0.5*csqrt(s) -+ 0.5*csqrt(...)
    let left = Complex::new(q0 - half_s.re, -half_s.im);
    let right = Complex::new(q0 + half_s.re, half_s.im);
    let f1 = Complex::new(left.re - half_minus.re, left.im - half_minus.im);
    let f2 = Complex::new(left.re + half_minus.re, left.im + half_minus.im);
    let f3 = Complex::new(right.re - half_plus.re, right.im - half_plus.im);
    let f4 = Complex::new(right.re + half_plus.re, right.im + half_plus.im);
    [f1, f2, f3, f4]
        .into_iter()
        .any(|f| real_root_inside(f, fstart, fend))
}

/// `CrossZeroP5`: whether the quintic polynomial crosses zero in
/// `[fstart, fend]`, from the roots of `gsl_poly_complex_solve` (packed as
/// `z[2k]` real, `z[2k + 1]` imaginary part).
///
/// LAL tests the third to fifth roots with shifted indices (`z[4]` as the
/// imaginary part of the third root, `z[5]` and `z[7]` for the fourth and
/// fifth) and so can miss a real root or take a complex one for a crossing;
/// every root is tested with its own parts here (docs/lal-quirks.md).
fn cross_zero_p5(a: &[f64; 6], fstart: f64, fend: f64) -> bool {
    let z = match poly_roots::complex_solve(a) {
        Some(z) => z,
        // gsl_poly_complex_solve failed: GSL's error handler aborts.
        None => panic!("IMRPhenomXHM: gsl_poly_complex_solve failed in CrossZeroP5"),
    };
    let test = |im: f64, re: f64| im.abs() < ROOT_THRESHOLD && inside_interval(re, fstart, fend);
    (0..5).any(|k| test(z[2 * k + 1], z[2 * k]))
}

/// `ChoosePolOrder`: lower the order of the intermediate polynomial until
/// it does not cross zero in `[fstart, fend]` (its inverse is the
/// amplitude). Updates the reconstruction version and the coefficients and
/// returns `InterAmpPolOrder`.
fn choose_pol_order(
    version: &mut i32,
    delta: &mut [f64; 6],
    inputs: &PolyInputs,
    fstart: f64,
    fend: f64,
) -> i32 {
    loop {
        match *version {
            105 => {
                if cross_zero_p5(delta, fstart, fend) {
                    *delta = poly::coefficients(inputs, 1042);
                    *version = 104;
                } else {
                    return 105;
                }
            }
            1042 | 104 => {
                if cross_zero_p4(delta, fstart, fend) {
                    *delta = poly::coefficients(inputs, 1032);
                    *version = 1032;
                } else {
                    return 104;
                }
            }
            1032 => {
                if cross_zero_p3(delta, fstart, fend) {
                    *delta = poly::coefficients(inputs, 102);
                    *version = 102;
                } else {
                    return 103;
                }
            }
            102 => {
                if cross_zero_p2(delta, fstart, fend) {
                    *delta = poly::coefficients(inputs, 101);
                    *version = 101;
                } else {
                    return 102;
                }
            }
            // The linear reconstruction connects two positive values.
            101 => return 101,
            other => panic!("IMRPhenomXHM: intermediate amplitude version {other} is not valid"),
        }
    }
}

/// `IMRPhenomXHM_Intermediate_Amp_Ansatz` (release 122019): the inverse of
/// the polynomial of order `order` (1042 for the first of two intermediate
/// regions), in strain units without `amp0`.
#[inline]
fn intermediate_ansatz(powers: &UsefulPowers, amp: &HMAmpCoefficients, order: i32) -> f64 {
    let c = &amp.v122019;
    let a = &c.delta;
    let ff = powers.itself;
    let polynomial = match order {
        101 => fma(a[1], ff, a[0]),
        102 => fma(a[2], powers.two, fma(a[1], ff, a[0])),
        103 => fma(
            a[3],
            powers.three,
            fma(a[2], powers.two, fma(a[1], ff, a[0])),
        ),
        104 => fma(
            a[4],
            powers.four,
            fma(
                a[3],
                powers.three,
                fma(a[2], powers.two, fma(a[1], ff, a[0])),
            ),
        ),
        1042 => {
            let a = &c.alpha;
            fma(
                a[4],
                powers.four,
                fma(
                    a[3],
                    powers.three,
                    fma(a[2], powers.two, fma(a[1], ff, a[0])),
                ),
            )
        }
        105 => fma(
            a[5],
            powers.five,
            fma(
                a[4],
                powers.four,
                fma(
                    a[3],
                    powers.three,
                    fma(a[2], powers.two, fma(a[1], ff, a[0])),
                ),
            ),
        ),
        _ => panic!("IMRPhenomXHM: InterAmpPolOrder {order} is not valid"),
    };
    amp.amp_norm / polynomial
}

/* ---------------------------------------------------------------------- */
/*                      Coefficients and amplitude                         */
/* ---------------------------------------------------------------------- */

/// `IMRPhenomXHM_GetAmplitudeCoefficients` for the 122019 release.
/// `rotate(powers, amp)` is `SpheroidalToSpherical` for the (3,2) with the
/// coefficients computed so far.
pub(crate) fn get_amplitude_coefficients<R>(
    amp: &mut HMAmpCoefficients,
    wfhm: &mut HMWaveformStruct,
    wf: &WaveformStruct,
    rotate: R,
) where
    R: Fn(&UsefulPowers, &HMAmpCoefficients) -> Complex<f64>,
{
    amp.amp_norm = wf.amp_norm;
    amp.pn_dominant = wf.amp_norm * (2. / wfhm.emm as f64).powf(-7. / 6.);
    wfhm.f_amp_rd_falloff = 0.;
    amp.n_coefficients_inter = 0;
    amp.n_coefficients_rd_aux = 0;
    let mode_tag = wfhm.mode_tag;
    let mode_int = wfhm.mode_int;
    let mut flags = wfhm.v122019;
    let pnr = flags.pnr;

    // Options for the extrapolation outside the calibration region.
    amp.v122019.use_insp_ansatz_ringdown =
        (mode_tag == 44 || mode_tag == 33) && wf.q > 7. && wf.chi1l > 0.95;
    let wavy_insp = mode_tag == 21 || mode_tag == 32;
    let wavy_int = mode_tag == 21;

    amp.f_amp_match_in = amplitude_fcut_insp(wfhm, wf);
    amp.f_amp_match_im = amplitude_fcut_rd(wfhm, wf);

    // Frequencies of the intermediate collocation points, computed before
    // the inspiral veto changes fAmpMatchIN.
    let df = amp.f_amp_match_im - amp.f_amp_match_in;
    let freqs_inter = [
        amp.f_amp_match_in + df / 3.,
        amp.f_amp_match_in + df * 2. / 3.,
    ];

    /* Inspiral */
    crate::imrphenomx::hm_inspiral::get_pn_amplitude_coefficients(amp, wfhm, wf);
    if flags.use_famp_pn {
        amp.v122019.pn21 = pn21_coefficients(wf);
        amp.pn = [Complex::new(0., 0.); 7];
    }

    let fcut = amp.f_amp_match_in;
    amp.colloc_freqs_insp = [fcut, 0.75 * fcut, 0.5 * fcut];
    let mut powers_of_fcut_insp = UsefulPowers::new(fcut);
    let mut powers_of_f1 = UsefulPowers::new(amp.colloc_freqs_insp[0]);
    let mut powers_of_f2 = UsefulPowers::new(amp.colloc_freqs_insp[1]);
    let powers_of_f3 = UsefulPowers::new(amp.colloc_freqs_insp[2]);
    // The PN amplitude is rescaled by the leading order of the 22 here.
    for (i, powers) in [&powers_of_f1, &powers_of_f2, &powers_of_f3]
        .into_iter()
        .enumerate()
    {
        amp.pn_amplitude_insp[i] = pn_amp_rescaled(powers, wfhm, amp);
    }
    let insp_fits = fits::inspiral_fits(mode_int);
    for (value, fit) in amp.colloc_values_insp.iter_mut().zip(insp_fits) {
        *value = fit(wf).abs();
    }
    let colloc = amp.colloc_values_insp;
    let mut iv1 = colloc[0] - amp.pn_amplitude_insp[0];
    let mut iv2 = colloc[1] - amp.pn_amplitude_insp[1];
    let mut iv3 = colloc[2] - amp.pn_amplitude_insp[2];

    // Vetoes: outside the calibration region the collocation points can be
    // wavy, remove some of them.
    let mut version = flags.inspiral_amp_version;
    if flags.inspiral_amp_veto {
        let threshold = 0.2 / wfhm.amp_norm;
        if colloc[0] < threshold * powers_of_f1.seven_sixths {
            iv1 = 0.;
            version = 2;
        }
        if colloc[1] < threshold * powers_of_f2.seven_sixths {
            iv2 = 0.;
            version -= 1;
        }
        if colloc[2] < threshold * powers_of_f3.seven_sixths {
            iv3 = 0.;
            version -= 1;
        }
    }
    let (q, chi1l, chi2l) = (wf.q, wf.chi1l, wf.chi2l);
    // The (3,2) uses only the PN amplitude in this corner.
    if mode_tag == 32 && q > 2.5 && chi1l < -0.9 && chi2l < -0.9 {
        version = 0;
    }
    // The (3,2) and (3,3) remove the last collocation point in these corners.
    if mode_tag == 32 && q > 2.5 && chi1l < -0.6 && chi2l > 0. && iv1 != 0. {
        version -= 1;
        iv1 = 0.;
    }
    if mode_tag == 33 && 1.2 > q && q > 1. && chi1l < -0.1 && chi2l > 0. && iv1 != 0. {
        version -= 1;
        iv1 = 0.;
    }
    // Check for wavy collocation points (f3 < f2 < f1).
    if version == 3
        && wavy_insp
        && wavy_points(
            colloc[0] * powers_of_f1.m_seven_sixths,
            colloc[1] * powers_of_f2.m_seven_sixths,
            colloc[2] * powers_of_f3.m_seven_sixths,
        )
    {
        iv2 = 0.;
        version -= 1;
    }
    // Rename the remaining collocation points and frequencies.
    if iv2 == 0. {
        iv2 = iv3;
        iv3 = 0.;
        powers_of_f2 = powers_of_f3;
    }
    if iv1 == 0. {
        iv1 = iv2;
        powers_of_f1 = powers_of_f2;
        powers_of_fcut_insp = powers_of_f1;
        iv2 = iv3;
        iv3 = 0.;
        powers_of_f2 = powers_of_f3;
    }
    if version == 0 {
        // PN only: the inspiral ends at the MECO frequency of the mode.
        powers_of_fcut_insp = UsefulPowers::new(wfhm.f_meco_lm);
    }
    amp.f_amp_match_in = powers_of_fcut_insp.itself;
    flags.inspiral_amp_version = version;
    amp.v122019.rho = rho_coefficients(
        version,
        [iv1, iv2, iv3],
        &powers_of_fcut_insp,
        [&powers_of_f1, &powers_of_f2, &powers_of_f3],
    );
    amp.fcut_insp_seven_thirds = powers_of_fcut_insp.seven_thirds;
    amp.fcut_insp_eight_thirds = powers_of_fcut_insp.eight_thirds;
    amp.fcut_insp_three = powers_of_fcut_insp.three;

    /* Ringdown */
    let [alambda, lambda, sigma] = fits::ringdown_fits(mode_int);
    let mut rd = [alambda(wf).abs(), lambda(wf), sigma(wf), 1. / 12.];
    if mode_tag == 33 {
        rd[0] = fma(pnr.dev_parameter, pnr.mu[2], rd[0]);
        rd[2] = fma(pnr.dev_parameter, pnr.mu[3], rd[2]);
    }
    amp.v122019.rd = rd;
    // With extreme spins the transition to the ringdown is very sharp: the
    // ringdown at fAmpMatchIM is 0.9 times the inspiral at its end.
    if amp.v122019.use_insp_ansatz_ringdown {
        let powers_of_frd = UsefulPowers::new(amp.f_amp_match_im);
        let insp = inspiral_ansatz_rescaled(&powers_of_fcut_insp, wfhm, amp);
        let ringdown = rd_amp_ansatz(&powers_of_frd, wfhm, amp);
        amp.v122019.rd[0] = 0.9 * (insp / (ringdown / amp.v122019.rd[0])).abs();
    }

    /* Intermediate */
    let [int1, int2, int0, dint0] = fits::intermediate_fits(mode_int);
    amp.colloc_freqs_inter[..2].copy_from_slice(&freqs_inter);
    amp.colloc_values_inter[0] = int1(wf).abs();
    amp.colloc_values_inter[1] = int2(wf).abs();

    let mut f1 = powers_of_fcut_insp.itself;
    let f2 = freqs_inter[0];
    let mut f3 = freqs_inter[1];
    let f4 = amp.f_amp_match_im;
    let powers_of_ff1 = UsefulPowers::new(f1);
    let powers_of_ff4 = UsefulPowers::new(f4);

    // Values and derivatives at the boundaries (rescaled ansatzes).
    let insp_f1 = inspiral_ansatz_rescaled(&powers_of_ff1, wfhm, amp);
    let rd_f4 = if wfhm.mixing_on {
        let h = rotate(&powers_of_ff4, amp);
        h.re.hypot(h.im)
    } else {
        rd_amp_ansatz(&powers_of_ff4, wfhm, amp)
    };
    let mut d1 = nd_derivative(f1, |x| {
        inspiral_ansatz_rescaled(&UsefulPowers::new(x), wfhm, amp)
    });
    let mut d4 = if wfhm.mixing_on {
        nd_derivative(f4, |x| {
            let h = rotate(&UsefulPowers::new(x), amp);
            h.re.hypot(h.im)
        })
    } else {
        rd_amp_dansatz(&powers_of_ff4, wfhm, amp)
    };
    // Derivatives of the inverse of the full amplitudes.
    d1 = 7.0 / 6.0 * f1.powf(1.0 / 6.0) / insp_f1 - f1.powf(7.0 / 6.0) * d1 / (insp_f1 * insp_f1);
    d4 = 7.0 / 6.0 * f4.powf(1.0 / 6.0) / rd_f4 - f4.powf(7.0 / 6.0) * d4 / (rd_f4 * rd_f4);

    // By default a fifth order polynomial.
    let mut version = 105;
    let mut v1 = powers_of_ff1.m_seven_sixths * insp_f1;
    let mut v2 = amp.colloc_values_inter[0];
    let mut v3 = amp.colloc_values_inter[1];
    let v4 = powers_of_ff4.m_seven_sixths * rd_f4;
    if amp.v122019.use_insp_ansatz_ringdown {
        // No intermediate collocation points for these extreme cases.
        v2 = 1.;
        v3 = 1.;
        version = 101;
    }
    // The reconstruction is for the inverse of the amplitude.
    v1 = 1.0 / v1;
    v2 = 1.0 / v2;
    v3 = 1.0 / v3;
    let v4 = 1.0 / v4;
    if mode_tag == 33 {
        v2 = fma(pnr.dev_parameter, pnr.mu[0], v2);
        v3 = fma(pnr.dev_parameter, pnr.mu[1], v3);
    }

    if flags.amp_emr {
        // Two intermediate regions: a fourth order polynomial from fcutInsp
        // to the extra collocation point F0 with the value and derivative
        // from fits, and the usual reconstruction from F0 to fcutRD.
        let f0 = f1 + (f2 - f1) / 3.;
        amp.v122019.f_amp_match_int12 = f0;
        let mut v0 = int0(wf);
        let mut d0 = dint0(wf);
        let f0_seven_sixths = f0.powf(7.0 / 6.0);
        d0 = 7.0 / 6.0 / (v0 * f0) - d0 / (v0 * v0 * f0_seven_sixths);
        v0 = 1. / v0;
        let first = PolyInputs {
            d1,
            d4: d0,
            v1,
            v2,
            v3,
            v4: v0,
            f1,
            f2,
            f3,
            f4: f0,
        };
        let alpha = poly::coefficients(&first, 104);
        amp.v122019.alpha = [alpha[0], alpha[1], alpha[2], alpha[3], alpha[4]];
        f1 = f0;
        v1 = v0;
        d1 = d0;
    } else {
        amp.v122019.f_amp_match_int12 = 0.;
        amp.v122019.alpha = [1.; 5];
        if flags.intermediate_amp_veto && version == 105 {
            intermediate_veto(&mut v2, &mut v3, &mut version, wfhm.amp_norm);
        }
        if flags.ringdown_amp_veto {
            ringdown_veto(&mut v2, &mut v3, v4, &mut version, wfhm);
        }
        if version == 105 && wavy_int && wavy_points(v2, v3, v4) {
            v3 = v2;
            f3 = f2;
            v2 = 1.;
            version = 1042;
        }
        if (q > 40. && chi1l > 0.9 && v2 != 1. && v3 != 1.)
            || (mode_tag == 32
                && version != 101
                && ((q > 2.5 && chi1l < -0.6 && chi2l > 0.) || (chi1l < -0.9 && chi2l < -0.9)))
            || (mode_tag == 21 && wf.eta < 0.23 && chi1l > 0.7 && chi2l < -0.5)
        {
            v2 = 1.;
            v3 = 1.;
            version = 1032;
        }
    }
    // The reconstruction expects the point with a value at F3.
    if v3 == 1. {
        v3 = v2;
        f3 = f2;
        v2 = 1.;
    }
    let inputs = PolyInputs {
        d1,
        d4,
        v1,
        v2,
        v3,
        v4,
        f1,
        f2,
        f3,
        f4,
    };
    amp.v122019.poly_inputs = inputs;
    let mut delta = poly::coefficients(&inputs, version);
    // The amplitude is the inverse of the polynomial: lower its order until
    // it does not cross zero.
    amp.v122019.inter_amp_pol_order = choose_pol_order(
        &mut version,
        &mut delta,
        &inputs,
        amp.f_amp_match_in,
        amp.f_amp_match_im,
    );
    amp.v122019.delta = delta;
    flags.intermediate_amp_version = version;
    wfhm.v122019 = flags;
}

/// `IMRPhenomXHM_Amplitude_noModeMixing` / `_ModeMixing` /
/// `_ModeMixingRecycle` for the 122019 release: the amplitude (strain,
/// without `amp0`) at the frequency of `powers`, with `rotation()` the
/// spherical (3,2) ringdown (`SpheroidalToSpherical`, rescaled by the
/// leading order of the 22). Unlike the 122022 release negative amplitudes
/// are kept.
#[inline]
pub(crate) fn amplitude<F: FnOnce() -> Complex<f64>>(
    powers: &UsefulPowers,
    wfhm: &HMWaveformStruct,
    amp: &HMAmpCoefficients,
    rotation: F,
) -> f64 {
    let mf = powers.itself;
    if mf < amp.f_amp_match_in {
        inspiral_ansatz(powers, wfhm, amp)
    } else if mf >= amp.f_amp_match_im {
        let rd = if wfhm.mixing_on {
            let h = rotation();
            h.re.hypot(h.im)
        } else {
            rd_amp_ansatz(powers, wfhm, amp)
        };
        rd * (powers.m_seven_sixths * wfhm.amp_norm)
    } else if wfhm.v122019.amp_emr && mf < amp.v122019.f_amp_match_int12 {
        intermediate_ansatz(powers, amp, 1042)
    } else {
        intermediate_ansatz(powers, amp, amp.v122019.inter_amp_pol_order)
    }
}

/* ---------------------------------------------------------------------- */
/*          Generated from the LAL sources (clang's contractions)          */
/* ---------------------------------------------------------------------- */

fn rd_amp_ansatz0(ff: f64, frd: f64, fda: f64, a0: f64, lambda: f64, sigma: f64, lc: f64) -> f64 {
    let dfr = ff - frd;
    let dfd = fda * sigma;
    fda * a0.abs() * sigma * (-dfr * lambda / dfd).exp() / fma(dfr, dfr, dfd * dfd) * ff.powf(-lc)
}

fn rd_amp_dansatz0(ff: f64, frd: f64, fda: f64, a0: f64, lambda: f64, sigma: f64, lc: f64) -> f64 {
    let dfd = fda * sigma;
    let numerator = a0.abs()
        * ff.powf(-1.0 - lc)
        * fma(
            ff * ff * ff,
            lambda,
            fma(
                ff * ff,
                fma(-2.0 * lambda, frd, dfd * (2.0 + lc)),
                fma(
                    dfd * lc,
                    fma(frd, frd, dfd * dfd),
                    ff * fma(
                        dfd * dfd,
                        lambda,
                        fma(frd * frd, lambda, -(2.0 * dfd * frd * (1.0 + lc))),
                    ),
                ),
            ),
        );
    let denom = fma(ff, ff, fma(-(ff * 2.0), frd, fma(frd, frd, dfd * dfd)));
    let expon = (-(ff - frd) * lambda / dfd).exp();
    -numerator * expon / (denom * denom)
}

fn xdot_t4(mf: &UsefulPowers, c: &Pn21Coefficients) -> f64 {
    let mut acc0 = mf.four * c.xdot6;
    acc0 = fma(mf.five_thirds * mf.five_thirds, c.xdot5, acc0);
    acc0 = fma(mf.eight_thirds * mf.five_thirds, c.xdot65, acc0);
    acc0 = fma(mf.seven_thirds * mf.seven_thirds, c.xdot7, acc0);
    acc0 = fma(mf.eight_thirds * mf.seven_thirds, c.xdot75, acc0);
    acc0 = fma(mf.eight_thirds * mf.eight_thirds, c.xdot8, acc0);
    acc0 = fma(
        (mf.log * 2.0 / 3.0 + c.log2pi_two_thirds) * mf.eight_thirds * mf.eight_thirds,
        c.xdot8_log,
        acc0,
    );
    acc0 = fma(
        mf.eight_thirds * mf.eight_thirds * mf.one_third,
        c.xdot85,
        acc0,
    );
    acc0
}

fn p2_discriminant(a0: f64, a1: f64, a2: f64) -> f64 {
    fma(-4.0 * a0, a2, a1 * a1)
}

fn p3_q2(a1: f64, a2: f64, a3: f64) -> f64 {
    fma(-a2, a2, 3.0 * a1 * a3)
}

fn p3_q3(a0: f64, a1: f64, a2: f64, a3: f64) -> f64 {
    fma(
        -(27.0 * a0 * a3),
        a3,
        fma(-2.0 * a2 * a2, a2, 9.0 * a1 * a2 * a3),
    )
}

fn p3_discri(q2: f64, q3: f64) -> f64 {
    fma(4.0 * q2 * q2, q2, q3 * q3)
}

fn p4_q1(a2: f64, a3: f64, a4: f64) -> f64 {
    a3 * a3 / (4.0 * a4 * a4) - 2.0 * a2 / (3.0 * a4)
}

fn p4_q2(a0: f64, a1: f64, a2: f64, a3: f64, a4: f64) -> f64 {
    fma(12.0 * a0, a4, fma(a2, a2, -(3.0 * a1 * a3)))
}

fn p4_q3(a0: f64, a1: f64, a2: f64, a3: f64, a4: f64) -> f64 {
    fma(
        -(72.0 * a0 * a2),
        a4,
        fma(
            27.0 * a1 * a1,
            a4,
            fma(
                27.0 * a0 * a3,
                a3,
                fma(2.0 * a2 * a2, a2, -(9.0 * a1 * a2 * a3)),
            ),
        ),
    )
}

fn p4_sqrt_arg(q2: f64, q3: f64) -> f64 {
    fma(-4.0 * q2 * q2, q2, q3 * q3)
}

fn p4_bignum(a1: f64, a2: f64, a3: f64, a4: f64) -> f64 {
    -a3 * a3 * a3 / (a4 * a4 * a4) + 4.0 * a2 * a3 / (a4 * a4) - 8.0 * a1 / a4
}

fn rho1_3(
    v1: f64,
    v2: f64,
    v3: f64,
    fcut_insp: &UsefulPowers,
    f1: &UsefulPowers,
    f2: &UsefulPowers,
    f3: &UsefulPowers,
) -> f64 {
    fcut_insp.seven_thirds
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
        / (f1.seven_thirds
            * (f1.one_third - f2.one_third)
            * f2.seven_thirds
            * (f1.one_third - f3.one_third)
            * (f2.one_third - f3.one_third)
            * f3.seven_thirds)
}

fn rho2_3(
    v1: f64,
    v2: f64,
    v3: f64,
    fcut_insp: &UsefulPowers,
    f1: &UsefulPowers,
    f2: &UsefulPowers,
    f3: &UsefulPowers,
) -> f64 {
    -(fcut_insp.eight_thirds
        * fma(
            f2.seven_thirds,
            fma(-f3.three, v1, f1.three * v3),
            fma(
                f2.three,
                fma(f3.seven_thirds, v1, -(f1.seven_thirds * v3)),
                fma(
                    -(f1.three * f3.seven_thirds),
                    v2,
                    f1.seven_thirds * f3.three * v2,
                ),
            ),
        ))
        / (f1.seven_thirds
            * (f1.one_third - f2.one_third)
            * f2.seven_thirds
            * (f1.one_third - f3.one_third)
            * (f2.one_third - f3.one_third)
            * f3.seven_thirds)
}

fn rho3_3(
    v1: f64,
    v2: f64,
    v3: f64,
    fcut_insp: &UsefulPowers,
    f1: &UsefulPowers,
    f2: &UsefulPowers,
    f3: &UsefulPowers,
) -> f64 {
    fcut_insp.three
        * fma(
            f2.eight_thirds,
            fma(f3.seven_thirds, v1, -(f1.seven_thirds * v3)),
            fma(
                f1.seven_thirds * (-f1.one_third + f3.one_third) * f3.seven_thirds,
                v2,
                f2.seven_thirds * fma(-f3.eight_thirds, v1, f1.eight_thirds * v3),
            ),
        )
        / (f1.seven_thirds
            * (f1.one_third - f2.one_third)
            * f2.seven_thirds
            * (f1.one_third - f3.one_third)
            * (f2.one_third - f3.one_third)
            * f3.seven_thirds)
}

fn rho1_2(v1: f64, v2: f64, fcut_insp: &UsefulPowers, f1: &UsefulPowers, f2: &UsefulPowers) -> f64 {
    fcut_insp.seven_thirds * fma(-f2.eight_thirds, v1, f1.eight_thirds * v2)
        / (f1.seven_thirds * (f1.one_third - f2.one_third) * f2.seven_thirds)
}

fn rho2_2(v1: f64, v2: f64, fcut_insp: &UsefulPowers, f1: &UsefulPowers, f2: &UsefulPowers) -> f64 {
    -(fcut_insp.eight_thirds * fma(-f2.seven_thirds, v1, f1.seven_thirds * v2))
        / (f1.seven_thirds * (f1.one_third - f2.one_third) * f2.seven_thirds)
}

fn rho1_1(v1: f64, fcut_insp: &UsefulPowers, f1: &UsefulPowers) -> f64 {
    fcut_insp.seven_thirds * v1 / f1.seven_thirds
}

fn xdot5(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    let m12 = m1 * m1;
    let m22 = m2 * m2;
    -(m1 * m2
        * fma(
            -419126400.0,
            m22,
            fma(-838252800.0 * m1, m2, -(419126400.0 * m12)),
        )
        / 3.274425e7)
        * p2.five_thirds
        * p2.five_thirds
        * lalpi.five_thirds
        * lalpi.five_thirds
}

fn xdot6(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    let m12 = m1 * m1;
    let m22 = m2 * m2;
    let m13 = m12 * m1;
    -(m1 * m2
        * fma(
            1247400.0 * m12,
            fma(1848.0, m22, 743.0),
            fma(
                2494800.0 * m1 * m2,
                fma(462.0, m22, 743.0),
                fma(1152597600.0 * m2, m13, 926818200.0 * m22),
            ),
        )
        / 3.274425e7)
        * p2.four
        * lalpi.four
}

fn xdot65(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    let m12 = m1 * m1;
    let m22 = m2 * m2;
    let m13 = m12 * m1;
    let m14 = m13 * m1;
    let chi1 = wf.chi1l;
    let chi2 = wf.chi2l;
    -(m1 * m2
        * fma(
            -(34927200.0 * fma(chi2 * fma(-188.0, m2, 75.0), m2, 48.0 * lalpi.itself)),
            m22,
            fma(
                6566313600.0 * chi1,
                m14,
                fma(
                    13132627200.0 * chi1 * m2,
                    m13,
                    fma(
                        -(2619540000.0 * chi1),
                        m13,
                        fma(
                            -34927200.0 * m1 * m2,
                            fma(-m2, fma(75.0, chi1, 376.0 * chi2 * m2), 96.0 * lalpi.itself),
                            -(34927200.0
                                * fma(
                                    -m2,
                                    fma(75.0, chi2, 188.0 * (chi1 + chi2) * m2),
                                    48.0 * lalpi.itself,
                                )
                                * m12),
                        ),
                    ),
                ),
            ),
        )
        / 3.274425e7)
        * p2.eight_thirds
        * p2.five_thirds
        * lalpi.eight_thirds
        * lalpi.five_thirds
}

fn xdot7(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    let m12 = m1 * m1;
    let m22 = m2 * m2;
    let m13 = m12 * m1;
    let m23 = m22 * m2;
    let m14 = m13 * m1;
    let m24 = m23 * m2;
    let chi1 = wf.chi1l;
    let chi2 = wf.chi2l;
    let chi12 = chi1 * chi1;
    let chi22 = chi2 * chi2;
    -(m1 * m2
        * fma(
            -(23100.0 * m12),
            fma(
                59472.0,
                m24,
                fma(
                    -(18.0
                        * fma(
                            10206.0,
                            chi22,
                            fma(10206.0, chi12, fma(-(19908.0 * chi1), chi2, -13661.0)),
                        )),
                    m22,
                    fma(91854.0, chi12, 34103.0),
                ),
            ),
            fma(
                -(2747606400.0 * m13),
                m23,
                fma(
                    23100.0 * m1 * m2,
                    fma(
                        -2.0,
                        fma(45927.0, chi22, fma(45927.0, chi12, 34103.0)),
                        9.0 * fma(
                            10206.0,
                            chi22,
                            fma(10206.0, chi12, fma(-(19908.0 * chi1), chi2, -13661.0)),
                        ) * m22,
                    ),
                    fma(
                        -(1373803200.0 * m14),
                        m22,
                        fma(
                            207900.0
                                * m2
                                * fma(
                                    10206.0,
                                    chi22,
                                    fma(10206.0, chi12, fma(-(19908.0 * chi1), chi2, -13661.0)),
                                ),
                            m13,
                            -(23100.0 * fma(91854.0, chi22, 34103.0) * m22),
                        ),
                    ),
                ),
            ),
        )
        / 3.274425e7)
        * p2.seven_thirds
        * p2.seven_thirds
        * lalpi.seven_thirds
        * lalpi.seven_thirds
}

fn xdot75(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    let m12 = m1 * m1;
    let m22 = m2 * m2;
    let m13 = m12 * m1;
    let m23 = m22 * m2;
    let m14 = m13 * m1;
    let m24 = m23 * m2;
    let m15 = m14 * m1;
    let chi1 = wf.chi1l;
    let chi2 = wf.chi2l;
    let mut acc0 = 5821200.0 * m2 * fma(5861.0, chi1, 1701.0 * lalpi.itself) * m13;
    acc0 = fma(-4036586400.0 * chi1, m13, acc0);
    acc0 = fma(17059026600.0 * chi1, m14, acc0);
    acc0 = fma(14721814800.0 * chi1 * m2, m14, acc0);
    acc0 = fma(-(34962127200.0 * chi1 * m2), m15, acc0);
    acc0 = fma(
        207900.0
            * fma(
                2.0 * chi2 * m2,
                fma(41027.0, m2, -9708.0),
                12477.0 * lalpi.itself,
            ),
        m22,
        acc0,
    );
    acc0 = fma(-(14721814800.0 * chi2 * m13), m22, acc0);
    acc0 = fma(-(69924254400.0 * chi1 * m14), m22, acc0);
    acc0 = fma(-(34962127200.0 * (chi1 + chi2) * m13), m23, acc0);
    acc0 = fma(
        207900.0 * m12,
        fma(
            3.0 * lalpi.itself,
            fma(31752.0, m22, 4159.0),
            2.0 * m2
                * fma(
                    -(168168.0 * chi2),
                    m23,
                    fma(
                        -(35406.0 * chi1),
                        m22,
                        fma(9708.0, chi2, 41027.0 * (chi1 + chi2) * m2),
                    ),
                ),
        ),
        acc0,
    );
    acc0 = fma(
        415800.0 * m1 * m2,
        fma(
            -(84084.0 * chi2),
            m24,
            fma(
                35406.0 * chi2,
                m23,
                fma(
                    3.0 * lalpi.itself,
                    fma(7938.0, m22, 4159.0),
                    fma(9708.0 * chi1, m2, 82054.0 * chi2 * m22),
                ),
            ),
        ),
        acc0,
    );
    -(m1 * m2 * acc0 / 3.274425e7)
        * p2.eight_thirds
        * p2.seven_thirds
        * lalpi.eight_thirds
        * lalpi.seven_thirds
}

fn xdot8(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    let m12 = m1 * m1;
    let m22 = m2 * m2;
    let m13 = m12 * m1;
    let m23 = m22 * m2;
    let m14 = m13 * m1;
    let m24 = m23 * m2;
    let m15 = m14 * m1;
    let m25 = m24 * m2;
    let m16 = m15 * m1;
    let chi1 = wf.chi1l;
    let chi2 = wf.chi2l;
    let chi12 = chi1 * chi1;
    let chi22 = chi2 * chi2;
    let log4 = 1.3862943611198906;
    let euler_gamma = 0.5772156649015329;
    let mut acc0 = -(63392868000.0 * chi1 * chi2 * m2 * m14);
    acc0 = fma(-10548014400.0 * chi1 * lalpi.itself, m13, acc0);
    acc0 = fma(
        34927200.0 * chi1 * fma(-375.0, chi1, 752.0 * lalpi.itself),
        m14,
        acc0,
    );
    acc0 = fma(63392868000.0 * chi12, m15, acc0);
    acc0 = fma(-(153213984000.0 * m2 * chi12), m15, acc0);
    acc0 = fma(-(76606992000.0 * chi12), m16, acc0);
    acc0 = fma(-(63392868000.0 * chi1 * (chi1 - chi2) * m13), m22, acc0);
    acc0 = fma(
        -(51975.0
            * fma(
                228508.0,
                chi22,
                fma(1702428.0, chi12, fma(2711352.0 * chi1, chi2, 4869.0)),
            )
            * m14),
        m22,
        acc0,
    );
    acc0 = fma(
        -(103950.0
            * fma(
                228508.0,
                chi22,
                fma(228508.0, chi12, fma(2711352.0 * chi1, chi2, 4869.0)),
            )
            * m13),
        m23,
        acc0,
    );
    acc0 = fma(906328500.0 * m15, m23, acc0);
    acc0 = fma(1812657000.0 * m14, m24, acc0);
    acc0 = fma(906328500.0 * m13, m25, acc0);
    acc0 = fma(
        1925.0 * m2 * m13,
        fma(
            -2045736.0,
            lalpi.two,
            fma(
                3617892.0,
                chi22,
                fma(
                    30746952.0,
                    chi12,
                    fma(
                        27288576.0 * chi1,
                        lalpi.itself,
                        fma(13635864.0 * chi1, chi2, 56198689.0),
                    ),
                ),
            ),
        ),
        acc0,
    );
    acc0 = fma(
        -(3.0 * m22),
        fma(
            745113600.0,
            lalpi.two,
            fma(
                25535664000.0 * chi22,
                m24,
                fma(
                    -(21130956000.0 * chi22),
                    m23,
                    fma(
                        4365900000.0 * chi22,
                        m22,
                        fma(
                            2321480700.0,
                            chi22,
                            fma(
                                -(23284800.0 * chi2 * m2 * fma(376.0, m2, -151.0)),
                                lalpi.itself,
                                fma(-2277918720.0, euler_gamma, 16447322263.0),
                            ) - 2277918720.0 * log4,
                        ),
                    ),
                ),
            ),
        ),
        acc0,
    );
    acc0 = fma(
        m12,
        fma(
            -3.0,
            fma(
                745113600.0,
                lalpi.two,
                fma(2321480700.0, chi12, fma(-2277918720.0, log4, 16447322263.0)),
            ),
            fma(
                -(3850.0 * m22),
                fma(
                    2045736.0,
                    lalpi.two,
                    fma(
                        -6976422.0,
                        chi22,
                        fma(
                            -6976422.0,
                            chi12,
                            fma(
                                -(6822144.0 * (chi1 + chi2)),
                                lalpi.itself,
                                fma(13580136.0 * chi1, chi2, -56198689.0),
                            ),
                        ),
                    ),
                ),
                fma(
                    -(51975.0
                        * fma(
                            1702428.0,
                            chi22,
                            fma(228508.0, chi12, fma(2711352.0 * chi1, chi2, 4869.0)),
                        )),
                    m24,
                    fma(
                        63392868000.0 * (chi1 - chi2) * chi2,
                        m23,
                        6833756160.0 * euler_gamma + 10548014400.0 * chi2 * m2 * lalpi.itself,
                    ),
                ),
            ),
        ),
        acc0,
    );
    acc0 = fma(
        m1 * m2,
        fma(
            -6.0,
            fma(
                745113600.0,
                lalpi.two,
                fma(
                    1160740350.0,
                    chi22,
                    fma(1160740350.0, chi12, fma(-2277918720.0, log4, 16447322263.0)),
                ),
            ),
            fma(
                -(1925.0 * m22),
                fma(
                    2045736.0,
                    lalpi.two,
                    fma(
                        -30746952.0,
                        chi22,
                        fma(
                            -3617892.0,
                            chi12,
                            fma(
                                -(27288576.0 * chi2),
                                lalpi.itself,
                                fma(-(13635864.0 * chi1), chi2, -56198689.0),
                            ),
                        ),
                    ),
                ),
                fma(
                    -(153213984000.0 * chi22),
                    m24,
                    fma(
                        -(63392868000.0 * chi1 * chi2),
                        m23,
                        13667512320.0 * euler_gamma + 10548014400.0 * chi1 * m2 * lalpi.itself,
                    ),
                ),
            ),
        ),
        acc0,
    );
    -(m1 * m2 * acc0 / 3.274425e7)
        * p2.eight_thirds
        * p2.eight_thirds
        * lalpi.eight_thirds
        * lalpi.eight_thirds
}

fn xdot8_log(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    -(m1 * m2 * 3416878080.0) / 3.274425e7
        * p2.eight_thirds
        * p2.eight_thirds
        * lalpi.eight_thirds
        * lalpi.eight_thirds
}

fn xdot85(wf: &WaveformStruct, p2: &UsefulPowers, lalpi: &UsefulPowers) -> f64 {
    let m1 = wf.m1;
    let m2 = wf.m2;
    let m12 = m1 * m1;
    let m22 = m2 * m2;
    let m13 = m12 * m1;
    let m23 = m22 * m2;
    let m14 = m13 * m1;
    let m24 = m23 * m2;
    let m15 = m14 * m1;
    let m25 = m24 * m2;
    let m16 = m15 * m1;
    let chi1 = wf.chi1l;
    let chi2 = wf.chi2l;
    let mut acc0 = 1925.0
        * m2
        * fma(
            -12912300.0,
            lalpi.itself,
            fma(97151928.0, chi1, 6613488.0 * chi2),
        )
        * m13;
    acc0 = fma(-14891068500.0 * chi1, m13, acc0);
    acc0 = fma(87143248500.0 * chi1, m14, acc0);
    acc0 = fma(33313480200.0 * chi1 * m2, m14, acc0);
    acc0 = fma(-(198816225300.0 * chi1 * m2), m15, acc0);
    acc0 = fma(
        57750.0
            * fma(
                2.0 * chi2 * m2,
                fma(754487.0, m2, -128927.0),
                7947.0 * lalpi.itself,
            ),
        m22,
        acc0,
    );
    acc0 = fma(-(33313480200.0 * chi2 * m13), m22, acc0);
    acc0 = fma(
        -(138600.0
            * fma(
                182990.0,
                lalpi.itself,
                fma(3399633.0, chi1, 530712.0 * chi2),
            )
            * m14),
        m22,
        acc0,
    );
    acc0 = fma(-(35665037100.0 * chi1 * m15), m22, acc0);
    acc0 = fma(84184254000.0 * chi1 * m16, m22, acc0);
    acc0 = fma(
        -(23100.0 * m1 * m2),
        fma(
            15.0 * lalpi.itself,
            fma(71735.0, m22, -2649.0),
            m2 * fma(
                -chi1,
                fma(551124.0, m2, 644635.0),
                chi2 * m2 * fma(8606763.0, m22, fma(-1442142.0, m2, -8095994.0)),
            ),
        ),
        acc0,
    );
    acc0 = fma(
        -(69300.0 * fma(4991769.0, chi1 + chi2, 731960.0 * lalpi.itself) * m13),
        m23,
        acc0,
    );
    acc0 = fma(35665037100.0 * chi2 * m14, m23, acc0);
    acc0 = fma(170726094000.0 * chi1 * m15, m23, acc0);
    acc0 = fma(2357586000.0 * chi2 * m15, m23, acc0);
    acc0 = fma(35665037100.0 * chi1 * m13, m24, acc0);
    acc0 = fma(88899426000.0 * (chi1 + chi2) * m14, m24, acc0);
    acc0 = fma(
        9702000.0 * fma(243.0, chi1, 17597.0 * chi2) * m13,
        m25,
        acc0,
    );
    acc0 = fma(
        -(11550.0 * m12),
        fma(
            15.0 * lalpi.itself,
            fma(146392.0, m24, fma(286940.0, m22, -2649.0)),
            2.0 * m2
                * fma(
                    -(3644340.0 * chi2),
                    m25,
                    fma(
                        1543941.0 * chi2,
                        m24,
                        fma(
                            54.0 * fma(58968.0, chi1, 377737.0 * chi2),
                            m23,
                            fma(
                                1442142.0 * chi1,
                                m22,
                                fma(-644635.0, chi2, -(4874683.0 * (chi1 + chi2) * m2)),
                            ),
                        ),
                    ),
                ),
        ),
        acc0,
    );
    -(m1 * m2 * acc0 / 3.274425e7)
        * p2.eight_thirds
        * p2.eight_thirds
        * p2.one_third
        * lalpi.eight_thirds
        * lalpi.eight_thirds
        * lalpi.one_third
}

#[cfg(test)]
mod tests {
    use num_complex::Complex;

    use super::RELEASE_122019;
    use crate::imrphenomx::hm_internals::HMMode;
    use crate::imrphenomx::internals::{UsefulPowers, WaveformStruct};
    use crate::imrphenomx::phase::PhaseCoefficients;

    /// The decomposition of the (3,2) ringdown used by the grid paths
    /// reproduces the mode with the original ringdown ansatz.
    #[test]
    fn mixing_ringdown_terms_match_the_mode() {
        for (m1, m2, chi1, chi2) in [(36.0, 29.0, 0.3, -0.4), (79.0, 0.27, -0.89, 0.96)] {
            let mut wf =
                WaveformStruct::new(m1, m2, chi1, chi2, 0.0, 20.0, 0.4, 10.0, 0.0, 500.0, 0.0)
                    .unwrap();
            wf.hm_release = RELEASE_122019;
            let mode = HMMode::new(3, 2, &wf, &PhaseCoefficients::new(&wf));
            let start = mode.mixing_ringdown_start().unwrap();
            let [c1, c2] = mode.mixing_factors();
            for i in 0..50 {
                let mf = start * (1.0 + 0.04 * i as f64);
                let powers = UsefulPowers::new(mf);
                let (amp, phi) = mode.amplitude_phase(&powers, None);
                let h = Complex::from_polar(amp, phi);
                let (lin, [(a1, p1), (a2, p2)]) = mode.mixing_ringdown_terms(&powers);
                let g = Complex::from_polar(1.0, lin)
                    * (c1 * Complex::from_polar(a1, p1) + c2 * Complex::from_polar(a2, p2));
                assert!((g - h).norm() <= 1e-12 * h.norm(), "{mf}: {g} != {h}");
            }
        }
    }

    /// The root tests of `ChoosePolOrder`; `CrossZeroP5` finds every real
    /// root, including the third that LAL's indexing misses.
    #[test]
    fn polynomial_roots_in_interval() {
        // (f - 0.05)(f - 0.2)(f + 1) = f^3 + 0.75 f^2 - 0.24 f + 0.01
        let cubic = [0.01, -0.24, 0.75, 1.0, 0.0, 0.0];
        assert!(super::cross_zero_p3(&cubic, 0.04, 0.06));
        assert!(!super::cross_zero_p3(&cubic, 0.06, 0.15));
        // times (f^2 + 1): GSL returns the roots -1, 0.2, 0.05, i, -i (LAL
        // misses the third, 0.05)
        let quintic = [0.01, -0.24, 0.76, 0.76, 0.75, 1.0];
        assert!(super::cross_zero_p5(&quintic, 0.15, 0.25));
        assert!(super::cross_zero_p5(&quintic, 0.04, 0.06));
        assert!(!super::cross_zero_p5(&quintic, 0.06, 0.15));
        // cubic times (f - 0.3)
        let quartic = [-0.003, 0.082, -0.465, 0.45, 1.0, 0.0];
        assert!(super::cross_zero_p4(&quartic, 0.28, 0.32));
        assert!(!super::cross_zero_p4(&quartic, 0.21, 0.28));
        // f^2 - 0.01 has a root at 0.1
        let quadratic = [-0.01, 0.0, 1.0, 0.0, 0.0, 0.0];
        assert!(super::cross_zero_p2(&quadratic, 0.05, 0.15));
        assert!(!super::cross_zero_p2(&quadratic, 0.15, 0.25));
    }
}
