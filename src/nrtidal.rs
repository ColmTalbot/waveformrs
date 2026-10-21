// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The NRTidal corrections for binary neutron stars, NRTidalv2 (Dietrich et
//! al. 2019, arXiv:1905.06011) and NRTidalv3 (Abac et al. 2024,
//! arXiv:2311.07456), following `LALSimNRTunedTides.c`, with the
//! quadrupole-monopole parameter of `XLALSimInspiralEOSQfromLambda`
//! (`LALSimInspiralEOS.c`) and the spin-induced octupole of
//! `XLALSimUniversalRelationSpinInducedOctupoleVSSpinInducedQuadrupole`
//! (`LALSimUniversalRelations.c`).
//!
//! The tidal deformabilities are LAL's `lambda1`, `lambda2`: the
//! dimensionless `Lambda_i = (2 / 3) k_2 / C_i^5`. The models adding these
//! corrections are `IMRPhenomD_NRTidalv2`
//! ([`crate::imrphenomd::IMRPhenomD::with_nrtidal_v2`]),
//! `IMRPhenomXAS_NRTidalv2/v3` and `IMRPhenomXP_NRTidalv2/v3`.

use std::f64::consts::PI;

use crate::constants::{MSUN_SI, MTSUN_SI};
use crate::spintaylor::fp::fma;

/// The version of the NRTidal corrections (LAL's `NRTidal_version_type`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NRTidalVersion {
    /// NRTidalv2 (`NRTidalv2_V`, LAL's `PhenomXTidalFlag` 2).
    V2,
    /// NRTidalv3 (`NRTidalv3_V`, LAL's `PhenomXTidalFlag` 3).
    V3,
}

/// `PlanckTaper`: 0 below `t1`, 1 above `t2`, smooth in between.
pub(crate) fn planck_taper(t: f64, t1: f64, t2: f64) -> f64 {
    if t <= t1 {
        0.
    } else if t >= t2 {
        1.
    } else {
        1. / (((t2 - t1) / (t - t1) + (t2 - t1) / (t - t2)).exp() + 1.)
    }
}

/// `XLALSimNRTunedTidesComputeKappa2T`: the tidal coupling constant
/// `kappa_2^T` (Eq. 2 of arXiv:1706.02969) for masses in kg, in either
/// order.
pub fn kappa2t(m1_si: f64, m2_si: f64, lambda1: f64, lambda2: f64) -> f64 {
    let (m1_si, m2_si, lambda1, lambda2) = if m1_si >= m2_si {
        (m1_si, m2_si, lambda1, lambda2)
    } else {
        (m2_si, m1_si, lambda2, lambda1)
    };
    let m1 = m1_si / MSUN_SI;
    let m2 = m2_si / MSUN_SI;
    let mtot = m1 + m2;
    let xa = m1 / mtot;
    let xb = m2 / mtot;
    // 2 k_2^A / C_A^5 = 3 lambda1
    let term1 = (1.0 + 12.0 * xb / xa) * xa.powf(5.0) * lambda1;
    let term2 = (1.0 + 12.0 * xa / xb) * xb.powf(5.0) * lambda2;
    (3.0 / 13.0) * (term1 + term2)
}

/// `XLALSimNRTunedTidesMergerFrequency`: the merger frequency (Hz) of
/// NRTidal and NRTidalv2 for the total mass in solar masses and the mass
/// ratio `q >= 1`.
pub fn merger_frequency(mtot: f64, kappa2t: f64, q: f64) -> f64 {
    let a_0 = 0.3586;
    let n_1 = 3.35411203e-2;
    let n_2 = 4.31460284e-5;
    let d_1 = 7.54224145e-2;
    let d_2 = 2.23626859e-4;
    let kappa2t2 = kappa2t * kappa2t;
    let num = 1.0 + n_1 * kappa2t + n_2 * kappa2t2;
    let den = 1.0 + d_1 * kappa2t + d_2 * kappa2t2;
    let q_0 = a_0 / q.sqrt();
    // dimensionless angular frequency of merger
    let momega_merger = q_0 * (num / den);
    momega_merger / (2.0 * PI) / (mtot * MTSUN_SI)
}

/// `XLALSimNRTunedTidesMergerFrequency_v3`: the merger frequency (Hz) of
/// NRTidalv3 (Eq. 23 of arXiv:2210.16366) for the total mass in solar
/// masses, the tidal deformabilities, the mass ratio `q >= 1` and the
/// aligned spins.
pub fn merger_frequency_v3(
    mtot: f64,
    lambda1: f64,
    lambda2: f64,
    q: f64,
    chi1: f64,
    chi2: f64,
) -> f64 {
    let xa = q / (1.0 + q);
    let xb = 1.0 - xa;
    // EnforcePrimaryMassIsm1_v3 swaps nothing for q >= 1
    let nu = xa * xb;
    let xa_2 = xa * xa;
    let xa_3 = xa_2 * xa;
    let xb_2 = xb * xb;
    let xb_3 = xb_2 * xb;
    let kappa2eff = 3.0 * nu * (xa_3 * lambda1 + xb_3 * lambda2);

    let a_0 = 0.22;
    let a_1m = 0.80;
    let a_1s = 0.25;
    // not in arXiv:2311.07456, from Table 2 of arXiv:2210.16366
    let b_1s = -1.99;
    let a_1t = 0.0485;
    let a_2t = 0.00000586;
    let a_3t = 0.10;
    let a_4t = 0.000186;
    let b_1t = 1.80;
    let b_2t = 599.99;
    let b_3t = 7.80;
    let b_4t = 84.76;

    let xval = 1.0 - 4.0 * nu;
    let p_1s = a_1s * (1.0 + b_1s * xval);
    let p_1t = a_1t * (1.0 + b_1t * xval);
    let p_2t = a_2t * (1.0 + b_2t * xval);
    let p_3t = a_3t * (1.0 + b_3t * xval);
    let p_4t = a_4t * (1.0 + b_4t * xval);
    let kappa2eff2 = kappa2eff * kappa2eff;
    let sval = xa_2 * chi1 + xb_2 * chi2;
    let qm = 1.0 + a_1m * xval;
    let qs = 1.0 + p_1s * sval;
    let num = 1.0 + p_1t * kappa2eff + p_2t * kappa2eff2;
    let den = 1.0 + p_3t * kappa2eff + p_4t * kappa2eff2;
    let qt = num / den;
    let qfit = a_0 * qm * qs * qt;
    let momega_merger = nu * qfit * (2.0 * PI);
    momega_merger / (2.0 * PI) / (mtot * MTSUN_SI)
}

/// `SimNRTunedTidesFDTidalAmplitude`: the NRTidalv2 amplitude correction
/// (Eq. 24 of arXiv:1905.06011) at `f_hz`.
pub(crate) fn tidal_amplitude(f_hz: f64, mtot: f64, kappa2t: f64) -> f64 {
    let m_sec = mtot * MTSUN_SI;
    let prefac = 9.0 * kappa2t;
    let x = (PI * m_sec * f_hz).powf(2.0 / 3.0);
    let n1 = 4.157407407407407;
    let n289 = 2519.111111111111;
    let d = 13477.8073677;
    let poly = (1.0 + n1 * x + n289 * x.powf(2.89)) / (1. + d * x.powf(4.));
    -prefac * x.powf(3.25) * poly
}

/// `XLALSimNRTunedTidesSetFDTidalPhase_v2_Coeffs`: `c_Newt, n_1, n_3/2,
/// n_2, n_5/2, n_3, d_1, d_3/2, d_2`.
pub(crate) const NRTIDALV2_COEFFS: [f64; 9] = [
    2.4375,
    -12.615214237993088,
    19.0537346970349,
    -21.166863146081035,
    90.55082156324926,
    -60.25357801943598,
    -15.111207827736678,
    22.195327350624694,
    8.064109635305156,
];

/// `SimNRTunedTidesFDTidalPhase_v2` (Eq. 22 of arXiv:1905.06011).
pub(crate) fn tidal_phase_v2(f_hz: f64, xa: f64, xb: f64, mtot: f64, kappa2t: f64) -> f64 {
    let m_omega = PI * f_hz * (mtot * MTSUN_SI);
    let pn_x = m_omega.powf(2.0 / 3.0);
    let pn_x_2 = pn_x * pn_x;
    let pn_x_3 = pn_x * pn_x_2;
    let pn_x_3over2 = pn_x.powf(3.0 / 2.0);
    let pn_x_5over2 = pn_x.powf(5.0 / 2.0);
    let [c_newt, n_1, n_3over2, n_2, n_5over2, n_3, d_1, d_3over2, d_2] = NRTIDALV2_COEFFS;
    let tidal_phase = -kappa2t * c_newt / (xa * xb) * pn_x_5over2;
    let num = 1.0
        + (n_1 * pn_x)
        + (n_3over2 * pn_x_3over2)
        + (n_2 * pn_x_2)
        + (n_5over2 * pn_x_5over2)
        + (n_3 * pn_x_3);
    let den = 1.0 + (d_1 * pn_x) + (d_3over2 * pn_x_3over2) + (d_2 * pn_x_2);
    tidal_phase * (num / den)
}

/// `XLALSimNRTunedTidesSetFDTidalPhase_PN_Coeffs`: the 7.5PN tidal phase
/// coefficients of each star (Eq. 45 of arXiv:2311.07456), `c_Newt, c_1,
/// c_3/2, c_2, c_5/2` of star A then of star B.
pub(crate) fn pn_coeffs(xa: f64) -> [f64; 10] {
    let xb = 1.0 - xa;
    let xa_2 = xa * xa;
    let xa_3 = xa_2 * xa;
    let xa_4 = xa_3 * xa;
    let xa_5 = xa_4 * xa;
    let xb_2 = xb * xb;
    let xb_3 = xb_2 * xb;
    let xb_4 = xb_3 * xb;
    let xb_5 = xb_4 * xb;
    let den_a = 11.0 * xa - 12.0;
    let den_b = 11.0 * xb - 12.0;
    [
        -3.0 * den_a / (16. * xa * xb_2),
        (-1300.0 * xa_3 + 11430.0 * xa_2 + 4595.0 * xa - 15895.0) / (672.0 * den_a),
        -PI,
        (22861440.0 * xa_5 - 102135600.0 * xa_4
            + 791891100.0 * xa_3
            + 874828080.0 * xa_2
            + 216234195.0 * xa
            - 1939869350.0)
            / (27433728.0 * den_a),
        -PI * (10520.0 * xa_3 - 7598.0 * xa_2 + 22415.0 * xa - 27719.0) / (672.0 * den_a),
        -3.0 * den_b / (16. * xb * xa_2),
        (-1300.0 * xb_3 + 11430.0 * xb_2 + 4595.0 * xb - 15895.0) / (672.0 * den_b),
        -PI,
        (22861440.0 * xb_5 - 102135600.0 * xb_4
            + 791891100.0 * xb_3
            + 874828080.0 * xb_2
            + 216234195.0 * xb
            - 1939869350.0)
            / (27433728.0 * den_b),
        -PI * (10520.0 * xb_3 - 7598.0 * xb_2 + 22415.0 * xb - 27719.0) / (672.0 * den_b),
    ]
}

/// `XLALSimNRTunedTidesSetFDTidalPhase_v3_Coeffs`: the coefficients of the
/// effective Love number enhancement (`s_1, s_2, s_3, exp(s_2 s_3)`), the
/// effective deformabilities `kappa_A, kappa_B` and the Padé coefficients
/// of NRTidalv3 (`n_5/2, n_3, d_1` and the PN-constrained `n_1, n_3/2, n_2,
/// d_3/2` of each star), in LAL's order.
pub(crate) fn v3_coeffs(
    xa: f64,
    mtot: f64,
    lambda1: f64,
    lambda2: f64,
    pn: &[f64; 10],
) -> [f64; 20] {
    let xb = 1.0 - xa;
    let q = xa / xb;
    let s10 = 1.273000423;
    let s11 = 3.64169971e-03;
    let s12 = 1.76144380e-03;
    let s20 = 2.78793291e+01;
    let s21 = 1.18175396e-02;
    let s22 = -5.39996790e-03;
    let s30 = 1.42449682e-01;
    let s31 = -1.70505852e-05;
    let s32 = 3.38040594e-05;

    let m1_si = xa * mtot * MSUN_SI;
    let m2_si = xb * mtot * MSUN_SI;
    let kappa2t = kappa2t(m1_si, m2_si, lambda1, lambda2);

    let mut c = [0.0; 20];
    c[0] = s10 + s11 * kappa2t + s12 * q * kappa2t;
    c[1] = s20 + s21 * kappa2t + s22 * q * kappa2t;
    c[2] = s30 + s31 * kappa2t + s32 * q * kappa2t;
    let s2s3 = c[1] * c[2];
    c[3] = s2s3.cosh() + s2s3.sinh();
    c[4] = 3.0 * xb * xa * xa * xa * xa * lambda1;
    c[5] = 3.0 * xa * xb * xb * xb * xb * lambda2;

    let alpha = -8.08155404e-03;
    let beta = -1.13695919e+00;
    let kappa_a_alpha = (c[4] + 1.).powf(alpha);
    let kappa_b_alpha = (c[5] + 1.).powf(alpha);
    let xa_beta = xa.powf(beta);
    let xb_beta = xb.powf(beta);

    let n_5over20 = -9.40654388e+02;
    let n_5over21 = 6.26517157e+02;
    let n_5over22 = 5.53629706e+02;
    let n_5over23 = 8.84823087e+01;
    let n_30 = 4.05483848e+02;
    let n_31 = -4.25525054e+02;
    let n_32 = -1.92004957e+02;
    let n_33 = -5.10967553e+01;
    let d_10 = 3.80343306e+00;
    let d_11 = -2.52026996e+01;
    let d_12 = -3.08054443e+00;

    c[6] = n_5over20 + n_5over21 * xa + n_5over22 * kappa_a_alpha + n_5over23 * xa_beta;
    c[7] = n_30 + n_31 * xa + n_32 * kappa_a_alpha + n_33 * xa_beta;
    c[8] = d_10 + d_11 * xa + d_12 * xa_beta;
    c[9] = n_5over20 + n_5over21 * xb + n_5over22 * kappa_b_alpha + n_5over23 * xb_beta;
    c[10] = n_30 + n_31 * xb + n_32 * kappa_b_alpha + n_33 * xb_beta;
    c[11] = d_10 + d_11 * xb + d_12 * xb_beta;

    let [_, c_1a, c_3over2a, c_2a, c_5over2a, _, c_1b, c_3over2b, c_2b, c_5over2b] = *pn;
    let inv_c1_a = 1.0 / c_1a;
    c[12] = c_1a + c[8];
    c[13] = ((c_1a * c_3over2a) - c_5over2a - (c_3over2a) * c[8] + c[6]) * inv_c1_a;
    c[14] = c_2a + c_1a * c[8];
    c[15] = -(c_5over2a + c_3over2a * c[8] - c[6]) * inv_c1_a;
    let inv_c1_b = 1.0 / c_1b;
    c[16] = c_1b + c[11];
    c[17] = ((c_1b * c_3over2b) - c_5over2b - (c_3over2b) * c[11] + c[9]) * inv_c1_b;
    c[18] = c_2b + c_1b * c[11];
    c[19] = -(c_5over2b + c_3over2b * c[11] - c[9]) * inv_c1_b;
    c
}

/// `SimNRTunedTidesFDTidalPhase_v3` (Eqs. 27, 30 of arXiv:2311.07456).
fn tidal_phase_v3(f_hz: f64, mtot: f64, c: &[f64; 20], pn: &[f64; 10]) -> f64 {
    let m_omega = PI * f_hz * (mtot * MTSUN_SI);
    let s1 = c[0];
    let s2 = c[1];
    let exps2s3 = c[3];
    let s2mom = -s2 * m_omega * 2.0;
    let exps2mom = s2mom.cosh() + s2mom.sinh();
    // the effective Love number enhancement factor, Eq. (27)
    let dynk2bar = 1.0 + (s1 - 1.) * (1.0 / (1.0 + exps2mom * exps2s3))
        - ((s1 - 1.0) / (1.0 + exps2s3))
        - 2.0 * m_omega * (s1 - 1.) * s2 * exps2s3 / ((1.0 + exps2s3) * (1.0 + exps2s3));

    let pn_x = m_omega / m_omega.cbrt();
    let pn_x_2 = pn_x * pn_x;
    let pn_x_3 = pn_x * pn_x_2;
    let pn_x_3over2 = pn_x * pn_x.sqrt();
    let pn_x_5over2 = pn_x_3over2 * pn_x;

    let dynkappa_a = c[4] * dynk2bar;
    let dynkappa_b = c[5] * dynk2bar;
    let factor_a = -pn[0] * pn_x_5over2 * dynkappa_a;
    let factor_b = -pn[5] * pn_x_5over2 * dynkappa_b;
    // Padé approximant, Eq. (32)
    let num_a = 1.0
        + (c[12] * pn_x)
        + (c[13] * pn_x_3over2)
        + (c[14] * pn_x_2)
        + (c[6] * pn_x_5over2)
        + (c[7] * pn_x_3);
    let den_a = 1.0 + (c[8] * pn_x) + (c[15] * pn_x_3over2);
    let num_b = 1.0
        + (c[16] * pn_x)
        + (c[17] * pn_x_3over2)
        + (c[18] * pn_x_2)
        + (c[9] * pn_x_5over2)
        + (c[10] * pn_x_3);
    let den_b = 1.0 + (c[11] * pn_x) + (c[19] * pn_x_3over2);
    factor_a * (num_a / den_a) + factor_b * (num_b / den_b)
}

/// `SimNRTunedTidesFDTidalPhase_PN`: the 7.5PN tidal phase (Eqs. 22, 45 of
/// arXiv:2311.07456).
fn tidal_phase_pn(
    f_hz: f64,
    xa: f64,
    mtot: f64,
    lambda1: f64,
    lambda2: f64,
    pn: &[f64; 10],
) -> f64 {
    let m_omega = PI * f_hz * (mtot * MTSUN_SI);
    let xb = 1.0 - xa;
    let pn_x = m_omega / m_omega.cbrt();
    let pn_x_2 = pn_x * pn_x;
    let pn_x_3over2 = pn_x * pn_x.sqrt();
    let pn_x_5over2 = pn_x_3over2 * pn_x;
    let kappa_a = 3.0 * xb * xa * xa * xa * xa * lambda1;
    let kappa_b = 3.0 * xa * xb * xb * xb * xb * lambda2;
    let factor_a = -pn[0] * pn_x_5over2 * kappa_a;
    let factor_b = -pn[5] * pn_x_5over2 * kappa_b;
    let phase_a = factor_a
        * (1.0 + (pn[1] * pn_x) + (pn[2] * pn_x_3over2) + (pn[3] * pn_x_2) + (pn[4] * pn_x_5over2));
    let phase_b = factor_b
        * (1.0 + (pn[6] * pn_x) + (pn[7] * pn_x_3over2) + (pn[8] * pn_x_2) + (pn[9] * pn_x_5over2));
    phase_a + phase_b
}

/// The tidal phase, amplitude correction and Planck window of
/// [`tidal_series`] on a frequency grid.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TidalSeries {
    /// The tidal phase `phi_T(f)`, subtracted from the phase of the
    /// (2,2) mode.
    pub phase: Vec<f64>,
    /// The tidal amplitude correction `A_T(f)` (Eq. 24 of
    /// arXiv:1905.06011).
    pub amplitude: Vec<f64>,
    /// `1 - PlanckTaper(f)`, which switches the waveform off after merger.
    pub taper: Vec<f64>,
}

/// `XLALSimNRTunedTidesFDTidalPhaseFrequencySeries` with `NRTidalv2_V` or
/// `NRTidalv3_V` at `frequencies` (Hz, increasing), for masses in kg and
/// the aligned spins (used by NRTidalv3), in either order.
///
/// The NRTidalv3 phase depends on the whole grid, as in LAL: from 0.9 times
/// the merger frequency it is held at its first minimum on the grid.
#[allow(clippy::too_many_arguments)]
pub fn tidal_series(
    frequencies: &[f64],
    m1_si: f64,
    m2_si: f64,
    lambda1: f64,
    lambda2: f64,
    chi1: f64,
    chi2: f64,
    version: NRTidalVersion,
) -> TidalSeries {
    // EnforcePrimaryMassIsm1_v3
    let (m1_si, m2_si, lambda1, lambda2, chi1, chi2) = if m1_si >= m2_si {
        (m1_si, m2_si, lambda1, lambda2, chi1, chi2)
    } else {
        (m2_si, m1_si, lambda2, lambda1, chi2, chi1)
    };
    let m1 = m1_si / MSUN_SI;
    let m2 = m2_si / MSUN_SI;
    let mtot = m1 + m2;
    let q = m1 / m2;
    let xa = m1 / mtot;
    let xb = m2 / mtot;
    let kappa2t = kappa2t(m1_si, m2_si, lambda1, lambda2);
    let n = frequencies.len();
    let mut series = TidalSeries {
        phase: Vec::with_capacity(n),
        amplitude: Vec::with_capacity(n),
        taper: Vec::with_capacity(n),
    };
    match version {
        NRTidalVersion::V2 => {
            let f_mrg = merger_frequency(mtot, kappa2t, q);
            let f_end_taper = 1.2 * f_mrg;
            for &f in frequencies {
                series.phase.push(tidal_phase_v2(f, xa, xb, mtot, kappa2t));
                series.amplitude.push(tidal_amplitude(f, mtot, kappa2t));
                series.taper.push(1.0 - planck_taper(f, f_mrg, f_end_taper));
            }
        }
        NRTidalVersion::V3 => {
            let f_mrg = merger_frequency_v3(mtot, lambda1, lambda2, q, chi1, chi2);
            let f_end_taper = 1.2 * f_mrg;
            let pn = pn_coeffs(xa);
            let c = v3_coeffs(xa, mtot, lambda1, lambda2, &pn);
            // From 0.9 f_merger, the phase stays at its first minimum.
            // (LAL compares the first point with the one before it, out of
            // bounds, if the grid starts above 0.9 f_merger; the first point
            // is never a minimum here.)
            let f_check = 0.9 * f_mrg;
            let mut minimum = None;
            for (i, &f) in frequencies.iter().enumerate() {
                let phase = tidal_phase_v3(f, mtot, &c, &pn);
                if i > 0 && f >= f_check && phase >= series.phase[i - 1] {
                    minimum = Some(series.phase[i - 1]);
                    break;
                }
                series.phase.push(phase);
            }
            if let Some(minimum) = minimum {
                series.phase.resize(n, minimum);
            }
            for (phase, &f) in series.phase.iter_mut().zip(frequencies) {
                // the smooth connection to the PN phase after merger, Eq. (45)
                let planck = planck_taper(f, 1.15 * f_mrg, 1.35 * f_mrg);
                *phase = *phase * (1.0 - planck)
                    + tidal_phase_pn(f, xa, mtot, lambda1, lambda2, &pn) * planck;
                series.amplitude.push(tidal_amplitude(f, mtot, kappa2t));
                series.taper.push(1.0 - planck_taper(f, f_mrg, f_end_taper));
            }
        }
    }
    series
}

/// `XLALSimInspiralEOSQfromLambda`: the quadrupole-monopole parameter
/// (`quadparam = 1 + dQuadMon`) of a neutron star from its tidal
/// deformability (arXiv:1303.1528), 1 (a black hole) below `Lambda = 0.5`.
pub fn quadparam_from_lambda(lambda: f64) -> f64 {
    if lambda < 5e-1 {
        1.0
    } else {
        // 0.194 + 0.0936 L + 0.0474 L^2 - 0.00421 L^3 + 0.000123 L^4 in
        // Horner's form (the value enters the SpinTaylor angles, whose start
        // is sensitive to its last bit)
        let loglam = lambda.ln();
        let q = fma(0.000123, loglam, -0.00421);
        let q = fma(q, loglam, 0.0474);
        let q = fma(q, loglam, 0.0936);
        fma(q, loglam, 0.194).exp()
    }
}

/// `XLALSimInspiralSetQuadMonParamsFromLambdas` for an empty `dQuadMon`:
/// the quadrupole-monopole parameter for `lambda > 0`, 1 otherwise.
pub(crate) fn default_quadparam(lambda: f64) -> f64 {
    // LAL stores dQuadMon = Q - 1 and the models read 1 + dQuadMon
    let dquadmon = if lambda > 0.0 {
        quadparam_from_lambda(lambda) - 1.
    } else {
        0.0
    };
    1. + dquadmon
}

/// `XLALSimUniversalRelationSpinInducedOctupoleVSSpinInducedQuadrupole`
/// (Yagi & Yunes, arXiv:1608.02582).
fn octupole_from_quadrupole(qm_def: f64) -> f64 {
    let coeffs = [0.003131, 2.071, -0.7152, 0.2458, -0.03309];
    let lnx = qm_def.ln();
    // XLALSimUniversalRelation
    let lny = coeffs[0]
        + coeffs[1] * lnx
        + coeffs[2] * lnx * lnx
        + coeffs[3] * lnx * lnx * lnx
        + coeffs[4] * lnx * lnx * lnx * lnx;
    lny.exp()
}

/// `XLALSimInspiralGetHOSpinTerms`: the 3.5PN spin-squared (tail) and
/// spin-cubed phase terms of the spin-induced quadrupole and octupole
/// moments (Eq. 27 of arXiv:1905.06011), `(SS_3p5PN, SSS_3p5PN)`.
pub(crate) fn ho_spin_terms(
    x_a: f64,
    x_b: f64,
    chi1: f64,
    chi2: f64,
    quadparam1: f64,
    quadparam2: f64,
) -> (f64, f64) {
    let x_asq = x_a * x_a;
    let x_bsq = x_b * x_b;
    let chi1_sq = chi1 * chi1;
    let chi2_sq = chi2 * chi2;
    // remove 1 for the black-hole baseline
    let octparam1 = octupole_from_quadrupole(quadparam1) - 1.;
    let octparam2 = octupole_from_quadrupole(quadparam2) - 1.;
    let ss = -400. * PI * (quadparam1 - 1.) * chi1_sq * x_asq
        - 400. * PI * (quadparam2 - 1.) * chi2_sq * x_bsq;
    let sss = 10.
        * ((x_asq + 308. / 3. * x_a) * chi1 + (x_bsq - 89. / 3. * x_b) * chi2)
        * (quadparam1 - 1.)
        * x_asq
        * chi1_sq
        + 10.
            * ((x_bsq + 308. / 3. * x_b) * chi2 + (x_asq - 89. / 3. * x_a) * chi1)
            * (quadparam2 - 1.)
            * x_bsq
            * chi2_sq
        - 440. * octparam1 * x_a * x_asq * chi1_sq * chi1
        - 440. * octparam2 * x_b * x_bsq * chi2_sq * chi2;
    (ss, sss)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn black_hole_limits() {
        assert_eq!(quadparam_from_lambda(0.0), 1.0);
        assert_eq!(default_quadparam(0.0), 1.0);
        // no spin-induced quadrupole term for black holes; the octupole fit
        // of LAL is not exactly 1 there, so the spin-cubed term is not zero
        let (ss, sss) = ho_spin_terms(0.6, 0.4, 0.3, 0.2, 1.0, 1.0);
        assert_eq!(ss, 0.0);
        assert!(sss != 0.0 && sss.abs() < 1e-2);
        assert_eq!(kappa2t(1.4 * MSUN_SI, 1.3 * MSUN_SI, 0.0, 0.0), 0.0);
    }

    #[test]
    fn kappa2t_is_symmetric() {
        let a = kappa2t(1.5 * MSUN_SI, 1.2 * MSUN_SI, 300.0, 800.0);
        let b = kappa2t(1.2 * MSUN_SI, 1.5 * MSUN_SI, 800.0, 300.0);
        assert_eq!(a, b);
    }

    #[test]
    fn v3_phase_is_held_at_its_minimum() {
        let frequencies: Vec<f64> = (0..4000).map(|i| 20.0 + i as f64).collect();
        let series = tidal_series(
            &frequencies,
            1.4 * MSUN_SI,
            1.3 * MSUN_SI,
            400.0,
            600.0,
            0.0,
            0.0,
            NRTidalVersion::V3,
        );
        assert_eq!(series.phase.len(), frequencies.len());
        assert_eq!(series.taper.last(), Some(&0.0));
        assert_eq!(series.taper[0], 1.0);
    }
}
