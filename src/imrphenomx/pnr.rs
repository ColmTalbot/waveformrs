// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The PNR ("phenomenological and numerical-relativity") machinery of
//! IMRPhenomXPNR in lalsimulation 6.2.0, following
//! `LALSimIMRPhenomX_PNR.c`, `LALSimIMRPhenomX_PNR_internals.c`,
//! `LALSimIMRPhenomX_PNR_alpha.c`, `LALSimIMRPhenomX_PNR_beta.c`,
//! `LALSimIMRPhenomX_AntisymmetricWaveform.c` and the PNR-related parts of
//! `LALSimIMRPhenomX_precession.c` (`IMRPhenomXGetAndSetPrecessionVariables`
//! and `IMRPhenomX_SetPrecessingRemnantParams` with `PhenomXPrecVersion`
//! 330 and `PhenomXPFinalSpinMod` 7, and the SpinTaylor angles of version
//! 330), optionally with the changes of lalsuite merge request !2593 to
//! the last two files (`parity_fix`):
//!
//! - [`PnrPrecession::new`] sets up the precession struct of XPNR: the
//!   SpinTaylor integration starting at the frequency the PNR angle
//!   interpolants need, the effective single-spin parameters, the
//!   PhenomXCP co-precessing deviations ([`CoprecDeviations`]), the
//!   precessing final spin with the sign of `cos(beta)` in the ringdown and
//!   the effective ringdown frequency, and the J frame;
//! - [`PnrPrecession::initialize_euler_angles`] builds the SpinTaylor
//!   `alpha` and `cos(beta)` of version 330 with their merger-ringdown
//!   continuation;
//! - [`PnrPrecession::angle_interpolants`] computes the PNR `alpha`, `beta`
//!   on a uniform grid covering the (mapped) frequencies of all modes and
//!   `gamma` from the minimal rotation condition, and returns their splines;
//! - [`PnrPrecession::remap_theta_jsf`] moves the line of sight to the J
//!   frame defined with the PNR `beta` at the reference frequency;
//! - [`PnrPrecession::antisymmetric_waveform`] is the antisymmetric (2,2)
//!   co-precessing mode of Ghosh et al. 2024 (arXiv:2310.16980).
//!
//! As in the rest of the crate, products and sums are fused where the
//! hardware has fused multiply-adds ([`crate::spintaylor::fp`]).

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use crate::constants::{C_SI, G_SI, MTSUN_SI};
use crate::grid::{FrequencyGrid, ScalePowers};
use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::hm_qnm;
use crate::imrphenomx::internals::{IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::msa::{MsaInputs, MsaSystem};
use crate::imrphenomx::phase::{
    phase_22, phase_22_connection_coefficients, time_shift_22, PhaseCoefficients,
};
use crate::imrphenomx::pnr_angles::{
    alpha_interpolation_0, alpha_interpolation_1, alpha_interpolation_2, alpha_interpolation_3,
    alphadot_cosbeta, angles_window, arctan_window, attach_mr_beta, beta_connection,
    beta_connection_frequencies, beta_rescaling_1, beta_rescaling_2, boole_step, chi_calc,
    coprec_window, intermediate_alpha, linear_frequency_map,
    linear_frequency_map_transition_frequencies, mr_alpha, mr_beta, pn_waveform_beta, rescale_beta,
    sqrt_one_minus_square, st_beta_coeffs_0, st_beta_coeffs_1, st_beta_coeffs_2, st_beta_coeffs_3,
    AlphaParams, BetaParams, PNR_CHI_WINDOW_LOWER, PNR_CHI_WINDOW_UPPER, PNR_DMF,
    PNR_Q_WINDOW_LOWER, PNR_Q_WINDOW_UPPER,
};
use crate::imrphenomx::pnr_fits;
use crate::imrphenomx::precession::{
    atan2_tol, delta_mf_from_wfstruct, hz_to_mf, mf_to_hz, nnlo_coefficients,
    precessing_final_spin_2017, rotate_y, rotate_z, NnloCoefficients, OrbitalAngularMomentum,
    MAX_TOL_ATAN,
};
use crate::imrphenomx::qnm::{fdamp22, fring22};
use crate::imrphenomx::spintaylor_angles::{
    inspiral_angles_spintaylor_with, unwrap_in_place, AlphaMrd, PnArrays,
};
use crate::imrphenomx::utilities::{chi_eff, final_spin_2017};
use crate::imrphenomx::wigner::SphericalHarmonics;
use crate::imrphenomx::{inspiral, intermediate, ringdown};
use crate::interpolate::CubicSpline;
use crate::spintaylor::fp::{fma, fms, fnms};

/// The PhenomXCP co-precessing deviations of IMRPhenomXPNR, stored in the
/// waveform struct: `PNR_DEV_PARAMETER` and the fits it multiplies (zero for
/// the other models), whether the tuned co-precessing model is in use
/// (`IMRPhenomXPNRUseTunedCoprec`: the higher modes then use the
/// precessing final spin and the effective ringdown frequency) and
/// `fRINGEffShiftDividedByEmm`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct CoprecDeviations {
    /// `PNR_DEV_PARAMETER = window a1 sin(theta_LS)`.
    pub(crate) dev: f64,
    /// Ringdown amplitude `v1RD`.
    pub(crate) mu1: f64,
    /// Ringdown amplitude `gamma3`.
    pub(crate) mu2: f64,
    /// Intermediate amplitude collocation value `V2`.
    pub(crate) mu3: f64,
    /// Ringdown phase `cL`.
    pub(crate) nu4: f64,
    /// (2,2) ringdown frequency.
    pub(crate) nu5: f64,
    /// (2,2) damping frequency.
    pub(crate) nu6: f64,
    /// Intermediate phase `b4`.
    pub(crate) zeta1: f64,
    /// Intermediate phase `b1`.
    pub(crate) zeta2: f64,
    /// `IMRPhenomXPNRUseTunedCoprec`.
    pub(crate) tuned: bool,
    /// `pnr_window`: one in the calibration region of the co-precessing
    /// model, zero far outside.
    pub(crate) pnr_window: f64,
    /// `fRINGEffShiftDividedByEmm`.
    pub(crate) f_ring_eff_shift_divided_by_emm: f64,
}

/// `IMRPhenomXHM_GenerateRingdownFrequency`: the ringdown frequency of mode
/// `(l, m)` (with the effective shift of the tuned co-precessing model).
pub(crate) fn ringdown_frequency(ell: u32, emm: u32, wf: &WaveformStruct) -> f64 {
    let afinal = if wf.coprec.tuned {
        wf.a_final_prec
    } else {
        wf.a_final
    };
    let f_ring = match (ell, emm) {
        (2, 1) => hm_qnm::fring21(afinal) / wf.m_final,
        (2, 2) => wf.f_ring,
        (3, 3) => hm_qnm::fring33(afinal) / wf.m_final,
        (3, 2) => hm_qnm::fring32(afinal) / wf.m_final,
        (4, 4) => hm_qnm::fring44(afinal) / wf.m_final,
        _ => unreachable!("mode ({ell},{emm}) is not available"),
    };
    if wf.coprec.tuned && (ell, emm) != (2, 2) {
        fnms(
            emm as f64,
            wf.coprec.f_ring_eff_shift_divided_by_emm,
            f_ring,
        )
    } else {
        f_ring
    }
}

/// The merger-ringdown continuation of the SpinTaylor `beta` of version 330
/// (`PhenomXPbetaMRD` with the offset set by the PNR `beta` in the
/// ringdown).
#[derive(Clone, Copy, Debug, Default)]
struct BetaMrd330 {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    cosbeta_sign: f64,
    flat: bool,
    dfdamp: f64,
}

impl BetaMrd330 {
    /// `betaMRD_coeff` for version 330: `off` is
    /// `IMRPhenomX_PNR_GenerateRingdownPNRBeta`.
    fn new(spline: &CubicSpline, fmax_pn: f64, wf: &WaveformStruct, off: f64) -> Self {
        let fdamp = wf.f_damp;
        let dfdamp = hm_qnm::fdamp21(wf.a_final) / wf.m_final - fdamp;
        let kappa = 2. * PI * dfdamp;
        let f1 = 0.97 * fmax_pn;
        let f2 = 0.98 * fmax_pn;
        let f1sq = f1 * f1;
        let f2sq = f2 * f2;
        let ef1 = (kappa * f1).exp();
        let ef2 = (kappa * f2).exp();
        let cosbeta1 = spline.eval(f1);
        let dcosbeta2 = if spline.contains(f2) {
            spline.derivative(f2)
        } else {
            f64::NAN
        };
        let cosbeta2 = spline.eval(f2);
        // As in LAL, only the status of the last evaluation is checked.
        let success = spline.contains(fmax_pn);
        let cosbetamax = spline.eval(fmax_pn);
        let mut beta = Self {
            dfdamp,
            cosbeta_sign: 1.0f64.copysign(cosbetamax),
            ..Default::default()
        };
        if cosbeta1.abs() > 1. || cosbeta2.abs() > 1. || cosbetamax.abs() > 1. || !success {
            beta.flat = true;
        } else {
            let beta1 = cosbeta1.acos();
            let beta2 = cosbeta2.acos();
            let sqrtarg = fnms(cosbeta2, cosbeta2, 1.);
            let dbeta2 = -dcosbeta2 / (if sqrtarg <= 0. { 1. } else { sqrtarg }).sqrt();
            let den = (f1 - f2).powf(2.);
            let f1_4 = f1.powf(4.);
            let f2_3 = f2.powf(3.);
            let f2_4 = f2.powf(4.);
            let k2 = fma(f2, kappa, 2.);
            let k3 = fma(f2, kappa, 3.);
            let k4 = fma(f2, kappa, 4.);
            let w_a = fma(-f2, k3, f1 * k4);
            let z_a = fma(f2 * (-f1 + f2), dbeta2, w_a * (off - beta2));
            beta.a = fma(-(ef1 * f1_4), off - beta1, ef2 * f2_3 * z_a) / den;
            let w_b = fma(-f2sq, k2, f1sq * k4);
            let z_b = fms((f1 - f2) * f2 * (f1 + f2), dbeta2, w_b * (off - beta2));
            beta.b = fma(2. * ef1 * f1_4 * f2, off - beta1, ef2 * f2_3 * z_b) / den;
            let w_c = fma(-f2, k2, f1 * k3);
            let z_c = fma(f2 * (-f1 + f2), dbeta2, w_c * (off - beta2));
            beta.c = fma(-(ef1 * f1_4 * f2sq), off - beta1, ef2 * f1 * f2_4 * z_c) / den;
            beta.d = off;
        }
        beta
    }

    /// `betaMRD`.
    #[inline]
    fn beta(&self, mf: f64) -> f64 {
        if self.flat {
            self.cosbeta_sign.acos()
        } else {
            let kappa = 2. * PI * self.dfdamp;
            let mf2 = mf * mf;
            let mf3 = mf2 * mf;
            fma(
                (-mf * kappa).exp() / mf,
                self.a / mf + self.b / mf2 + self.c / mf3,
                self.d,
            )
        }
    }
}

/// The SpinTaylor angles of version 330
/// (`IMRPhenomX_InterpolateAlphaBeta_SpinTaylor`).
#[derive(Clone, Debug)]
struct St330 {
    alpha_spline: CubicSpline,
    cosbeta_spline: CubicSpline,
    ftrans_mrd: f64,
    alpha_mrd: AlphaMrd,
    beta_mrd: BetaMrd330,
}

/// The splines of the PNR angles (`IMRPhenomX_PNR_angle_spline`) against
/// the frequency in Hz.
#[derive(Clone, Debug)]
pub(crate) struct PnrSplines {
    pub(crate) alpha: CubicSpline,
    pub(crate) beta: CubicSpline,
    pub(crate) gamma: CubicSpline,
}

impl PnrSplines {
    /// `gsl_spline_eval` of the three angles at `f` (Hz), which fails
    /// outside the knots.
    #[inline]
    pub(crate) fn angles(&self, f: f64) -> Result<(f64, f64, f64), IMRPhenomXError> {
        if !self.alpha.contains(f) {
            return Err(IMRPhenomXError::PnrAngles(
                "interpolation of the PNR angles out of range",
            ));
        }
        Ok((self.alpha.eval(f), self.beta.eval(f), self.gamma.eval(f)))
    }

    /// [`PnrSplines::angles`] searching the interval of `f` from `hint`,
    /// which it updates (the three splines share their knots; the result
    /// does not depend on the hint).
    #[inline]
    pub(crate) fn angles_from(
        &self,
        f: f64,
        hint: &mut usize,
    ) -> Result<(f64, f64, f64), IMRPhenomXError> {
        if !self.alpha.contains(f) {
            return Err(IMRPhenomXError::PnrAngles(
                "interpolation of the PNR angles out of range",
            ));
        }
        let index = self.alpha.locate(f, *hint);
        *hint = index;
        Ok((
            self.alpha.eval_in(index, f),
            self.beta.eval_in(index, f),
            self.gamma.eval_in(index, f),
        ))
    }
}

/// The modes the SpinTaylor and PNR angles must cover: `M_MIN`, `M_MAX`
/// and `L_MAX_PNR` of LAL's precession struct.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) struct ModeRange {
    pub(crate) m_min: u32,
    pub(crate) m_max: u32,
    pub(crate) l_max: u32,
}

impl ModeRange {
    /// `IMRPhenomX_GetandSetModes` and `L_MAX_PNR` for a co-precessing mode
    /// array `(l, m' > 0)` in the LALDict (both signs of `m'`).
    pub(crate) fn of_mode_array(modes: &[(u32, u32)]) -> Self {
        let l_max = if modes.contains(&(4, 4)) {
            4
        } else if modes.contains(&(3, 3)) || modes.contains(&(3, 2)) {
            3
        } else {
            2
        };
        Self {
            m_min: modes.iter().map(|&(_, m)| m).fold(4, u32::min),
            m_max: modes.iter().map(|&(_, m)| m).fold(1, u32::max),
            l_max,
        }
    }

    /// The single-mode functions without a mode array in the LALDict:
    /// `M_MIN = max(1, |m|)`, `M_MAX = L_MAX_PNR = l` for the J-frame mode
    /// `(l, m)`.
    pub(crate) fn of_jframe_mode(ell: u32, m: i32) -> Self {
        Self {
            m_min: m.unsigned_abs().max(1),
            m_max: ell,
            l_max: ell,
        }
    }
}

/// The quantities of `IMRPhenomXPrecessionStruct` used by IMRPhenomXPNR,
/// built by [`PnrPrecession::new`] (`IMRPhenomXGetAndSetPrecessionVariables`).
#[derive(Clone, Debug)]
pub(crate) struct PnrPrecession {
    /// `IMRPhenomXPrecVersion` after the fallbacks: 330 (SpinTaylor), 223
    /// (MSA, when the SpinTaylor integration fails) or 102 (NNLO, for
    /// negligible in-plane spins or when the MSA system fails).
    pub(crate) version: u32,
    /// First digit of the version when the precessing variables are set.
    pub(crate) precessing_tag: u32,
    /// `IMRPhenomXPNRUseTunedAngles` and `IMRPhenomXAntisymmetricWaveform`.
    pub(crate) use_tuned_angles: bool,
    pub(crate) antisymmetric: bool,
    /// Masses in units of the total mass and their squares.
    m1_2: f64,
    m2_2: f64,
    /// `pi G M / c^3` (s).
    pi_gm: f64,
    eta: f64,
    chi1: [f64; 3],
    chi2: [f64; 3],
    chi1_norm: f64,
    chi2_norm: f64,
    chi1_perp: f64,
    chi2_perp: f64,
    s1_perp: f64,
    s2_perp: f64,
    chi_tot_perp: f64,
    chi_p: f64,
    sl: f64,
    sperp: f64,
    /// The spins at the end of the SpinTaylor integration in the L frame.
    chi1_evolved: [f64; 3],
    chi2_evolved: [f64; 3],
    /// PNR single-spin quantities (`IMRPhenomX_PNR_GetAndSetPNRVariables`).
    chi_single_spin: f64,
    costheta_single_spin: f64,
    costheta_final_single_spin: f64,
    chi_single_spin_antisymmetric: f64,
    theta_antisymmetric: f64,
    pub(crate) pnr_inspiral_scaling: bool,
    pub(crate) pnr_hm_mflow: f64,
    pub(crate) pnr_hm_mfhigh: f64,
    /// `UseMRbeta` (state of the SpinTaylor `beta`).
    use_mr_beta: bool,
    /// Mode content: smallest and largest `|m|`, largest `l` for the PNR
    /// frequency range.
    m_min: u32,
    m_max: u32,
    l_max_pnr: u32,
    /// SpinTaylor integration (version 330).
    integration_buffer: f64,
    fmin_integration: f64,
    mfmin_integration: f64,
    pn_arrays: Option<PnArrays>,
    st: Option<St330>,
    /// Orbital angular momentum and J frame.
    l: OrbitalAngularMomentum,
    l_ref: f64,
    j0_norm: f64,
    theta_j_sf: f64,
    phi_j_sf: f64,
    kappa: f64,
    alpha0: f64,
    pub(crate) epsilon0: f64,
    pub(crate) theta_jn: f64,
    pub(crate) zeta_polarization: f64,
    nnlo: NnloCoefficients,
    /// The MSA system of the binary (version 223).
    msa: Option<MsaSystem>,
    /// Offsets of `alpha` and `epsilon` (updated with the PNR angles).
    pub(crate) alpha_offset: f64,
    pub(crate) epsilon_offset: f64,
    pub(crate) ylm: SphericalHarmonics,
    /// `pWF->betaRD`.
    beta_rd: f64,
    /// The code of lalsuite !2593 (LALSimulation master, not 6.2.0 or
    /// 6.2.1): the phase of the antisymmetric (2,2) mode referred to
    /// `alpha(f_ref)` and the azimuth of `J` rather than to `alpha_offset`
    /// and the polarization angle, and the evolved spins rotated into the
    /// L frame with the components of `L` rather than `atan2` and `acos`.
    parity_fix: bool,
}

/// The effective single-spin binary of the two-spin map
/// (`pPrec_SingleSpin`, set up with `PhenomXPrecVersion` 223): only its MSA
/// angles are used. If its MSA system sets `MSA_ERROR` the struct carries
/// version 102, and `cos(theta_L)` then uses the 3PN `L` of that version.
struct SingleSpin {
    msa: MsaSystem,
    l_102: Option<OrbitalAngularMomentum>,
}

impl SingleSpin {
    /// `acos(cos(theta_L))` of `IMRPhenomX_Return_phi_zeta_costhetaL_MSA`.
    fn beta(&self, v: f64) -> f64 {
        match &self.l_102 {
            None => self.msa.phi_zeta_costhetal(v).cos_theta_l.acos(),
            Some(l) => self.msa.cos_theta_l_with(v, l).acos(),
        }
    }
}

/// The single-spin system of the two-spin map and the PNR angle parameters
/// (`IMRPhenomX_PNR_PopulateStructs`).
struct PnrStructs {
    single: Option<SingleSpin>,
    alpha: AlphaParams,
    beta: BetaParams,
}

/// `sqrt(x^2 + y^2)` as contracted by clang.
#[inline]
fn hypot2(x: f64, y: f64) -> f64 {
    fma(x, x, y * y).sqrt()
}

/// The in-plane spin magnitude `sqrt(chi1x^2 + chi1y^2 + chi2x^2 + chi2y^2)`.
#[inline]
fn chi_in_plane(chi1: [f64; 3], chi2: [f64; 3]) -> f64 {
    fma(
        chi2[1],
        chi2[1],
        fma(chi2[0], chi2[0], fma(chi1[0], chi1[0], chi1[1] * chi1[1])),
    )
    .sqrt()
}

/// The spins `s1`, `s2` in the frame of `L = (lx, ly, lz)` along `z`, as
/// lalsuite !2593 rotates the evolved spins of version 330: by `-phi` about
/// `z` and `-theta` about `y`, with the cosines and sines taken from the
/// components of `L` (`phi = 0` for `L` along `z`). Unlike `atan2`, this
/// commutes exactly with the reflection `(x, y) -> (-x, -y)` of `L` and the
/// spins.
fn rotate_to_l_frame(
    lx: f64,
    ly: f64,
    lz: f64,
    s1: [f64; 3],
    s2: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let lperp_sq = fma(lx, lx, ly * ly);
    let lperp = lperp_sq.sqrt();
    let lnorm = fma(lz, lz, lperp_sq).sqrt();
    let (cos_theta, sin_theta) = (lz / lnorm, lperp / lnorm);
    let (cos_phi, sin_phi) = if lperp > 0.0 {
        (lx / lperp, ly / lperp)
    } else {
        (1.0, 0.0)
    };
    let rotate = |s: [f64; 3]| {
        let xp = fma(s[0], cos_phi, s[1] * sin_phi);
        let yp = fnms(s[0], sin_phi, s[1] * cos_phi);
        let zp = s[2];
        [
            fma(xp, cos_theta, -(zp * sin_theta)),
            yp,
            fma(xp, sin_theta, cos_theta * zp),
        ]
    };
    (rotate(s1), rotate(s2))
}

fn pnr_error(why: &'static str) -> IMRPhenomXError {
    IMRPhenomXError::PnrAngles(why)
}

impl PnrPrecession {
    /// `IMRPhenomXGetAndSetPrecessionVariables` for the flags set by
    /// IMRPhenomXPNR (`PhenomXPrecVersion` 330, `PhenomXPFinalSpinMod` 7,
    /// `PNRUseTunedAngles`, `PNRUseTunedCoprec`, `AntisymmetricWaveform`,
    /// `PhenomXPConvention` 1, multibanding off), for the range of modes
    /// `range`. Modifies `wf` as LAL
    /// modifies `pWF` (`phi0`, `deltaMF`, the co-precessing deviations, the
    /// final spin and the (2,2) ringdown frequencies).
    ///
    /// `m1_si`, `m2_si` (kg) and the spins must be ordered with the heavier
    /// body first (`XLALIMRPhenomXPCheckMassesAndSpins`). `parity_fix`
    /// selects the code of lalsuite !2593 (LALSimulation master) over that
    /// of LALSimulation 6.2.0.
    pub(crate) fn new(
        wf: &mut WaveformStruct,
        m1_si: f64,
        m2_si: f64,
        chi1: [f64; 3],
        chi2: [f64; 3],
        range: ModeRange,
        parity_fix: bool,
    ) -> Result<Self, IMRPhenomXError> {
        let [chi1x, chi1y, chi1z] = chi1;
        let [chi2x, chi2y, chi2z] = chi2;

        let chi_in_plane = chi_in_plane(chi1, chi2);
        let version = if chi_in_plane < 1e-6 { 102 } else { 330 };
        // disable the tuned angles, the tuned co-precessing model and the
        // mode asymmetries for negligible in-plane spins
        let pnr_on = chi_in_plane >= 1e-7;
        // IMRPhenomXSetWaveformVariables reads PNRUseTunedCoprec
        wf.coprec = CoprecDeviations {
            tuned: true,
            ..Default::default()
        };

        let m1 = m1_si / wf.m_tot_si;
        let m2 = m2_si / wf.m_tot_si;
        let big_m = m1 + m2;
        let m1_2 = m1 * m1;
        let m1_3 = m1 * m1_2;
        let m1_4 = m1 * m1_3;
        let m1_5 = m1 * m1_4;
        let m1_6 = m1 * m1_5;
        let m1_7 = m1 * m1_6;
        let m1_8 = m1 * m1_7;
        let m2_2 = m2 * m2;
        let q = m1 / m2;
        let eta = wf.eta;
        let delta = wf.delta;
        let pi_gm = PI * (m1_si + m2_si) * (G_SI / C_SI) / (C_SI * C_SI);

        let chi1_norm = fma(chi1z, chi1z, fma(chi1x, chi1x, chi1y * chi1y)).sqrt();
        let chi2_norm = fma(chi2z, chi2z, fma(chi2x, chi2x, chi2y * chi2y)).sqrt();
        for norm in [chi1_norm, chi2_norm] {
            if norm.is_nan() || norm.abs() > 1.0 {
                return Err(IMRPhenomXError::UnphysicalSpinMagnitude(norm));
            }
        }
        let s1x = chi1x * m1_2;
        let s1y = chi1y * m1_2;
        let s2x = chi2x * m2_2;
        let s2y = chi2y * m2_2;
        let s_tot_perp = hypot2(s1x + s2x, s1y + s2y);
        let chi_tot_perp = s_tot_perp * (big_m * big_m) / m1_2;

        // effective precessing spin (Schmidt et al. 2015)
        let s1_perp = m1_2 * hypot2(chi1x, chi1y);
        let s2_perp = m2_2 * hypot2(chi2x, chi2y);
        let a1 = 2.0 + (3.0 * m2) / (2.0 * m1);
        let a2 = 2.0 + (3.0 * m1) / (2.0 * m2);
        let asp1 = a1 * s1_perp;
        let asp2 = a2 * s2_perp;
        let num = if asp2 > asp1 { asp2 } else { asp1 };
        let den = if m2 > m1 { a2 * m2_2 } else { a1 * m1_2 };
        let chip = num / den;

        let ModeRange {
            m_min,
            m_max,
            l_max: l_max_pnr,
        } = range;

        let mut prec = Self {
            version,
            precessing_tag: version / 100,
            use_tuned_angles: pnr_on,
            antisymmetric: pnr_on,
            m1_2,
            m2_2,
            pi_gm,
            eta,
            chi1,
            chi2,
            chi1_norm,
            chi2_norm,
            chi1_perp: hypot2(chi1x, chi1y),
            chi2_perp: hypot2(chi2x, chi2y),
            s1_perp,
            s2_perp,
            chi_tot_perp,
            chi_p: chip,
            sl: fma(chi1z, m1_2, chi2z * m2_2),
            sperp: chip * m1_2,
            chi1_evolved: chi1,
            chi2_evolved: chi2,
            chi_single_spin: 0.0,
            costheta_single_spin: 0.0,
            costheta_final_single_spin: 0.0,
            chi_single_spin_antisymmetric: 0.0,
            theta_antisymmetric: 0.0,
            pnr_inspiral_scaling: false,
            pnr_hm_mflow: 0.0,
            pnr_hm_mfhigh: 0.0,
            use_mr_beta: false,
            m_min,
            m_max,
            l_max_pnr,
            integration_buffer: 0.0,
            fmin_integration: 0.0,
            mfmin_integration: 0.0,
            pn_arrays: None,
            st: None,
            l: OrbitalAngularMomentum::three_pn(eta, delta, chi1z, chi2z),
            l_ref: 0.0,
            j0_norm: 0.0,
            theta_j_sf: 0.0,
            phi_j_sf: 0.0,
            kappa: 0.0,
            alpha0: 0.0,
            epsilon0: 0.0,
            theta_jn: 0.0,
            zeta_polarization: 0.0,
            nnlo: NnloCoefficients::default(),
            msa: None,
            alpha_offset: 0.0,
            epsilon_offset: 0.0,
            ylm: SphericalHarmonics::default(),
            beta_rd: 0.0,
            parity_fix,
        };

        // ---- SpinTaylor set-up (version 330) ----
        if prec.version / 100 == 3 {
            prec.spintaylor_setup(wf)?;
        }

        // update the precessing version to catch the fallbacks
        prec.precessing_tag = prec.version / 100;

        // the two-spin to single-spin map used by PNR and XCP
        prec.set_pnr_variables(wf, prec.version);

        // IMRPhenomX_PNR_GetAndSetCoPrecParams
        prec.set_coprec_params(wf, pnr_on);

        if prec.version == 223 {
            let msa = MsaSystem::new_unchecked(&prec.msa_inputs(wf, chi1, chi2));
            if msa.msa_error() {
                // default to the NNLO angles with the 3PN aligned-spin L
                prec.version = 102;
            } else {
                prec.msa = Some(msa);
            }
        }

        // IMRPhenomX_SetPrecessingRemnantParams
        prec.set_remnant(wf);

        // orbital angular momentum
        prec.l = if prec.version == 223 {
            OrbitalAngularMomentum::three_pn_msa(eta, delta, chi1z, chi2z)
        } else {
            OrbitalAngularMomentum::three_pn(eta, delta, chi1z, chi2z)
        };
        prec.l_ref = big_m * big_m * prec.l.ansatz(wf.v_ref, wf.eta / wf.v_ref);

        // source frame J = L + S1 + S2
        let j0 = [
            fma(m1_2, chi1x, m2_2 * chi2x),
            fma(m1_2, chi1y, m2_2 * chi2y),
            fma(m1_2, chi1z, m2_2 * chi2z) + prec.l_ref,
        ];
        prec.j0_norm = fma(j0[2], j0[2], fma(j0[0], j0[0], j0[1] * j0[1])).sqrt();
        prec.theta_j_sf = if prec.j0_norm < 1e-10 {
            0.0
        } else {
            (j0[2] / prec.j0_norm).acos()
        };
        // convention 1
        prec.phi_j_sf = if j0[0].abs() < MAX_TOL_ATAN && j0[1].abs() < MAX_TOL_ATAN {
            0.0
        } else {
            j0[1].atan2(j0[0])
        };
        wf.phi0 = 0.0;

        let n_sf = Self::line_of_sight(wf);
        let tmp = rotate_y(-prec.theta_j_sf, rotate_z(-prec.phi_j_sf, n_sf));
        prec.kappa = atan2_tol(tmp[1], tmp[0], MAX_TOL_ATAN);
        prec.alpha0 = PI - prec.kappa;
        let j0_dot_n = fma(j0[2], n_sf[2], fma(j0[0], n_sf[0], j0[1] * n_sf[1]));
        prec.theta_jn = (j0_dot_n / prec.j0_norm).acos();
        prec.zeta_polarization =
            prec.polarization_angle(wf, prec.theta_jn.sin(), prec.theta_jn.cos());

        if prec.version / 100 == 1 {
            prec.nnlo = nnlo_coefficients(
                eta, delta, q, wf.chi_eff, chip, m1, m1_2, m1_3, m1_4, m1_5, m1_6, m1_8,
            );
        }
        prec.epsilon0 = prec.phi_j_sf - PI;

        // Get_alphaepsilon_atfref with mprime = 2
        let omega_ref = wf.pi_m * wf.f_ref * 2. / 2.;
        (prec.alpha_offset, prec.epsilon_offset) = match &prec.msa {
            Some(msa) => {
                let angles = msa.phi_zeta_costhetal(omega_ref.cbrt());
                (angles.phiz - prec.alpha0, angles.zeta - prec.epsilon0)
            }
            None => {
                let logomega_ref = omega_ref.ln();
                let omega_ref_cbrt = omega_ref.cbrt();
                let omega_ref_cbrt2 = omega_ref_cbrt * omega_ref_cbrt;
                let (alpha, epsilon) =
                    prec.nnlo
                        .evaluate(omega_ref, logomega_ref, omega_ref_cbrt, omega_ref_cbrt2);
                (alpha - prec.alpha0, epsilon - prec.epsilon0)
            }
        };
        prec.ylm = SphericalHarmonics::new(prec.theta_jn);
        let _ = m1_7;
        Ok(prec)
    }

    /// The direction of the line of sight in the source frame.
    fn line_of_sight(wf: &WaveformStruct) -> [f64; 3] {
        let inclination = wf.inclination;
        let phi_ref = wf.phi_ref_in;
        [
            inclination.sin() * (FRAC_PI_2 - phi_ref).cos(),
            inclination.sin() * (FRAC_PI_2 - phi_ref).sin(),
            inclination.cos(),
        ]
    }

    /// The polarization angle `zeta` from the current J frame (`thetaJ_Sf`,
    /// `phiJ_Sf`, `kappa`) and the components `(Nx, Nz)` of the line of
    /// sight in the J frame (convention 1).
    fn polarization_angle(&self, wf: &WaveformStruct, nx_jf: f64, nz_jf: f64) -> f64 {
        let inclination = wf.inclination;
        let phi_ref = wf.phi_ref_in;
        let x_sf = [
            -inclination.cos() * phi_ref.sin(),
            -inclination.cos() * phi_ref.cos(),
            inclination.sin(),
        ];
        let tmp = rotate_z(
            -self.kappa,
            rotate_y(-self.theta_j_sf, rotate_z(-self.phi_j_sf, x_sf)),
        );
        let p_arun = [nz_jf, 0.0, -nx_jf];
        let q_arun = [0.0, 1.0, 0.0];
        let x_dot_p = fma(
            tmp[2],
            p_arun[2],
            fma(tmp[0], p_arun[0], tmp[1] * p_arun[1]),
        );
        let x_dot_q = fma(
            tmp[2],
            q_arun[2],
            fma(tmp[0], q_arun[0], tmp[1] * q_arun[1]),
        );
        x_dot_q.atan2(x_dot_p)
    }

    /// Inputs of the MSA system for spins `chi1`, `chi2`.
    fn msa_inputs(&self, wf: &WaveformStruct, chi1: [f64; 3], chi2: [f64; 3]) -> MsaInputs {
        MsaInputs {
            eta: wf.eta,
            m1: wf.m1,
            m2: wf.m2,
            m1_2: self.m1_2,
            m2_2: self.m2_2,
            chi1,
            chi2,
            pi_gm: self.pi_gm,
            f_ref: wf.f_ref,
        }
    }

    /// The SpinTaylor part of `IMRPhenomXGetAndSetPrecessionVariables` with
    /// the PNR angles: the start of the integration needed by the PNR angle
    /// interpolants, the integration and the evolved spins; falls back to
    /// version 223 if the integration fails.
    fn spintaylor_setup(&mut self, wf: &mut WaveformStruct) -> Result<(), IMRPhenomXError> {
        let mut flow = wf.f_min;
        if wf.delta_f == 0. {
            wf.delta_mf = delta_mf_from_wfstruct(wf);
        }
        if !self.use_tuned_angles {
            self.integration_buffer = if wf.delta_f > 0. {
                3. * wf.delta_f
            } else {
                0.5
            };
            flow = (wf.f_min - self.integration_buffer) * 2. / self.m_max as f64;
        } else {
            if wf.delta_f != 0. {
                let i_start_here = (wf.f_min / wf.delta_f) as usize;
                flow = i_start_here as f64 * wf.delta_f;
            }
            let fmin_hm_inspiral = flow * 2.0 / self.m_max as f64;
            // the PNR variables and the alpha and beta coefficients with a
            // fake version 223 to catch edge cases
            self.set_pnr_variables(wf, 223);
            let alpha_params = self.precompute_alpha_coefficients(wf, 223);
            let mut beta_params = self.precompute_beta_coefficients(wf, 223);
            beta_connection_frequencies(&mut beta_params);
            let mf_alpha_upper = alpha_params.a4 / 3.0;
            let mut mf_low_cut = (3.0 / 3.5) * mf_alpha_upper;
            let mut mf_high_cut = beta_params.mf_beta_lower;
            if mf_high_cut > wf.f_cut_def || mf_high_cut < 0.1 * wf.f_ring {
                mf_high_cut = wf.f_ring;
            }
            if mf_low_cut > wf.f_cut_def || mf_high_cut < mf_low_cut {
                mf_low_cut = mf_high_cut / 2.0;
            }
            let flow_alpha = mf_to_hz(mf_low_cut * 0.65 * self.m_max as f64 / 2.0, wf.m_tot);
            if flow_alpha < flow {
                // flow is approximately in the intermediate region of the
                // frequency map: conservatively reduce it
                flow = fmin_hm_inspiral / 1.5;
            } else {
                let mf_rd_22 = wf.f_ring;
                let mf_rd_lm = ringdown_frequency(self.l_max_pnr, self.m_max, wf);
                let fmin_hm_ringdown =
                    mf_to_hz(hz_to_mf(flow, wf.m_tot) - (mf_rd_lm - mf_rd_22), wf.m_tot);
                flow = if fmin_hm_ringdown < fmin_hm_inspiral && fmin_hm_ringdown > 0.0 {
                    fmin_hm_ringdown
                } else {
                    fmin_hm_inspiral
                };
            }
            let pnr_interpolation_deltaf = self.interpolation_delta_f(wf, flow);
            self.integration_buffer = 1.4 * pnr_interpolation_deltaf;
            flow = if fnms(2.0, pnr_interpolation_deltaf, flow) < 0. {
                flow / 2.0
            } else {
                fnms(2.0, pnr_interpolation_deltaf, flow)
            };
            let i_start_here = (flow / pnr_interpolation_deltaf) as usize;
            flow = i_start_here as f64 * pnr_interpolation_deltaf;
        }
        if flow.is_nan() || flow <= 0. {
            return Err(pnr_error(
                "starting frequency for SpinTaylor angles must be positive",
            ));
        }
        self.integrate_spintaylor(wf, flow);
        Ok(())
    }

    /// `IMRPhenomX_InspiralAngles_SpinTaylor` from `flow` (Hz) to
    /// `fRING + 8 fDAMP` of the current waveform struct; on success stores
    /// the PN solution, the start of the integration and the evolved spins,
    /// otherwise switches to version 223.
    fn integrate_spintaylor(&mut self, wf: &WaveformStruct, flow: f64) {
        match inspiral_angles_spintaylor_with(
            wf,
            self.chi1,
            self.chi2,
            330,
            wf.f_ring + 8. * wf.f_damp,
            flow,
        ) {
            Ok(arrays) => {
                self.fmin_integration = arrays.v[0].powf(3.) / self.pi_gm;
                self.mfmin_integration = hz_to_mf(self.fmin_integration, wf.m_tot);
                // the spins at the end of the integration in the L frame
                let last = arrays.len() - 1;
                let (lx, ly, lz) = (
                    arrays.lnhat[0][last],
                    arrays.lnhat[1][last],
                    arrays.lnhat[2][last],
                );
                (self.chi1_evolved, self.chi2_evolved) = if self.parity_fix {
                    rotate_to_l_frame(lx, ly, lz, arrays.s1_end, arrays.s2_end)
                } else {
                    let phi = ly.atan2(lx);
                    let theta = (lz / fma(lz, lz, fma(lx, lx, ly * ly)).sqrt()).acos();
                    let evolve = |s: [f64; 3]| rotate_y(-theta, rotate_z(-phi, s));
                    (evolve(arrays.s1_end), evolve(arrays.s2_end))
                };
                self.pn_arrays = Some(arrays);
            }
            Err(_) => {
                // the model defaults to the MSA angles
                self.pn_arrays = None;
                self.version = 223;
            }
        }
    }

    /// `IMRPhenomX_PNR_CheckTwoSpin`.
    fn two_spin(&self) -> bool {
        self.chi1_norm != 0.0 && self.chi2_norm >= 1.0e-3 && self.version >= 200
    }

    /// `IMRPhenomX_PNR_GetAndSetPNRVariables` with `IMRPhenomXPrecVersion`
    /// equal to `version`.
    fn set_pnr_variables(&mut self, wf: &WaveformStruct, version: u32) {
        let m1 = wf.m1 * wf.m_tot;
        let m2 = wf.m2 * wf.m_tot;
        let q = wf.q;
        let (chi1, chi2) = if version == 330 {
            (self.chi1_evolved, self.chi2_evolved)
        } else {
            (self.chi1, self.chi2)
        };
        let [chi1x, chi1y, chi1z] = chi1;
        let [chi2x, chi2y, chi2z] = chi2;
        let chieff = chi_eff(wf.eta, chi1z, chi2z);
        let chipar = wf.m_tot * chieff / m1;
        let m1sq = m1 * m1;
        let m2sq = m2 * m2;
        let chis = || {
            hypot2(
                fma(m1sq, chi1x, m2sq * chi2x),
                fma(m1sq, chi1y, m2sq * chi2y),
            ) / (m1 * m1)
        };
        let s = ((q - 1.0) * PI).sin();
        let c = ((q - 1.0) * PI).cos();
        let chiperp = if version == 330 {
            chis()
        } else if q <= 1.5 {
            fma(s * s, self.chi_p, c * c * chis())
        } else {
            self.chi_p
        };
        let chiperp_antisymmetric = if q <= 1.5 {
            let antisymmetric_chis = hypot2(
                fma(m1sq, chi1x, -(m2sq * chi2x)),
                fma(m1sq, chi1y, -(m2sq * chi2y)),
            ) / (m1 * m1);
            fma(s * s, self.chi_p, c * c * antisymmetric_chis)
        } else {
            self.chi_p
        };
        let chi_mag = hypot2(chipar, chiperp);
        self.chi_single_spin = chi_mag;
        let chi_mag_antisymmetric = hypot2(chipar, chiperp_antisymmetric);
        self.chi_single_spin_antisymmetric = chi_mag_antisymmetric;
        let costheta = if chi_mag >= 1.0e-6 {
            chipar / chi_mag
        } else {
            0.
        };
        self.theta_antisymmetric = if chi_mag_antisymmetric >= 1.0e-6 {
            (chipar / chi_mag_antisymmetric).acos()
        } else {
            0.0
        };
        self.costheta_single_spin = costheta;

        // approximate final spin from the single-spin mapping
        let xfparr = final_spin_2017(wf.eta, chi_mag * costheta, 0.0);
        let qfactor = q / (1.0 + q);
        let xfperp = qfactor * qfactor * chi_mag * sqrt_one_minus_square(costheta);
        let xf = hypot2(xfparr, xfperp);
        self.costheta_final_single_spin = if xf > 1.0e-6 { xfparr / xf } else { 0. };

        self.pnr_hm_mflow = 0.0;
        self.pnr_hm_mfhigh = 0.0;
        self.pnr_inspiral_scaling =
            q > PNR_Q_WINDOW_UPPER || self.chi_single_spin > PNR_CHI_WINDOW_UPPER;
    }

    /// `IMRPhenomX_PNR_GetAndSetCoPrecParams` (`PNRUseTunedCoprec` equal to
    /// `tuned`).
    fn set_coprec_params(&self, wf: &mut WaveformStruct, tuned: bool) {
        let a1 = self.chi_single_spin;
        let theta_ls = self.costheta_single_spin.acos();
        let pnr_window = if tuned { coprec_window(wf.q) } else { 0.0 };
        let apply = if tuned { 1.0 } else { 0.0 };
        let mut dev = a1 * theta_ls.sin() * apply;
        if tuned {
            dev *= pnr_window;
        }
        let mut coprec = CoprecDeviations {
            dev,
            tuned,
            pnr_window,
            ..Default::default()
        };
        if tuned {
            let coprec_eta;
            let mut coprec_a1;
            if self.version == 330 {
                // flatten the dependences to limit extrapolation artefacts
                coprec_eta = if wf.eta >= 0.09876 {
                    wf.eta
                } else {
                    fma(-(0.09876 - wf.eta), 0.1641, 0.09876)
                };
                coprec_a1 = if a1 <= 0.8 {
                    a1
                } else {
                    0.8 + (a1 - 0.8) / 12.0
                };
            } else {
                coprec_eta = if wf.eta >= 0.09876 { wf.eta } else { 0.09876 };
                coprec_a1 = if a1 <= 0.8 { a1 } else { 0.8 };
            }
            coprec_a1 = if coprec_a1 >= 0.2 { coprec_a1 } else { 0.2 };
            let (e, a) = (coprec_eta, coprec_a1);
            coprec.mu1 = pnr_fits::mu1_l2m2(theta_ls, e, a);
            coprec.mu2 = pnr_fits::mu2_l2m2(theta_ls, e, a);
            coprec.mu3 = pnr_fits::mu3_l2m2(theta_ls, e, a);
            coprec.nu4 = pnr_fits::nu4_l2m2(theta_ls, e, a);
            coprec.nu5 = pnr_fits::nu5_l2m2(theta_ls, e, a);
            coprec.nu6 = pnr_fits::nu6_l2m2(theta_ls, e, a);
            coprec.zeta1 = pnr_fits::zeta1_l2m2(theta_ls, e, a);
            coprec.zeta2 = pnr_fits::zeta2_l2m2(theta_ls, e, a);
        }
        wf.coprec = coprec;
    }

    /// `IMRPhenomX_PNR_GenerateRingdownPNRBeta`: the PNR `beta` in the
    /// ringdown.
    fn ringdown_pnr_beta(&self, wf: &WaveformStruct) -> f64 {
        arctan_window(
            self.costheta_final_single_spin.acos()
                - pnr_fits::beta_bf(wf.eta, self.chi_single_spin, self.costheta_single_spin),
        )
    }

    /// `IMRPhenomX_SetPrecessingRemnantParams` with `PhenomXPFinalSpinMod`
    /// 7: the precessing final spin with the sign of `cos(betaRD)`, the
    /// final spin of the co-precessing model (non-precessing in its
    /// calibration region), the deviations of the (2,2) QNM frequencies and
    /// the effective ringdown frequency.
    fn set_remnant(&mut self, wf: &mut WaveformStruct) {
        let tuned = wf.coprec.tuned;
        if tuned {
            self.beta_rd = self.ringdown_pnr_beta(wf);
        }
        let sign = 1.0f64.copysign(self.beta_rd.cos());
        let afinal_prec =
            precessing_final_spin_2017(wf.eta, self.chi1[2], self.chi2[2], self.chi_tot_perp);
        wf.a_final_prec = sign * afinal_prec.abs();
        let window = wf.coprec.pnr_window;
        wf.a_final = if tuned {
            fma(window, wf.a_final_nonprec, (1.0 - window) * wf.a_final_prec)
        } else {
            wf.a_final_prec
        };
        if wf.a_final.abs() > 1.0 {
            wf.a_final = 1.0f64.copysign(wf.a_final);
        }
        if wf.a_final_prec.abs() > 1.0 {
            wf.a_final_prec = 1.0f64.copysign(wf.a_final_prec);
        }
        wf.f_ring = fring22(wf.a_final) / wf.m_final;
        wf.f_damp = fdamp22(wf.a_final) / wf.m_final;
        if tuned {
            // APPLY_PNR_DEVIATIONS
            let c = wf.coprec;
            wf.f_ring = fnms(c.dev, c.nu5, wf.f_ring);
            wf.f_damp = fma(c.dev, c.nu6, wf.f_damp);
            // the effective ringdown frequency
            let f_ring22_prec = fring22(wf.a_final_prec) / wf.m_final;
            let f_ring21_prec = hm_qnm::fring21(wf.a_final_prec) / wf.m_final;
            let shift = (1.0 - self.beta_rd.cos().abs()) * (f_ring22_prec - f_ring21_prec);
            wf.coprec.f_ring_eff_shift_divided_by_emm = shift;
            wf.f_ring = fnms((1.0 - window) * 2.0, shift, wf.f_ring);
        }
    }

    /// `IMRPhenomX_PNR_HMInterpolationDeltaF`: the frequency spacing (Hz)
    /// of the uniform grid of the PNR angle interpolants.
    fn interpolation_delta_f(&self, wf: &WaveformStruct, f_min: f64) -> f64 {
        let error_tolerance = 0.01;
        // aligned spins: no resolution needed for the (vanishing) angles
        if self.chi1_perp == 0.0 && self.chi2_perp == 0.0 {
            return if wf.delta_f != 0. { wf.delta_f } else { 0.1 };
        }
        // leading-order scaling of alpha
        let mf = hz_to_mf(f_min, wf.m_tot);
        let eta_term = fnms(4.0, wf.eta, 1.0).sqrt();
        let numerator = 3.0 * PI * mf * mf * mf * mf * mf * error_tolerance * (1.0 + eta_term);
        let denominator = fma(13.0, eta_term, 7.0);
        let constant = 4.0 * (2.0f64 / 5.0).sqrt();
        let delta_f_alpha = mf_to_hz(constant * (numerator / denominator).sqrt().sqrt(), wf.m_tot);

        // two-spin oscillations: sample the period of the spin magnitude
        if self.two_spin() {
            let (l, terms) = if self.version / 100 == 3 {
                // temporarily initialise the MSA system with the 3PN L
                let msa = MsaSystem::new_unchecked(&self.msa_inputs(wf, self.chi1, self.chi2));
                (
                    OrbitalAngularMomentum::three_pn(
                        self.eta,
                        wf.delta,
                        self.chi1[2],
                        self.chi2[2],
                    ),
                    msa.psi_terms(),
                )
            } else {
                let msa = self.msa.as_ref().expect("MSA system");
                (self.l, msa.psi_terms())
            };
            let (g0, deltam, psi1, psi2) = terms;
            let v0 = (PI * mf).cbrt();
            let v02 = v0 * v0;
            let v06 = v02 * v02 * v02;
            let dpsi = g0 * deltam * PI / (4.0 * v06) * fma(psi2, v02, fma(2.0 * psi1, v0, 3.0));
            let mut dpsi_inv = (1.0 / dpsi).abs();
            let l_fmin = wf.m_tot * wf.m_tot * l.ansatz(v0, wf.eta / v0);
            let beta_min = (self.s1_perp - self.s2_perp).abs().atan2(l_fmin + self.sl);
            let beta_max = (self.s1_perp + self.s2_perp).atan2(l_fmin + self.sl);
            if beta_min < 0.01 && beta_min / beta_max < 0.55 {
                dpsi_inv /= 4.0;
            }
            // four points per oscillation
            let delta_f_twospin = mf_to_hz(dpsi_inv / 4.0, wf.m_tot);
            if delta_f_twospin < delta_f_alpha && !dpsi.is_nan() {
                return if delta_f_twospin < 1e-2 {
                    1e-2
                } else {
                    delta_f_twospin
                };
            }
        }
        if delta_f_alpha < 1e-2 {
            1e-2
        } else {
            delta_f_alpha
        }
    }

    /// The single-spin `eta` and `chi` of the MR fits.
    fn fit_parameters(&self, wf: &WaveformStruct, version: u32) -> (f64, f64) {
        if version == 330 {
            let eta = if wf.eta >= 0.09 { self.eta } else { 0.09 };
            let chiboundary = fma(-0.20, (-((wf.q - 6.0) / 1.5).powf(8.0)).exp(), 0.80);
            let chi = if self.chi_single_spin <= chiboundary {
                self.chi_single_spin
            } else {
                chiboundary
            };
            (eta, chi)
        } else {
            (wf.eta, self.chi_single_spin)
        }
    }

    /// `IMRPhenomX_PNR_precompute_alpha_coefficients`.
    fn precompute_alpha_coefficients(&self, wf: &WaveformStruct, version: u32) -> AlphaParams {
        let (eta, chi) = self.fit_parameters(wf, version);
        let costheta = self.costheta_single_spin;
        let mut a2 = pnr_fits::alpha_a2(eta, chi, costheta);
        if a2 > 0.0 {
            a2 = 0.0;
        }
        let mut a3 = pnr_fits::alpha_a3(eta, chi, costheta).abs();
        if a3 < 0.00001 {
            a3 = 0.00001;
        }
        let max_a2 = PI * PI * a3.sqrt();
        if a2 < -max_a2 {
            a2 = -max_a2;
        }
        AlphaParams {
            a1: pnr_fits::alpha_a1(eta, chi, costheta).abs(),
            a2,
            a3,
            a4: pnr_fits::alpha_a4(eta, chi, costheta),
            ..Default::default()
        }
    }

    /// `IMRPhenomX_PNR_precompute_beta_coefficients`.
    fn precompute_beta_coefficients(&self, wf: &WaveformStruct, version: u32) -> BetaParams {
        let (eta, chi) = self.fit_parameters(wf, version);
        let costheta = self.costheta_single_spin;
        let costhetaf = self.costheta_final_single_spin;
        let mut b4 = pnr_fits::beta_b4(eta, chi, costheta);
        if b4 <= 175.0 {
            b4 = 175.0;
        }
        let b2 = pnr_fits::beta_b2(eta, chi, costheta);
        BetaParams {
            b0: costhetaf.acos() - pnr_fits::beta_b0(eta, chi, costheta),
            b1: pnr_fits::beta_b1(eta, chi, costheta),
            b2,
            b3: b2 * pnr_fits::beta_b3(eta, chi, costheta),
            b4,
            b5: pnr_fits::beta_b5(eta, chi, costheta),
            ..Default::default()
        }
    }

    // -----------------------------------------------------------------------
    // PN angles
    // -----------------------------------------------------------------------

    /// `alpha_SpinTaylor_IMR` (version 330).
    fn st_alpha(&self, mf: f64) -> Result<f64, IMRPhenomXError> {
        let st = self.st.as_ref().expect("SpinTaylor angles");
        if mf < st.ftrans_mrd {
            if !st.alpha_spline.contains(mf) {
                return Err(pnr_error(
                    "error in the interpolation of the SpinTaylor alpha",
                ));
            }
            Ok(st.alpha_spline.eval(mf) + self.alpha_offset)
        } else {
            Ok(st.alpha_mrd.alpha(mf) + self.alpha_offset)
        }
    }

    /// `beta_SpinTaylor_IMR` (version 330).
    fn st_beta(&self, mf: f64, bp: &BetaParams) -> Result<f64, IMRPhenomXError> {
        let st = self.st.as_ref().expect("SpinTaylor angles");
        let ftrans = if self.use_mr_beta {
            if st.ftrans_mrd <= 0.9 * bp.mf_beta_lower {
                st.ftrans_mrd - 0.0005
            } else {
                fma(0.9, bp.mf_beta_lower, -0.0005)
            }
        } else {
            st.ftrans_mrd
        };
        if mf < ftrans {
            if !st.cosbeta_spline.contains(mf) {
                return Err(pnr_error(
                    "error in the interpolation of the SpinTaylor beta",
                ));
            }
            let cosbeta = st.cosbeta_spline.eval(mf);
            let cosbeta = if 1.0 < cosbeta { 1.0 } else { cosbeta };
            let cosbeta = if -1.0 > cosbeta { -1.0 } else { cosbeta };
            Ok(cosbeta.acos())
        } else if self.use_mr_beta {
            Ok(beta_connection(mf, bp))
        } else {
            Ok(st.beta_mrd.beta(mf))
        }
    }

    /// `IMRPhenomX_PNR_GetPNAlphaAtFreq`: the PN `alpha` (NNLO, MSA or
    /// SpinTaylor).
    fn pn_alpha(&self, mf: f64) -> Result<f64, IMRPhenomXError> {
        let omega = PI * mf;
        let omega_cbrt = omega.cbrt();
        match self.version {
            102 => {
                let omega_cbrt2 = omega_cbrt * omega_cbrt;
                let logomega = omega.ln();
                let (alpha, _) = self.nnlo.evaluate(omega, logomega, omega_cbrt, omega_cbrt2);
                Ok(alpha - self.alpha_offset)
            }
            223 => {
                let msa = self.msa.as_ref().expect("MSA system");
                Ok(msa.phi_zeta_costhetal(omega_cbrt).phiz)
            }
            _ => {
                let alpha = self.st_alpha(mf)?;
                if !alpha.is_finite() {
                    return Err(pnr_error("alpha_SpinTaylor_IMR returned an invalid value"));
                }
                Ok(alpha)
            }
        }
    }

    /// `IMRPhenomX_PNR_GetPNBetaAtFreq`: the PN `beta`, with the two-spin
    /// oscillations tapered away below the lower connection frequency.
    fn pn_beta(
        &self,
        mf: f64,
        bp: &BetaParams,
        single: Option<&SingleSpin>,
    ) -> Result<f64, IMRPhenomXError> {
        let omega = PI * mf;
        let omega_cbrt = omega.cbrt();
        let beta_full = match self.version {
            102 => return Ok(self.nnlo_beta(omega_cbrt)),
            223 => {
                let msa = self.msa.as_ref().expect("MSA system");
                msa.phi_zeta_costhetal(omega_cbrt).cos_theta_l.acos()
            }
            _ => {
                let beta = self.st_beta(mf, bp)?;
                if !beta.is_finite() {
                    return Err(pnr_error("beta_SpinTaylor_IMR returned an invalid value"));
                }
                beta
            }
        };
        if self.two_spin() {
            let single = single.expect("single-spin system");
            let beta_single_spin = single.beta(omega_cbrt);
            if mf <= bp.mf_beta_lower {
                let oscillations = beta_full - beta_single_spin;
                let e = (2.0 * PI * mf / (4.0 * bp.mf_beta_lower)).cos();
                Ok(fma(oscillations, e * e, beta_single_spin))
            } else {
                Ok(beta_single_spin)
            }
        } else {
            Ok(beta_full)
        }
    }

    /// `IMRPhenomX_PNR_GetPNBetaAtFreq_fulltwospin`.
    fn pn_beta_fulltwospin(&self, mf: f64, bp: &BetaParams) -> Result<f64, IMRPhenomXError> {
        let omega = PI * mf;
        let omega_cbrt = omega.cbrt();
        match self.version {
            102 => Ok(self.nnlo_beta(omega_cbrt)),
            223 => {
                let msa = self.msa.as_ref().expect("MSA system");
                Ok(msa.phi_zeta_costhetal(omega_cbrt).cos_theta_l.acos())
            }
            _ => {
                let beta = self.st_beta(mf, bp)?;
                if !beta.is_finite() {
                    return Err(pnr_error("beta_SpinTaylor_IMR returned an invalid value"));
                }
                Ok(beta)
            }
        }
    }

    /// The NNLO opening angle at `v = (pi Mf)^(1/3)` (sign of `L + SL`
    /// ignored in `s = Sp / (L + SL)`).
    fn nnlo_beta(&self, v: f64) -> f64 {
        let l = self.l.ansatz(v, self.eta / v);
        let s = self.sperp / (l + self.sl);
        let s2 = s * s;
        (1.0f64.copysign(l + self.sl) / (1.0 + s2).sqrt()).acos()
    }

    /// `IMRPhenomX_PNR_PNWaveformBetaWrapper`: the "waveform" `beta` of the
    /// PN `beta` (Eq. 41 of arXiv:2107.08876).
    fn waveform_beta(&self, wf: &WaveformStruct, mf: f64, pn_beta: f64) -> f64 {
        let big_m = wf.m_tot;
        let m1 = wf.m1 * big_m;
        let m2 = wf.m2 * big_m;
        let chi_parr = (m1 + m2) * wf.chi_eff / m1;
        let v = (PI * mf).cbrt();
        let l_norm = wf.eta / v;
        let l_3pn = big_m * big_m * self.l.ansatz(v, l_norm);
        let chi_temp = chi_calc(m1, l_3pn, self.j0_norm, self.l_ref, chi_parr, pn_beta);
        let costheta_temp = chi_parr / chi_temp;
        pn_waveform_beta(mf, pn_beta, m1, m2, chi_temp, costheta_temp)
    }

    // -----------------------------------------------------------------------
    // connections of the PNR angles
    // -----------------------------------------------------------------------

    /// `IMRPhenomX_PNR_alpha_connection_parameters`.
    fn alpha_connection(
        &self,
        wf: &WaveformStruct,
        ap: &mut AlphaParams,
    ) -> Result<(), IMRPhenomXError> {
        let dmf = PNR_DMF;
        let mut mf_alpha_upper = ap.a4 / 3.0;
        let mut mf_alpha_lower = (3.0 / 3.5) * mf_alpha_upper;
        // just use the PN alpha below 2 Hz
        if mf_alpha_upper < hz_to_mf(2.0, wf.m_tot) {
            mf_alpha_lower = 100.0;
            mf_alpha_upper = 100.0;
        }
        // or below the start of the integration of the PN equations
        if self.version == 330 && mf_alpha_lower - 2. * dmf < self.mfmin_integration {
            mf_alpha_lower = 100.0;
            mf_alpha_upper = 100.0;
        }
        ap.mf_alpha_lower = mf_alpha_lower;
        ap.mf_alpha_upper = mf_alpha_upper;
        let a1 = self.pn_alpha(mf_alpha_lower - dmf)?;
        let alpha_lower = self.pn_alpha(mf_alpha_lower)?;
        let a3 = self.pn_alpha(mf_alpha_lower + dmf)?;
        let a4 = mr_alpha(mf_alpha_upper - dmf, ap);
        let alpha_upper = self.pn_alpha(mf_alpha_upper)?;
        let a6 = mr_alpha(mf_alpha_upper + dmf, ap);
        let derivative_alpha_lower = (a3 - a1) / (2.0 * dmf);
        let derivative_alpha_upper = (a6 - a4) / (2.0 * dmf);
        let args = (
            mf_alpha_lower,
            mf_alpha_upper,
            alpha_lower,
            alpha_upper,
            derivative_alpha_lower,
            derivative_alpha_upper,
        );
        let interp = [
            alpha_interpolation_0(args.0, args.1, args.2, args.3, args.4, args.5),
            alpha_interpolation_1(args.0, args.1, args.2, args.3, args.4, args.5),
            alpha_interpolation_2(args.0, args.1, args.2, args.3, args.4, args.5),
            alpha_interpolation_3(args.0, args.1, args.2, args.3, args.4, args.5),
        ];
        ap.alpha_interp = interp.map(|x| if x.is_nan() { 0.0 } else { x });
        let mr_alpha_at_mf_upper = mr_alpha(mf_alpha_upper, ap);
        ap.alpha_mr_offset = alpha_upper - mr_alpha_at_mf_upper;
        ap.alpha_lower = alpha_lower;
        ap.alpha_upper = alpha_upper;
        ap.derivative_alpha_lower = derivative_alpha_lower;
        ap.derivative_alpha_upper = derivative_alpha_upper;
        Ok(())
    }

    /// `IMRPhenomX_PNR_beta_connection_parameters` (versions 102, 223).
    fn beta_connection(
        &self,
        wf: &WaveformStruct,
        bp: &mut BetaParams,
        single: Option<&SingleSpin>,
    ) -> Result<(), IMRPhenomXError> {
        let dmf = PNR_DMF;
        beta_connection_frequencies(bp);
        let mf_beta_lower = bp.mf_beta_lower;
        let waveform_beta_at = |mf: f64| -> Result<f64, IMRPhenomXError> {
            let beta = self.pn_beta(mf, bp, single)?;
            Ok(self.waveform_beta(wf, mf, beta))
        };
        let b1 = waveform_beta_at(mf_beta_lower - dmf)?;
        let beta_lower = waveform_beta_at(mf_beta_lower)?;
        let b3 = waveform_beta_at(mf_beta_lower + dmf)?;
        let b4 = mr_beta(mf_beta_lower - dmf, bp);
        let beta_upper = mr_beta(mf_beta_lower, bp);
        let b6 = mr_beta(mf_beta_lower + dmf, bp);
        let derivative_beta_lower = (b3 - b1) / (2.0 * dmf);
        let derivative_beta_upper = (b6 - b4) / (2.0 * dmf);
        let r1 = beta_rescaling_1(
            mf_beta_lower,
            beta_lower,
            beta_upper,
            derivative_beta_lower,
            derivative_beta_upper,
        );
        let r2 = beta_rescaling_2(
            mf_beta_lower,
            beta_lower,
            beta_upper,
            derivative_beta_lower,
            derivative_beta_upper,
        );
        bp.beta_lower = beta_lower;
        bp.beta_upper = beta_upper;
        bp.derivative_beta_lower = derivative_beta_lower;
        bp.derivative_beta_upper = derivative_beta_upper;
        bp.beta_rescale_1 = if r1.is_nan() { 0.0 } else { r1 };
        bp.beta_rescale_2 = if r2.is_nan() { 0.0 } else { r2 };
        Ok(())
    }

    /// `IMRPhenomX_ST_PNR_beta_connection_parameters` (version 330): the
    /// cubic connecting the SpinTaylor `beta` to the (single-spin) MSA
    /// `beta`, and the rescaling of the inspiral `beta`.
    fn st_beta_connection(
        &self,
        wf: &WaveformStruct,
        bp: &mut BetaParams,
        single: Option<&SingleSpin>,
    ) -> Result<(), IMRPhenomXError> {
        let dmf = PNR_DMF;
        beta_connection_frequencies(bp);
        let st = self.st.as_ref().expect("SpinTaylor angles");
        // frequencies marking the beginning and end of the interpolation
        let ftrans = if st.ftrans_mrd <= 0.9 * bp.mf_beta_lower {
            st.ftrans_mrd
        } else {
            0.9 * bp.mf_beta_lower
        };
        let mfa = if ftrans - dmf < self.mfmin_integration + 2. * dmf {
            self.mfmin_integration + 2. * dmf
        } else {
            ftrans - dmf
        };
        let mfb = bp.mf_beta_lower;

        // SpinTaylor angles around the lower connection frequency
        let st_beta = |mf: f64| -> Result<f64, IMRPhenomXError> {
            if !st.cosbeta_spline.contains(mf) {
                return Err(pnr_error("error in beta interpolation"));
            }
            let cosbeta = st.cosbeta_spline.eval(mf);
            let cosbeta = if 1.0 < cosbeta { 1.0 } else { cosbeta };
            let cosbeta = if -1.0 > cosbeta { -1.0 } else { cosbeta };
            Ok(cosbeta.acos())
        };
        let ba1 = st_beta(mfa - dmf)?;
        let ba = st_beta(mfa)?;
        let ba2 = st_beta(mfa + dmf)?;

        // single-spin MSA angles at the upper connection frequency
        let omega = PI * mfb;
        let domega = PI * dmf;
        let main_msa = (!self.two_spin())
            .then(|| MsaSystem::new_unchecked(&self.msa_inputs(wf, self.chi1, self.chi2)));
        let msa_beta = |x: f64| match &main_msa {
            // evaluated with the precession struct of version 330
            Some(msa) => msa.cos_theta_l_with(x.cbrt(), &self.l).acos(),
            None => single.expect("single-spin system").beta(x.cbrt()),
        };
        let bb1 = msa_beta(omega - domega);
        let bb = msa_beta(omega);
        let bb2 = msa_beta(omega + domega);

        let dba = (ba2 - ba1) / (2.0 * dmf);
        let dbb = (bb2 - bb1) / (2.0 * dmf);
        // (sic: the SpinTaylor beta is passed for both ends)
        bp.beta_interp = [
            st_beta_coeffs_0(mfa, mfb, ba, ba, dba, dbb),
            st_beta_coeffs_1(mfa, mfb, ba, ba, dba, dbb),
            st_beta_coeffs_2(mfa, mfb, ba, ba, dba, dbb),
            st_beta_coeffs_3(mfa, mfb, ba, ba, dba, dbb),
        ];

        // rescaling
        let bbnew = if self.two_spin() { bb } else { ba };
        let bc1 = self.waveform_beta(wf, mfb - dmf, bb1);
        let bc = self.waveform_beta(wf, mfb, bbnew);
        let bc2 = self.waveform_beta(wf, mfb + dmf, bb2);
        let bd1 = mr_beta(mfb - dmf, bp);
        let bd = mr_beta(mfb, bp);
        let bd2 = mr_beta(mfb + dmf, bp);
        let dbc = (bc2 - bc1) / (2.0 * dmf);
        let dbd = (bd2 - bd1) / (2.0 * dmf);
        let r1 = beta_rescaling_1(mfb, bc, bd, dbc, dbd);
        let r2 = beta_rescaling_2(mfb, bc, bd, dbc, dbd);
        bp.beta_lower = bc;
        bp.beta_upper = bd;
        bp.derivative_beta_lower = dbc;
        bp.derivative_beta_upper = dbd;
        bp.beta_rescale_1 = if r1.is_nan() { 0.0 } else { r1 };
        bp.beta_rescale_2 = if r2.is_nan() { 0.0 } else { r2 };
        Ok(())
    }

    /// `IMRPhenomX_PNR_PopulateStructs`: the single-spin system of the
    /// two-spin map and the alpha and beta parameters.
    fn populate_structs(&self, wf: &WaveformStruct) -> Result<PnrStructs, IMRPhenomXError> {
        let single = if self.two_spin() {
            let chi = self.chi_single_spin;
            let costheta = self.costheta_single_spin;
            // IMRPhenomXGetAndSetPrecessionVariables of the single-spin
            // binary with PhenomXPrecVersion 223: only its MSA system is used
            let chi1 = [chi * costheta.acos().sin(), 0.0, chi * costheta];
            let msa = MsaSystem::new_unchecked(&self.msa_inputs(wf, chi1, [0.0, 0.0, 0.0]));
            let l_102 = msa
                .msa_error()
                .then(|| OrbitalAngularMomentum::three_pn(wf.eta, wf.delta, chi1[2], 0.0));
            Some(SingleSpin { msa, l_102 })
        } else {
            None
        };
        let mut alpha = self.precompute_alpha_coefficients(wf, self.version);
        self.alpha_connection(wf, &mut alpha)?;
        let mut beta = self.precompute_beta_coefficients(wf, self.version);
        if self.version == 330 {
            self.st_beta_connection(wf, &mut beta, single.as_ref())?;
        } else {
            self.beta_connection(wf, &mut beta, single.as_ref())?;
        }
        Ok(PnrStructs {
            single,
            alpha,
            beta,
        })
    }

    // -----------------------------------------------------------------------
    // PNR angles
    // -----------------------------------------------------------------------

    /// `IMRPhenomX_PNR_GeneratePNRAlphaAtMf`, Eq. 59 of arXiv:2107.08876.
    fn pnr_alpha_at(&self, mf: f64, ap: &AlphaParams) -> Result<f64, IMRPhenomXError> {
        if mf <= ap.mf_alpha_lower {
            return self.pn_alpha(mf);
        }
        if mf >= ap.mf_alpha_upper {
            return Ok(mr_alpha(mf, ap) + ap.alpha_mr_offset);
        }
        Ok(intermediate_alpha(mf, ap))
    }

    /// `IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf`: PNR and PN `alpha`
    /// blended in the transition region.
    fn merged_alpha_at(
        &self,
        wf: &WaveformStruct,
        mf: f64,
        ap: &AlphaParams,
    ) -> Result<f64, IMRPhenomXError> {
        let pnr_window = angles_window(wf.q, self.chi_single_spin);
        let msa_window = 1.0 - pnr_window;
        if mf <= ap.mf_alpha_lower {
            return self.pn_alpha(mf);
        }
        let alpha_pn = self.pn_alpha(mf)?;
        let alpha_pnr = if mf >= ap.mf_alpha_upper {
            mr_alpha(mf, ap) + ap.alpha_mr_offset
        } else {
            intermediate_alpha(mf, ap)
        };
        Ok(fma(pnr_window, alpha_pnr, msa_window * alpha_pn))
    }

    /// The `alpha` of the uniform grid: PNR, merged or PN depending on the
    /// region of parameter space.
    fn region_alpha_at(
        &self,
        wf: &WaveformStruct,
        mf: f64,
        ap: &AlphaParams,
    ) -> Result<f64, IMRPhenomXError> {
        let (q, chi) = (wf.q, self.chi_single_spin);
        if q <= PNR_Q_WINDOW_LOWER && chi <= PNR_CHI_WINDOW_LOWER {
            self.pnr_alpha_at(mf, ap)
        } else if q <= PNR_Q_WINDOW_UPPER && chi <= PNR_CHI_WINDOW_UPPER {
            self.merged_alpha_at(wf, mf, ap)
        } else {
            self.pn_alpha(mf)
        }
    }

    /// `IMRPhenomX_PNR_GeneratePNRBetaAtMf`, Eqs. 60 and 61 of
    /// arXiv:2107.08876.
    fn pnr_beta_at(
        &self,
        wf: &WaveformStruct,
        mf: f64,
        bp: &BetaParams,
        single: Option<&SingleSpin>,
    ) -> Result<f64, IMRPhenomXError> {
        if mf <= bp.mf_beta_lower {
            if self.version == 223 && bp.beta_lower < 0.01 * bp.beta_upper {
                // beta drops to almost zero before Mf_beta_lower
                return Ok(arctan_window(mr_beta(bp.mf_beta_lower, bp)));
            }
            let beta_pn = self.pn_beta(mf, bp, single)?;
            let beta_waveform = self.waveform_beta(wf, mf, beta_pn);
            return Ok(arctan_window(beta_waveform * rescale_beta(mf, bp)));
        }
        if mf >= bp.mf_beta_upper {
            return Ok(arctan_window(mr_beta(bp.mf_beta_upper, bp)));
        }
        Ok(arctan_window(mr_beta(mf, bp)))
    }

    /// `IMRPhenomX_PNR_GenerateMergedPNRBetaAtMf`: PNR and PN `beta` blended
    /// in the transition region (sets `UseMRbeta`).
    fn merged_beta_at(
        &mut self,
        wf: &WaveformStruct,
        mf: f64,
        bp: &BetaParams,
        single: Option<&SingleSpin>,
    ) -> Result<f64, IMRPhenomXError> {
        let pnr_window = angles_window(wf.q, self.chi_single_spin);
        let msa_window = 1.0 - pnr_window;
        if mf <= bp.mf_beta_lower {
            let (pnr_beta, pn_beta);
            if self.version != 330 {
                let beta_pn = self.pn_beta(mf, bp, single)?;
                let beta_waveform = self.waveform_beta(wf, mf, beta_pn);
                pnr_beta = beta_waveform * rescale_beta(mf, bp);
                pn_beta = beta_waveform;
            } else {
                self.use_mr_beta = true;
                let beta_pn1 = self.pn_beta(mf, bp, single)?;
                let beta_waveform1 = self.waveform_beta(wf, mf, beta_pn1);
                let rescaled_beta1 = beta_waveform1 * rescale_beta(mf, bp);
                self.use_mr_beta = false;
                let beta_pn2 = self.pn_beta(mf, bp, single)?;
                let beta_waveform2 = self.waveform_beta(wf, mf, beta_pn2);
                pnr_beta = rescaled_beta1;
                pn_beta = beta_waveform2;
            }
            return Ok(arctan_window(fma(
                pnr_window,
                pnr_beta,
                msa_window * pn_beta,
            )));
        }
        self.use_mr_beta = false;
        let beta_pn = self.pn_beta(mf, bp, single)?;
        let beta_waveform = self.waveform_beta(wf, mf, beta_pn);
        let mr = if mf >= bp.mf_beta_upper {
            mr_beta(bp.mf_beta_upper, bp)
        } else {
            mr_beta(mf, bp)
        };
        Ok(arctan_window(fma(
            pnr_window,
            mr,
            msa_window * beta_waveform,
        )))
    }

    /// `IMRPhenomX_PNR_GeneratePNRBetaNoMR`: the rescaled inspiral `beta`
    /// without the MR model.
    fn beta_no_mr(
        &self,
        wf: &WaveformStruct,
        mf: f64,
        bp: &BetaParams,
    ) -> Result<f64, IMRPhenomXError> {
        let pn_beta = self.pn_beta_fulltwospin(mf, bp)?;
        Ok(arctan_window(self.waveform_beta(wf, mf, pn_beta)))
    }

    /// `IMRPhenomX_PNR_GeneratePNRAngles_UniformFrequencies` without
    /// `gamma`: `alpha` and `beta` at `frequencies` (Hz).
    fn uniform_angles(
        &mut self,
        wf: &WaveformStruct,
        frequencies: &[f64],
        structs: &PnrStructs,
    ) -> Result<(Vec<f64>, Vec<f64>), IMRPhenomXError> {
        let (ap, bp) = (&structs.alpha, &structs.beta);
        let single = structs.single.as_ref();
        let (q, chi) = (wf.q, self.chi_single_spin);
        let n = frequencies.len();
        let mut alpha = Vec::with_capacity(n);
        let mut beta = Vec::with_capacity(n);
        let mf_of = |f: f64| hz_to_mf(f, wf.m_tot);
        if q <= PNR_Q_WINDOW_LOWER && chi <= PNR_CHI_WINDOW_LOWER {
            // inside the calibration region
            if attach_mr_beta(bp) {
                self.use_mr_beta = true;
                for &f in frequencies {
                    let mf = mf_of(f);
                    alpha.push(self.pnr_alpha_at(mf, ap)?);
                    beta.push(self.pnr_beta_at(wf, mf, bp, single)?);
                }
            } else {
                self.use_mr_beta = false;
                for &f in frequencies {
                    let mf = mf_of(f);
                    alpha.push(self.pnr_alpha_at(mf, ap)?);
                    beta.push(self.beta_no_mr(wf, mf, bp)?);
                }
            }
        } else if q <= PNR_Q_WINDOW_UPPER && chi <= PNR_CHI_WINDOW_UPPER {
            // inside the transition region
            if attach_mr_beta(bp) {
                for &f in frequencies {
                    let mf = mf_of(f);
                    alpha.push(self.merged_alpha_at(wf, mf, ap)?);
                    beta.push(self.merged_beta_at(wf, mf, bp, single)?);
                }
            } else {
                self.use_mr_beta = false;
                for &f in frequencies {
                    let mf = mf_of(f);
                    alpha.push(self.merged_alpha_at(wf, mf, ap)?);
                    beta.push(self.beta_no_mr(wf, mf, bp)?);
                }
            }
        } else {
            // outside the calibration region
            self.use_mr_beta = false;
            for &f in frequencies {
                let mf = mf_of(f);
                alpha.push(self.pn_alpha(mf)?);
                beta.push(self.beta_no_mr(wf, mf, bp)?);
            }
        }
        Ok((alpha, beta))
    }

    // -----------------------------------------------------------------------
    // SpinTaylor angles of version 330
    // -----------------------------------------------------------------------

    /// `IMRPhenomX_Initialize_Euler_Angles` for version 330 with
    /// multibanding off: the splines of the SpinTaylor `alpha` and
    /// `cos(beta)` and their merger-ringdown continuation. The spline of
    /// `gamma` is not needed with the PNR angles: only the conditions under
    /// which LAL fails to compute it are checked.
    pub(crate) fn initialize_euler_angles(
        &mut self,
        wf: &WaveformStruct,
    ) -> Result<(), IMRPhenomXError> {
        let buffer = self.integration_buffer;
        // start below fMin to avoid interpolation artefacts
        let fmin_angles = (wf.f_min - buffer) * 2. / self.m_max as f64;
        if fmin_angles.is_nan() || fmin_angles <= 0. {
            return Err(pnr_error(
                "fMin is too low and numerical angles could not be computed",
            ));
        }
        let mf_max = if wf.mf_max > fma(4., wf.f_damp, wf.f_ring) {
            wf.mf_max
        } else {
            fma(4., wf.f_damp, wf.f_ring)
        };
        let mf_max_angles = (mf_max + hz_to_mf(buffer, wf.m_tot)) * 2. / self.m_min as f64;
        let fmax_angles = mf_to_hz(mf_max_angles, wf.m_tot);

        // IMRPhenomX_SpinTaylorAnglesSplinesAll
        if fmax_angles <= 0.0 || fmax_angles <= fmin_angles || wf.f_ref < fmin_angles {
            return Err(pnr_error("invalid range of the SpinTaylor angles"));
        }
        let arrays = self
            .pn_arrays
            .take()
            .ok_or(pnr_error("no solution of the PN precession equations"))?;
        let st = self.interpolate_alpha_beta(wf, &arrays)?;

        // IMRPhenomX_InterpolateGamma_SpinTaylor: the failure conditions
        let mf_min = hz_to_mf(fmin_angles, wf.m_tot);
        let delta_mf = delta_mf_from_wfstruct(wf);
        let delta_f = mf_to_hz(delta_mf, wf.m_tot);
        let i_start = (fmin_angles / delta_f) as usize;
        let i_stop = (fmax_angles / delta_f) as usize + 1;
        let output_length = i_stop - i_start;
        if output_length < 4 {
            return Err(pnr_error(
                "number of points is insufficient for spline interpolation of gamma",
            ));
        }
        for i in 1..output_length {
            let mf = fma(i as f64, delta_mf, mf_min);
            if mf >= st.ftrans_mrd {
                continue;
            }
            // gamma_from_alpha_cosbeta evaluates the splines at five points
            let step = ((mf + delta_mf) - mf) * 0.25;
            for j in 0..5 {
                if !st.alpha_spline.contains(fma(j as f64, step, mf)) {
                    return Err(pnr_error(
                        "could not evaluate splines for alpha and/or gamma angles",
                    ));
                }
            }
        }
        let mf_last = fma((output_length - 1) as f64, delta_mf, mf_min);
        let gamma_contains = |mf: f64| mf >= mf_min && mf <= mf_last;
        if !gamma_contains(wf.mf_ref) {
            return Err(pnr_error(
                "gamma could not be evaluated at the reference frequency",
            ));
        }
        // alpha(f_ref)
        let mf_alpha_ref = if wf.mf_ref <= st.ftrans_mrd {
            wf.mf_ref
        } else {
            st.ftrans_mrd
        };
        if !st.alpha_spline.contains(mf_alpha_ref) {
            return Err(pnr_error("could not evaluate alpha(f_ref)"));
        }
        // the angles at the end of the inspiral
        if !(st.alpha_spline.contains(st.ftrans_mrd) && gamma_contains(st.ftrans_mrd)) {
            return Err(pnr_error(
                "could not compute the angles at the end of the inspiral",
            ));
        }
        self.st = Some(st);
        Ok(())
    }

    /// `IMRPhenomX_InterpolateAlphaBeta_SpinTaylor` for version 330 with the
    /// PNR angles and `PhenomXPFinalSpinMod` 7.
    fn interpolate_alpha_beta(
        &self,
        wf: &WaveformStruct,
        arrays: &PnArrays,
    ) -> Result<St330, IMRPhenomXError> {
        let len_pn = arrays.len();
        let mut fgw = Vec::with_capacity(len_pn);
        let mut alphaaux = Vec::with_capacity(len_pn);
        let mut cosbeta = Vec::with_capacity(len_pn);
        let mut mfmax_pn = 0.;
        for i in 0..len_pn {
            let lnhat = [arrays.lnhat[0][i], arrays.lnhat[1][i], arrays.lnhat[2][i]];
            let lnhat = rotate_z(
                -self.kappa,
                rotate_y(-self.theta_j_sf, rotate_z(-self.phi_j_sf, lnhat)),
            );
            let fgw_hz = arrays.v[i].powf(3.) / self.pi_gm;
            let fgw_mf = hz_to_mf(fgw_hz, wf.m_tot);
            if fgw_hz.is_nan() || fgw_hz <= 0. {
                return Err(pnr_error("non-positive frequency in the PN solution"));
            }
            alphaaux.push(lnhat[1].atan2(lnhat[0]));
            cosbeta.push(lnhat[2]);
            fgw.push(fgw_mf);
            mfmax_pn = fgw_mf;
        }
        if fgw
            .windows(2)
            .any(|w| w[0] >= w[1] || w[0].is_nan() || w[1].is_nan())
        {
            return Err(pnr_error("PN frequencies are not strictly increasing"));
        }
        let mut fmax_inspiral = if self.use_tuned_angles {
            mfmax_pn
        } else {
            mfmax_pn - wf.delta_mf
        };
        if fmax_inspiral > wf.f_ring - wf.f_damp {
            fmax_inspiral = 1.020 * wf.f_meco;
        }
        let ftrans_mrd = 0.98 * fmax_inspiral;
        let mut alpha = alphaaux;
        unwrap_in_place(&mut alpha);
        let alpha_spline = CubicSpline::new(&fgw, &alpha);
        let cosbeta_spline = CubicSpline::new(&fgw, &cosbeta);
        let alpha_mrd = AlphaMrd::new(&alpha_spline, fmax_inspiral)?;
        let beta_mrd = BetaMrd330::new(
            &cosbeta_spline,
            fmax_inspiral,
            wf,
            self.ringdown_pnr_beta(wf),
        );
        Ok(St330 {
            alpha_spline,
            cosbeta_spline,
            ftrans_mrd,
            alpha_mrd,
            beta_mrd,
        })
    }

    // -----------------------------------------------------------------------
    // interpolants of the PNR angles
    // -----------------------------------------------------------------------

    /// `IMRPhenomX_PNR_GeneratePNRAngleInterpolants` for the modes `modes`
    /// (`(l, m' > 0)`) with multibanding off: the PNR `alpha`, `beta` on a
    /// uniform grid covering the frequencies at which the modes use the
    /// (2,2) angles, and `gamma` integrated from them.
    pub(crate) fn angle_interpolants(
        &mut self,
        wf: &WaveformStruct,
        modes: &[(u32, u32)],
    ) -> Result<PnrSplines, IMRPhenomXError> {
        let mut structs = self.populate_structs(wf)?;

        // connection frequencies for the frequency map of the higher modes
        self.pnr_hm_mfhigh = structs.beta.mf_beta_lower;
        self.pnr_hm_mflow = structs.alpha.mf_alpha_lower;
        if self.pnr_hm_mfhigh > wf.f_cut_def || self.pnr_hm_mfhigh < 0.1 * wf.f_ring {
            self.pnr_hm_mfhigh = wf.f_ring;
        }
        if self.pnr_hm_mflow > wf.f_cut_def || self.pnr_hm_mfhigh < self.pnr_hm_mflow {
            self.pnr_hm_mflow = self.pnr_hm_mfhigh / 2.0;
        }

        let mut f_min_22 = wf.f_min;
        let mut f_max_22 = wf.f_max_prime;
        if wf.delta_f != 0. {
            // account for the finite frequency spacing in fMin
            let i_start = (f_min_22 / wf.delta_f) as usize;
            let f_min_from_delta_f = i_start as f64 * wf.delta_f;
            if f_min_from_delta_f < f_min_22 {
                f_min_22 = f_min_from_delta_f;
            }
        } else {
            f_max_22 = wf.f_max;
        }

        // the range of the mapped frequencies of the higher modes
        let mut f_min = f_min_22;
        let mut f_max = f_max_22;
        for ell in 2..=4u32 {
            for emm in 1..=ell {
                if (ell, emm) == (2, 2) || !modes.contains(&(ell, emm)) {
                    continue;
                }
                if wf.q == 1.0 && wf.chi1l == wf.chi2l && emm % 2 != 0 {
                    continue;
                }
                let mf_rd_22 = wf.f_ring;
                let mf_rd_lm = ringdown_frequency(ell, emm, wf);
                let (mf_low, mf_high) = linear_frequency_map_transition_frequencies(
                    emm,
                    mf_rd_22,
                    mf_rd_lm,
                    self.pnr_hm_mflow,
                    self.pnr_hm_mfhigh,
                );
                let map = |f: f64| {
                    linear_frequency_map(
                        hz_to_mf(f, wf.m_tot),
                        ell,
                        emm,
                        mf_low,
                        mf_high,
                        mf_rd_22,
                        mf_rd_lm,
                        self.pnr_inspiral_scaling,
                    )
                    .map(|mf| mf_to_hz(mf, wf.m_tot))
                    .ok_or(pnr_error(
                        "invalid transition frequencies of the frequency map",
                    ))
                };
                let fmin_mapped = map(f_min_22)?;
                let fmax_mapped = map(f_max_22)?;
                if fmin_mapped < f_min {
                    f_min = fmin_mapped;
                }
                if fmax_mapped > f_max {
                    f_max = fmax_mapped;
                }
            }
        }

        let delta_f_int = self.interpolation_delta_f(wf, f_min);
        // guarantee space to avoid extrapolation and boundary effects
        f_min = if fnms(2.0, delta_f_int, f_min) < 0. {
            f_min / 2.0
        } else {
            fnms(2.0, delta_f_int, f_min)
        };
        f_max = fma(2.0, delta_f_int, f_max);
        let i_start = (f_min / delta_f_int) as usize;
        let i_stop = (f_max / delta_f_int) as usize + 1;
        if i_stop < i_start + 3 {
            return Err(pnr_error("too few points for the PNR angle interpolants"));
        }
        let freqs_22: Vec<f64> = (i_start..i_stop).map(|i| i as f64 * delta_f_int).collect();

        // the SpinTaylor angles must start below the interpolation grid
        if self.version == 330 && self.fmin_integration > freqs_22[0] {
            self.integrate_spintaylor(wf, freqs_22[0]);
            if self.version != 330 {
                // LAL then uses freed memory
                return Err(pnr_error(
                    "the SpinTaylor angles failed when extending them to lower frequencies",
                ));
            }
            self.initialize_euler_angles(wf)?;
            structs = self.populate_structs(wf)?;
        }

        let (alpha, beta) = self.uniform_angles(wf, &freqs_22, &structs)?;
        let alpha_spline = CubicSpline::new(&freqs_22, &alpha);
        let beta_spline = CubicSpline::new(&freqs_22, &beta);
        // gamma from the minimal rotation condition with Boole's rule
        let mut gamma = Vec::with_capacity(freqs_22.len());
        gamma.push(0.0);
        // (the interval search starts from the previous one, as with GSL's
        // accelerators; the interval found does not depend on it)
        let hint = std::cell::Cell::new(0);
        let integrand = |t: f64| {
            let index = alpha_spline.locate(t, hint.get());
            hint.set(index);
            alphadot_cosbeta(
                alpha_spline.derivative_in(index, t),
                beta_spline.eval_in(index, t),
            )
        };
        for i in 1..freqs_22.len() {
            let previous = gamma[i - 1];
            gamma.push(boole_step(
                previous,
                freqs_22[i - 1],
                freqs_22[i],
                integrand,
            ));
        }
        let gamma_spline = CubicSpline::new(&freqs_22, &gamma);
        Ok(PnrSplines {
            alpha: alpha_spline,
            beta: beta_spline,
            gamma: gamma_spline,
        })
    }

    /// `IMRPhenomX_PNR_RemapThetaJSF`: recompute the J frame of the line of
    /// sight with the PNR `beta` at the reference frequency in place of the
    /// angle between J and L (updates `alpha_offset`, `thetaJN`, `zeta` and
    /// the spherical harmonics).
    pub(crate) fn remap_theta_jsf(&mut self, wf: &WaveformStruct, beta_ref: f64) {
        if self.j0_norm < 1e-10 {
            return;
        }
        if chi_in_plane(self.chi1, self.chi2) < 1e-3 {
            self.alpha_offset -= self.alpha0;
            return;
        }
        self.theta_j_sf = beta_ref;
        let n_sf = Self::line_of_sight(wf);
        let tmp = rotate_y(-self.theta_j_sf, rotate_z(-self.phi_j_sf, n_sf));
        self.kappa = atan2_tol(tmp[1], tmp[0], MAX_TOL_ATAN);
        // convention 1
        self.alpha0 = PI - self.kappa;
        self.alpha_offset -= self.alpha0;
        let n_jf = rotate_z(-self.kappa, tmp);
        self.theta_jn = n_jf[2].acos();
        self.zeta_polarization = self.polarization_angle(wf, n_jf[0], n_jf[2]);
        self.ylm = SphericalHarmonics::new(self.theta_jn);
    }

    /// The NNLO Euler angles `(alpha, cos(beta), epsilon)` at
    /// `omega = 2 pi Mf / m'` for the twisting up without the PNR angles
    /// (`Get_alpha_beta_epsilon`, including the reference offsets).
    pub(crate) fn nnlo_angles(&self, omega: f64) -> (f64, f64, f64) {
        let logomega = omega.ln();
        let omega_cbrt = omega.cbrt();
        let omega_cbrt2 = omega_cbrt * omega_cbrt;
        let (alpha, epsilon) = self.nnlo.evaluate(omega, logomega, omega_cbrt, omega_cbrt2);
        let v = omega_cbrt;
        let l = self.l.ansatz_3pn(v, self.eta / v);
        // the sign of L + SL is ignored in s := Sp / (L + SL)
        let s = self.sperp / (l + self.sl);
        let s2 = s * s;
        let cos_beta = 1.0f64.copysign(l + self.sl) / (1.0 + s2).sqrt();
        (
            alpha - self.alpha_offset,
            cos_beta,
            epsilon - self.epsilon_offset,
        )
    }

    // -----------------------------------------------------------------------
    // antisymmetric waveform
    // -----------------------------------------------------------------------

    /// The function `kappa(v)` of `IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio`
    /// before the smoothing: the ratio of the antisymmetric to the symmetric
    /// (2,2) amplitude as a function of `v = (pi Mf)^(1/3)` below the
    /// ringdown frequency, and its value at the ringdown frequency.
    fn antisymmetric_ratio(&self, wf: &WaveformStruct) -> (impl Fn(f64) -> f64, f64) {
        let m2 = wf.m2;
        let delta = fnms(2.0, m2, 1.0);
        let eta = wf.eta;
        let theta = self.theta_antisymmetric;
        let chi = self.chi_single_spin_antisymmetric;
        // the phenomenological fit of the amplitude ratio in the PN-NR
        // transition
        let (b0, b1, b2, b3) = (18.0387, 15.4509, 55.1140, -203.6290);
        let b = fma(b3 * eta, theta, fma(b2, theta, fma(b1, eta, b0)));
        let (sin_theta, cos_theta) = (theta.sin(), theta.cos());
        let kappa_at = move |v: f64| {
            let v2 = v * v;
            let v3 = v2 * v;
            let v5 = v3 * v2;
            let num = 21.0 * v2 * (1.0 + delta) * chi * sin_theta;
            let den = 2.0
                * fma(
                    -(28.0 * v3 * (1.0 + delta - eta) * chi),
                    cos_theta,
                    fma(v2, fma(55.0, eta, -107.0), fma(84.0 * PI, v3, 42.0)),
                );
            num / den * fma(b, v5, 1.0)
        };
        let kappa_rd = kappa_at((PI * wf.f_ring).cbrt());
        (kappa_at, kappa_rd)
    }

    /// `IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio`: the ratio of the
    /// antisymmetric to the symmetric (2,2) amplitude, smoothed over 80
    /// frequencies.
    fn antisymmetric_amp_ratio(&self, wf: &WaveformStruct, frequencies: &[f64]) -> Vec<f64> {
        let (kappa_at, kappa_rd) = self.antisymmetric_ratio(wf);
        let mut kappa: Vec<f64> = frequencies
            .iter()
            .map(|&f| {
                let mf = hz_to_mf(f, wf.m_tot);
                if mf < wf.f_ring {
                    kappa_at((PI * mf).cbrt())
                } else {
                    kappa_rd
                }
            })
            .collect();
        smooth_amp_ratio(&mut kappa, frequencies);
        kappa
    }

    /// [`PnrPrecession::antisymmetric_amp_ratio`] on a grid, with the powers
    /// of the frequencies from the grid (equal to rounding).
    fn antisymmetric_amp_ratio_on_grid(
        &self,
        wf: &WaveformStruct,
        grid: &FrequencyGrid,
    ) -> Vec<f64> {
        let (kappa_at, kappa_rd) = self.antisymmetric_ratio(wf);
        // v = (pi Mf)^(1/3) = (pi M)^(1/3) f^(1/3)
        let scale = ScalePowers::new(PI * (MTSUN_SI * wf.m_tot));
        let frequencies = grid.frequencies();
        let mut kappa: Vec<f64> = (0..frequencies.len())
            .map(|i| {
                if hz_to_mf(frequencies[i], wf.m_tot) < wf.f_ring {
                    kappa_at(grid.powers(&scale, i).cbrt)
                } else {
                    kappa_rd
                }
            })
            .collect();
        smooth_amp_ratio(&mut kappa, frequencies);
        kappa
    }

    /// `IMRPhenomX_PNR_GenerateAntisymmetricPhaseCoefficients` and the
    /// other frequency-independent quantities of
    /// `IMRPhenomX_PNR_GenerateAntisymmetricWaveform`.
    ///
    /// Below `Mf_T = 0.85 fRING` the phase is `phi_22 / 2 + alpha + A0 Mf +
    /// phi_A0`, above it follows `phi_22`. In LALSimulation 6.2.0 `phi_A0`
    /// is `alpha_offset = alpha(f_ref) - alpha0` and the polarization angle
    /// `zeta` is added; both depend on the direction of the line of sight
    /// (`alpha0 = pi - kappa`), and neither changes sign under the mirror
    /// reflection of the binary. With the fix of lalsuite !2593
    /// (`parity_fix`) the phase is referred to the frame of the binary:
    /// `phi_A0 = -alpha(f_ref) - pi / 2` and the azimuth `phiJ_Sf` of `J` is
    /// added, which moves by `pi` (a sign flip of the antisymmetric mode)
    /// when the in-plane spins are reversed.
    fn antisymmetric_coefficients(
        &self,
        wf: &WaveformStruct,
    ) -> Result<AntisymmetricCoefficients, IMRPhenomXError> {
        let mut phase22 = PhaseCoefficients::new(wf);
        phase_22_connection_coefficients(&mut phase22);
        let amp22 = AmpCoefficients::new(wf);
        let mf_t = 0.85 * wf.f_ring;
        let lina = 0.0;
        let linb = time_shift_22(&phase22, wf);
        let inveta = 1.0 / wf.eta;
        let phiref22 = {
            let powers = UsefulPowers::new(wf.mf_ref);
            fma(
                2.0,
                wf.phi0,
                fma(
                    -inveta,
                    phase_22(wf.mf_ref, &powers, &phase22),
                    -(linb * wf.mf_ref),
                ) - lina,
            ) + FRAC_PI_4
        };

        // IMRPhenomX_PNR_GenerateAntisymmetricPhaseCoefficients
        let mut alpha = self.precompute_alpha_coefficients(wf, self.version);
        self.alpha_connection(wf, &mut alpha)?;
        let delta_mf = 0.0005;
        let mf_right = mf_t + delta_mf;
        let mf_left = mf_t - delta_mf;
        let alpha_der_mf_t = (self.region_alpha_at(wf, mf_right, &alpha)?
            - self.region_alpha_at(wf, mf_left, &alpha)?)
            / 2.0
            / delta_mf;
        let alpha_mf_t = self.region_alpha_at(wf, mf_t, &alpha)?;
        let powers = UsefulPowers::new(mf_t);
        let (dphi, phi) = if mf_t < phase22.f_phase_match_in {
            (
                inspiral::phase_22_ansatz_exact(mf_t, &powers, &phase22),
                inspiral::phase_22_ansatz_int_exact(mf_t, &powers, &phase22),
            )
        } else if mf_t > phase22.f_phase_match_im {
            (
                ringdown::phase_22_ansatz_exact(mf_t, &powers, &phase22) + phase22.c2_mrd,
                fma(
                    phase22.c2_mrd,
                    mf_t,
                    ringdown::phase_22_ansatz_int_exact(mf_t, &powers, &phase22) + phase22.c1_mrd,
                ),
            )
        } else {
            (
                intermediate::phase_22_ansatz_exact(mf_t, &powers, &phase22) + phase22.c2_int,
                fma(
                    phase22.c2_int,
                    mf_t,
                    intermediate::phase_22_ansatz_int_exact(mf_t, &powers, &phase22)
                        + phase22.c1_int,
                ),
            )
        };
        let dphi = dphi * inveta + linb;
        let phi = phi * inveta + (fma(linb, mf_t, lina) + phiref22);
        let a0 = dphi / 2.0 - alpha_der_mf_t;
        let phi_a0 = if self.parity_fix {
            // alpha at the reference frequency (in the region of alpha_MfT)
            // instead of the offset, which includes the azimuth of the line
            // of sight
            let alpha_ref = self.region_alpha_at(wf, wf.mf_ref, &alpha)?;
            -alpha_ref - FRAC_PI_2
        } else {
            self.alpha_offset
        };
        let phi_b0 = fma(a0, mf_t, alpha_mf_t - phi / 2.0) + phi_a0;
        // the azimuth of J in the L frame instead of the polarization angle
        let phase_offset = if self.parity_fix {
            self.phi_j_sf
        } else {
            self.zeta_polarization
        };
        Ok(AntisymmetricCoefficients {
            phase22,
            amp22,
            alpha,
            mf_t,
            lina,
            linb,
            inveta,
            phiref22,
            a0,
            phi_a0,
            phi_b0,
            phase_offset,
            amp_prefactor: wf.amp0 * wf.amp_norm,
        })
    }

    /// `IMRPhenomX_PNR_GenerateAntisymmetricWaveform`: amplitude and phase
    /// of the antisymmetric (2,2) co-precessing mode at `frequencies` (Hz),
    /// `h = amp exp(i phase)` (see
    /// [`PnrPrecession::antisymmetric_coefficients`]).
    pub(crate) fn antisymmetric_waveform(
        &self,
        wf: &WaveformStruct,
        frequencies: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), IMRPhenomXError> {
        let kappa = self.antisymmetric_amp_ratio(wf, frequencies);
        let coefficients = self.antisymmetric_coefficients(wf)?;
        // IMRPhenomX_PNR_GeneratePNRAlphaForAntisymmetry
        let alpha_pnr = frequencies
            .iter()
            .map(|&f| self.region_alpha_at(wf, hz_to_mf(f, wf.m_tot), &coefficients.alpha))
            .collect::<Result<Vec<f64>, IMRPhenomXError>>()?;
        let mut amp = Vec::with_capacity(frequencies.len());
        let mut phase = Vec::with_capacity(frequencies.len());
        for (i, &f) in frequencies.iter().enumerate() {
            let mf = wf.m_sec * f;
            let (a, p) =
                coefficients.amplitude_phase(mf, &UsefulPowers::new(mf), kappa[i], alpha_pnr[i]);
            amp.push(a);
            phase.push(p);
        }
        Ok((amp, phase))
    }

    /// [`PnrPrecession::antisymmetric_waveform`] at the first `n` points of
    /// `grid`, with the powers of the frequencies from the grid (equal to
    /// rounding). The amplitude ratio is smoothed over the whole grid, as
    /// `antisymmetric_waveform` smooths it over all its frequencies; the PNR
    /// `alpha` is evaluated only below `Mf_T`, where it is used.
    pub(crate) fn antisymmetric_waveform_on_grid(
        &self,
        wf: &WaveformStruct,
        grid: &FrequencyGrid,
        n: usize,
    ) -> Result<(Vec<f64>, Vec<f64>), IMRPhenomXError> {
        let kappa = self.antisymmetric_amp_ratio_on_grid(wf, grid);
        let coefficients = self.antisymmetric_coefficients(wf)?;
        let frequencies = &grid.frequencies()[..n];
        let n_t = frequencies.partition_point(|&f| wf.m_sec * f < coefficients.mf_t);
        let alpha_pnr = frequencies[..n_t]
            .iter()
            .map(|&f| self.region_alpha_at(wf, hz_to_mf(f, wf.m_tot), &coefficients.alpha))
            .collect::<Result<Vec<f64>, IMRPhenomXError>>()?;
        let scale = ScalePowers::new(wf.m_sec);
        let mut amp = Vec::with_capacity(n);
        let mut phase = Vec::with_capacity(n);
        for (i, &kappa) in kappa[..n].iter().enumerate() {
            let powers = grid.powers(&scale, i);
            let alpha = if i < n_t { alpha_pnr[i] } else { 0.0 };
            let (a, p) = coefficients.amplitude_phase(
                powers.itself,
                &UsefulPowers::from_grid(&powers),
                kappa,
                alpha,
            );
            amp.push(a);
            phase.push(p);
        }
        Ok((amp, phase))
    }
}

/// The frequency-independent quantities of the antisymmetric (2,2) mode
/// ([`PnrPrecession::antisymmetric_coefficients`]).
struct AntisymmetricCoefficients {
    phase22: PhaseCoefficients,
    amp22: AmpCoefficients,
    /// The PNR `alpha` parameters (`IMRPhenomX_PNR_alpha_parameters`).
    alpha: AlphaParams,
    mf_t: f64,
    lina: f64,
    linb: f64,
    inveta: f64,
    phiref22: f64,
    a0: f64,
    phi_a0: f64,
    phi_b0: f64,
    /// `zeta_polarization` (6.2.0) or `phiJ_Sf` (with the fix).
    phase_offset: f64,
    amp_prefactor: f64,
}

impl AntisymmetricCoefficients {
    /// Amplitude and phase of the antisymmetric mode at `mf` with the
    /// powers `powers`, the smoothed amplitude ratio `kappa` and the PNR
    /// `alpha` (used below `Mf_T` only).
    #[inline]
    fn amplitude_phase(
        &self,
        mf: f64,
        powers: &UsefulPowers,
        kappa: f64,
        alpha: f64,
    ) -> (f64, f64) {
        // symmetric phase and amplitude
        let phi = self.phase22.phase_xas(mf, powers) * self.inveta;
        let phi = phi + (fma(self.linb, mf, self.lina) + self.phiref22);
        let a = self.amp22.ansatz(mf, powers);
        // antisymmetric amplitude and phase
        let amp_as = kappa * a;
        let phi_as = if mf < self.mf_t {
            fma(self.a0, mf, phi / 2.0 + alpha) + self.phi_a0
        } else {
            phi + self.phi_b0
        };
        (
            self.amp_prefactor * powers.m_seven_sixths * amp_as,
            self.phase_offset + phi_as,
        )
    }
}

/// The moving average of `IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio`,
/// in place as in LAL: `kappa[id + 40] = sum_j kappa[id + j] df_j / sum_j
/// df_j` over the 81 points `j = 0..=80` (fewer on short grids) with `df_j =
/// f[id + j + 1] - f[id + j]`, for `id = 0, 1, ...` in turn (so the later
/// windows include values already smoothed).
///
/// The operations are those of LAL, except that a window whose inputs (the
/// values of `kappa` it reads and the steps) are bit for bit those of the
/// previous window takes its result instead of recomputing it: the result
/// is the same. Above the ringdown frequency, where `kappa` is constant,
/// the smoothed values settle on a constant, so on a uniform grid most
/// windows there are not recomputed.
fn smooth_amp_ratio(kappa: &mut [f64], frequencies: &[f64]) {
    let len = kappa.len();
    let mut width = 80usize;
    if width > len.saturating_sub(1) {
        width = (len as f64 / 2.0).floor() as usize;
    }
    let half_width = (width as f64 / 2.0).floor() as usize;
    let end = len.saturating_sub(width + 1);
    if end == 0 {
        return;
    }
    let df: Vec<f64> = frequencies.windows(2).map(|w| w[1] - w[0]).collect();
    // the first index of the run of bit-for-bit equal values ending at each
    // index (of the input kappa and of the steps)
    let run_starts = |values: &[f64]| {
        let mut starts = vec![0usize; values.len()];
        for i in 1..values.len() {
            starts[i] = if values[i].to_bits() == values[i - 1].to_bits() {
                starts[i - 1]
            } else {
                i
            };
        }
        starts
    };
    let kappa_start = run_starts(kappa);
    let df_start = run_starts(&df);
    let mut previous = 0.0f64;
    // the number of equal results ending with the previous one
    let mut equal_results = 0usize;
    for id in 0..end {
        // Window id reads the results of windows id - hw ..= id - 1 at
        // j < hw (hw = half_width) and the input values at j >= hw. It
        // repeats window id - 1 if the results of windows id - hw - 1 ..= id
        // - 1 are equal, the input values are constant over id - 1 + hw ..=
        // id + width and the steps over id - 1 ..= id + width.
        let repeat = id > half_width
            && equal_results > half_width
            && kappa_start[id + width] < id + half_width
            && df_start[id + width] < id;
        let value = if repeat {
            previous
        } else {
            let mut smoothed_ratio = 0.0;
            let mut frequency_width = 0.0;
            for j in 0..=width {
                let step = df[id + j];
                smoothed_ratio = fma(kappa[id + j], step, smoothed_ratio);
                frequency_width += step;
            }
            smoothed_ratio * 1.0 / frequency_width
        };
        equal_results = if id > 0 && value.to_bits() == previous.to_bits() {
            equal_results + 1
        } else {
            1
        };
        previous = value;
        kappa[id + half_width] = value;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The moving average as LAL writes it.
    fn smooth_as_lal(kappa: &mut [f64], frequencies: &[f64]) {
        let len = kappa.len();
        let mut width = 80usize;
        if width > len.saturating_sub(1) {
            width = (len as f64 / 2.0).floor() as usize;
        }
        let half_width = (width as f64 / 2.0).floor() as usize;
        for id in 0..len.saturating_sub(width + 1) {
            let mut smoothed_ratio = 0.0;
            let mut frequency_width = 0.0;
            for j in 0..=width {
                let df = frequencies[id + j + 1] - frequencies[id + j];
                smoothed_ratio = fma(kappa[id + j], df, smoothed_ratio);
                frequency_width += df;
            }
            kappa[id + half_width] = smoothed_ratio * 1.0 / frequency_width;
        }
    }

    /// The windows skipped by `smooth_amp_ratio` do not change the result:
    /// uniform and geometric grids, a ratio constant above a "ringdown"
    /// frequency, constant everywhere, with plateaus, and short grids.
    #[test]
    fn smoothing_is_that_of_lal() {
        let uniform: Vec<f64> = (0..3000).map(|i| 20.0 + 0.125 * i as f64).collect();
        let geometric: Vec<f64> = (0..3000).map(|i| 20.0 * 1.001f64.powi(i)).collect();
        let short: Vec<f64> = (0..60).map(|i| 20.0 + 0.5 * i as f64).collect();
        let ratios: [&dyn Fn(f64) -> f64; 4] = [
            &|f| {
                if f < 150.0 {
                    0.1 / (1.0 + f / 300.0)
                } else {
                    0.1 / 1.5
                }
            },
            &|_| 0.07,
            &|f| ((f / 7.0).floor() * 0.01).min(0.2),
            &|f| 0.1 * (f / 50.0).sin(),
        ];
        for frequencies in [&uniform, &geometric, &short] {
            for ratio in ratios {
                let kappa: Vec<f64> = frequencies.iter().map(|&f| ratio(f)).collect();
                let (mut fast, mut lal) = (kappa.clone(), kappa);
                smooth_amp_ratio(&mut fast, frequencies);
                smooth_as_lal(&mut lal, frequencies);
                assert!(fast
                    .iter()
                    .zip(&lal)
                    .all(|(a, b)| a.to_bits() == b.to_bits()));
            }
        }
    }
}
