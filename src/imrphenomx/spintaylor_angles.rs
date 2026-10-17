//! Numerical "SpinTaylor" Euler angles of IMRPhenomXP/XPHM
//! (`PhenomXPrecVersion` 320, arXiv:2412.16721), following
//! `LALSimIMRPhenomX_precession.c`:
//!
//! - [`inspiral_angles_spintaylor`] ports `IMRPhenomX_InspiralAngles_SpinTaylor`,
//!   the integration of the orbit-averaged SpinTaylorT4 equations (coarse
//!   backward/forward runs around `f_ref`, then a fine run at high
//!   frequency). This is the only call into [`crate::spintaylor`].
//! - [`SpinTaylorAngles::new`] ports `IMRPhenomX_InterpolateAlphaBeta_SpinTaylor`:
//!   splines of `alpha` and `cos(beta)` against `Mf` in the J frame, the
//!   precessing final spin from the spins at the end of the integration
//!   (`PhenomXPFinalSpinMod` 4) and the coefficients of the merger-ringdown
//!   continuation (`alphaMRD_coeff`, `betaMRD_coeff`).
//! - [`SpinTaylorAngles::interpolate_gamma`] ports
//!   `IMRPhenomX_InterpolateGamma_SpinTaylor`: `gamma` from the minimal
//!   rotation condition integrated with Boole's rule
//!   (`gamma_from_alpha_cosbeta`) on a uniform grid.
//! - [`SpinTaylorAngles::angles_imr`] ports `IMRPhenomXPSpinTaylorAnglesIMR`,
//!   the angles on the frequency grid of the waveform.
//!
//! Like the SpinTaylor integrator, the arithmetic mirrors the fused
//! multiply-adds of the arm64 builds of LALSimulation and GSL (see
//! [`crate::spintaylor::fp`]): on that platform all intermediate quantities
//! (PN solution, splines, final spin, merger-ringdown coefficients, `gamma`)
//! agree with LALSimulation bit for bit.

use std::f64::consts::PI;

use crate::constants::{C_SI, G_SI, REAL4_EPS};
use crate::imrphenomx::internals::{IMRPhenomXError, WaveformStruct};
use crate::imrphenomx::precession::{
    delta_mf_from_wfstruct, hz_to_mf, mf_to_hz, precessing_final_spin_2017, rotate_y_cs,
    rotate_z_cs, PrecessionStruct,
};
use crate::interpolate::CubicSpline;
use crate::spintaylor::fp::{fma, fms, fnms};
use crate::spintaylor::{SpinTaylorOrbit, SpinTaylorT4};

/// Knots of the `gamma` grid integrated beyond the last frequency at which
/// the spline is evaluated, see [`SpinTaylorAngles::interpolate_gamma`].
const GAMMA_SPLINE_MARGIN: usize = 64;

/// The solution of the PN precession equations used for the angles
/// (`PhenomXPInspiralArrays`): `v = (pi M f)^(1/3)`, the dimensionless spins
/// and the Newtonian orbital angular momentum direction in the source frame
/// at `f_ref`.
#[derive(Clone, Debug, Default)]
pub(crate) struct PnArrays {
    pub(crate) v: Vec<f64>,
    pub(crate) s1: [Vec<f64>; 3],
    pub(crate) s2: [Vec<f64>; 3],
    pub(crate) lnhat: [Vec<f64>; 3],
}

impl PnArrays {
    pub(crate) fn len(&self) -> usize {
        self.v.len()
    }
}

/// `SpinTaylorOrbit` series in a fixed order for stitching.
fn series(orbit: &SpinTaylorOrbit) -> [&Vec<f64>; 13] {
    [
        &orbit.v,
        &orbit.s1_x,
        &orbit.s1_y,
        &orbit.s1_z,
        &orbit.s2_x,
        &orbit.s2_y,
        &orbit.s2_z,
        &orbit.lnhat_x,
        &orbit.lnhat_y,
        &orbit.lnhat_z,
        &orbit.e1_x,
        &orbit.e1_y,
        &orbit.e1_z,
    ]
}

/// `appendTS` of `LALSimIMRPhenomX_precession.c`: the second series is
/// copied over the last two samples of the first one and the result has
/// length `len1 + len2 - 1`, so that it ends with a zero.
fn append_ts(start: &[f64], end: &[f64]) -> Vec<f64> {
    let origlen = start.len();
    let mut out = Vec::with_capacity(origlen + end.len() - 1);
    out.extend_from_slice(&start[..origlen - 2]);
    out.extend_from_slice(end);
    out.push(0.0);
    out
}

/// `IMRPhenomX_InspiralAngles_SpinTaylor` for black-hole binaries: evolve the
/// SpinTaylorT4 equations from `f_min` (Hz, already including the buffer
/// below the first frequency of the waveform) through `f_ref` up to
/// `(fRING + 8 fDAMP) / M` of the aligned-spin model.
///
/// Every failure is reported as
/// [`IMRPhenomXError::SpinTaylorFallbackToMsa`]: LALSimulation then falls
/// back to the MSA angles.
pub(crate) fn inspiral_angles_spintaylor(
    wf: &WaveformStruct,
    prec: &PrecessionStruct,
    f_min: f64,
) -> Result<PnArrays, IMRPhenomXError> {
    let fallback = IMRPhenomXError::SpinTaylorFallbackToMsa;
    let f_ref = wf.f_ref;
    if f_min <= 0.0 {
        return Err(fallback("fmin must be positive"));
    }
    if f_ref < f_min {
        return Err(fallback("fRef must be >= fmin"));
    }
    let (m1_si, m2_si) = (wf.m1_si, wf.m2_si);
    let pi_gm = PI * (m1_si + m2_si) * (G_SI / C_SI) / (C_SI * C_SI);
    let binary = SpinTaylorT4::new(m1_si, m2_si, prec.chi1, prec.chi2);

    // For versions 32* to work, the integration must start in the inspiral.
    let f_meco_hz = mf_to_hz(wf.f_meco, wf.m_tot);
    let mut f_min = f_min;
    if f_min > f_meco_hz && (prec.version == 320 || prec.version == 321) {
        f_min = f_meco_hz;
    }
    let f_cut = mf_to_hz(prec.f_ring_aligned + 8. * prec.f_damp_aligned, wf.m_tot);
    // PhenomXPSpinTaylorCoarseFactor, read as a float.
    let coarse_fac = 10.0f32 as f64;
    let delta_t_coarse = 0.5 * coarse_fac / f_cut;

    let evolve = |binary: &SpinTaylorT4, dt: f64, f_start: f64, f_end: f64| {
        binary
            .evolve_orbit(dt, f_start, f_end)
            .map_err(|_| fallback("XLALSimInspiralSpinTaylorPNEvolveOrbit failed"))
    };

    let coarse: [Vec<f64>; 13] = if f_ref < REAL4_EPS || (f_ref - f_min).abs() < REAL4_EPS {
        // f_ref is (re)set to f_min: a single forward integration.
        let orbit = evolve(&binary, delta_t_coarse, f_min, f_cut)?;
        series(&orbit).map(|s| s.clone())
    } else if (f_ref - f_min) > REAL4_EPS {
        // Integrate backward to fStart, then forward to the end of the waveform.
        let backward = evolve(&binary, delta_t_coarse, f_ref, f_min - 0.5)?;
        if backward.len() <= 1 {
            return Err(fallback("backward PN integration returned a single sample"));
        }
        let forward = evolve(&binary, delta_t_coarse, f_ref, f_cut)?;
        let (b, f) = (series(&backward), series(&forward));
        std::array::from_fn(|k| append_ts(b[k], f[k]))
    } else {
        return Err(fallback("fRef must be >= fmin"));
    };

    // At high frequencies, perform a fine integration, starting some bins
    // before the last one.
    let len_low = coarse[0].len();
    let nbuffer = 9.min(len_low - 1);
    let copy_length = len_low - 1 - nbuffer;
    let at = |k: usize| coarse[k][copy_length];
    let f_trans = at(0).powf(3.) / pi_gm;
    if f_trans.is_nan() || f_trans <= 0. {
        return Err(fallback(
            "transition frequency in PN integration is not positive",
        ));
    }
    let fine_binary = SpinTaylorT4::new(m1_si, m2_si, [at(1), at(2), at(3)], [at(4), at(5), at(6)])
        .with_orientation([at(7), at(8), at(9)], [at(10), at(11), at(12)]);
    let fine = evolve(&fine_binary, 0.5 / f_cut, f_trans, f_cut)?;
    let len_pn = copy_length + fine.len();
    if len_pn < 4 {
        return Err(fallback("too few points for spline interpolation"));
    }
    let fine = series(&fine);
    let stitch = |k: usize| {
        let mut out = Vec::with_capacity(len_pn);
        out.extend_from_slice(&coarse[k][..copy_length]);
        out.extend_from_slice(fine[k]);
        out
    };
    let arrays = PnArrays {
        v: stitch(0),
        s1: [stitch(1), stitch(2), stitch(3)],
        s2: [stitch(4), stitch(5), stitch(6)],
        lnhat: [stitch(7), stitch(8), stitch(9)],
    };

    // The first frequency node must lie below the requested fmin to avoid
    // interpolation errors.
    let f_min_pn = arrays.v[0].powf(3.) / pi_gm;
    if f_min_pn < 0. || f_min_pn > f_min {
        return Err(fallback("PN solution starts above fmin"));
    }
    Ok(arrays)
}

/// `evaluate_QNMfit_fdamp21` of `LALSimIMRPhenomXHM_qnm.c`, the damping
/// frequency of the (2,1) mode (accumulated with fused multiply-adds like
/// the arm64 builds of LALSimulation, see `qnm.rs`).
fn fdamp21(af: f64) -> f64 {
    let x2 = af * af;
    let x3 = x2 * af;
    let x4 = x2 * x2;
    let x5 = x3 * x2;
    let num = [
        (-3.1358071947583093, af),
        (0.14456081596393977, x2),
        (1.2194717985037946, x3),
        (-0.2947372598589144, x4),
        (0.002943057145913646, x5),
    ]
    .iter()
    .fold(2.0696914454467294, |acc: f64, &(c, x): &(f64, f64)| {
        c.mul_add(x, acc)
    });
    let den = [
        (-219.81790388304876, af),
        (17.7141194900164, x2),
        (75.90115083917898, x3),
        (-18.975287709794745, x4),
    ]
    .iter()
    .fold(146.1779212636481, |acc: f64, &(c, x): &(f64, f64)| {
        c.mul_add(x, acc)
    });
    num / den
}

/// Coefficients of the merger-ringdown continuation of `alpha`,
/// `PhenomXPalphaMRD`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct AlphaMrd {
    a: f64,
    b: f64,
    c: f64,
}

impl AlphaMrd {
    /// `alphaMRD_coeff`: match `-(a + b / Mf^4 + c / Mf^2)` to the spline
    /// (value and derivative at `0.97 fmaxPN`, value at `0.99 fmaxPN`).
    fn new(spline: &CubicSpline, fmax_pn: f64) -> Result<Self, IMRPhenomXError> {
        let ftrans = fmax_pn;
        let f1 = 0.97 * ftrans;
        let f2 = 0.99 * ftrans;
        let f1sq = f1 * f1;
        let f2sq = f2 * f2;
        let f1cube = f1sq * f1;
        if !(spline.contains(f1) && spline.contains(f2)) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "alpha could not be interpolated for the MRD continuation",
            ));
        }
        let alpha1 = -spline.eval(f1);
        let dalpha1 = -spline.derivative(f1);
        let alpha2 = -spline.eval(f2);
        // The expressions of LAL with clang's contraction into fused
        // multiply-adds (see `crate::spintaylor::fp`).
        let f1_4 = f1.powf(4.);
        let f2_4 = f2.powf(4.);
        let den = 2. * (f1sq - f2sq).powf(2.);
        let a = fma(
            2. * f2_4,
            alpha2,
            fma(
                f1cube * (f1 - f2) * (f1 + f2),
                dalpha1,
                2. * fnms(2. * f1sq, f2sq, f1_4) * alpha1,
            ),
        ) / den;
        let b = (f1_4
            * f2sq
            * fma(
                f1 * (f1 - f2) * (f1 + f2),
                dalpha1,
                2. * f2sq * (-alpha1 + alpha2),
            ))
            / den;
        let c = (f1sq * fma(f1 * (-f1_4 + f2_4), dalpha1, 4. * f2_4 * (alpha1 - alpha2))) / den;
        Ok(Self { a, b, c })
    }

    /// `alphaMRD`.
    #[inline]
    pub(crate) fn alpha(&self, mf: f64) -> f64 {
        let mf2 = mf * mf;
        let mf4 = mf2 * mf2;
        -(self.a + self.b / mf4 + self.c / mf2)
    }

    /// `dalphaMRD`.
    #[inline]
    pub(crate) fn dalpha(&self, mf: f64) -> f64 {
        let mf2 = mf * mf;
        let mf3 = mf2 * mf;
        let mf5 = mf2 * mf3;
        4. * self.b / mf5 + 2. * self.c / mf3
    }
}

/// Coefficients of the merger-ringdown continuation of `beta`,
/// `PhenomXPbetaMRD`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct BetaMrd {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    cosbeta_sign: f64,
    flat: bool,
    dfdamp: f64,
}

impl BetaMrd {
    /// `betaMRD_coeff`: match
    /// `exp(-kappa Mf) / Mf (a / Mf + b / Mf^2 + c / Mf^3) + d` to the spline
    /// of `cos(beta)` at `0.97 fmaxPN` and `0.98 fmaxPN`, with
    /// `kappa = 2 pi (fdamp21 - fdamp22)` and `d = 0` or `pi`.
    fn new(spline: &CubicSpline, fmax_pn: f64, fmax_inspiral: f64, wf: &WaveformStruct) -> Self {
        let fdamp = wf.f_damp;
        let dfdamp = fdamp21(wf.a_final) / wf.m_final - fdamp;
        let kappa = 2. * PI * dfdamp;
        let f1 = 0.97 * fmax_pn;
        let f2 = 0.98 * fmax_pn;
        let f1sq = f1 * f1;
        let f2sq = f2 * f2;
        let ef1 = (kappa * f1).exp();
        let ef2 = (kappa * f2).exp();

        // As in LAL, only the status of the last evaluation is checked.
        let eval = |x: f64| {
            if spline.contains(x) {
                spline.eval(x)
            } else {
                f64::NAN
            }
        };
        let cosbeta1 = eval(f1);
        let dcosbeta2 = if spline.contains(f2) {
            spline.derivative(f2)
        } else {
            f64::NAN
        };
        let cosbeta2 = eval(f2);
        let success = spline.contains(fmax_inspiral);
        let cosbetamax = eval(fmax_inspiral);

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
            let off = if cosbetamax < 0. { PI } else { 0. };
            let den = (f1 - f2).powf(2.);
            // The expressions of LAL with clang's contraction into fused
            // multiply-adds (see `crate::spintaylor::fp`).
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
            beta.flat = false;
        }
        beta
    }

    /// `betaMRD`.
    #[inline]
    pub(crate) fn beta(&self, mf: f64) -> f64 {
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

/// The SpinTaylor Euler angles of one binary (the angle-related members of
/// `IMRPhenomXPrecessionStruct`): splines of `alpha` and `cos(beta)` over
/// the PN solution, the merger-ringdown continuation and, once
/// [`SpinTaylorAngles::interpolate_gamma`] is called, the spline of `gamma`.
#[derive(Clone, Debug)]
pub(crate) struct SpinTaylorAngles {
    /// `alpha` (unwrapped) and `cos(beta)` against `Mf`.
    pub(crate) alpha_spline: CubicSpline,
    pub(crate) cosbeta_spline: CubicSpline,
    /// `gamma` against `Mf` on a uniform grid.
    pub(crate) gamma_spline: Option<CubicSpline>,
    /// End of the inspiral angles and start of the MRD continuation.
    #[allow(dead_code)] // kept as in LAL
    pub(crate) fmax_inspiral: f64,
    pub(crate) ftrans_mrd: f64,
    pub(crate) alpha_mrd: AlphaMrd,
    pub(crate) beta_mrd: BetaMrd,
    /// `alpha` and `gamma` at the reference frequency.
    pub(crate) alpha_ref: f64,
    pub(crate) gamma_ref: f64,
    /// Precession version (320 or 321 continue the angles analytically).
    version: u32,
}

impl SpinTaylorAngles {
    /// `IMRPhenomX_InterpolateAlphaBeta_SpinTaylor` with
    /// `PhenomXPFinalSpinMod` 4: also sets the final spin of `wf` from the
    /// spins at the end of the PN integration and updates `fRING`, `fDAMP`.
    pub(crate) fn new(
        wf: &mut WaveformStruct,
        prec: &PrecessionStruct,
        arrays: &PnArrays,
    ) -> Result<Self, IMRPhenomXError> {
        let len_pn = arrays.len();
        let mut fgw = Vec::with_capacity(len_pn);
        let mut alphaaux = Vec::with_capacity(len_pn);
        let mut cosbeta = Vec::with_capacity(len_pn);
        let (cz1, sz1) = ((-prec.phi_j_sf).cos(), (-prec.phi_j_sf).sin());
        let (cy, sy) = ((-prec.theta_j_sf).cos(), (-prec.theta_j_sf).sin());
        let (cz2, sz2) = ((-prec.kappa).cos(), (-prec.kappa).sin());
        let mut mfmax_pn = 0.;
        for i in 0..len_pn {
            let lnhat = [arrays.lnhat[0][i], arrays.lnhat[1][i], arrays.lnhat[2][i]];
            let lnhat = rotate_z_cs(cz2, sz2, rotate_y_cs(cy, sy, rotate_z_cs(cz1, sz1, lnhat)));
            let fgw_hz = arrays.v[i].powf(3.) / prec.pi_gm;
            let fgw_mf = hz_to_mf(fgw_hz, wf.m_tot);
            if fgw_hz.is_nan() || fgw_hz <= 0. {
                // LAL leaves the remaining nodes uninitialised.
                return Err(IMRPhenomXError::SpinTaylorAngles(
                    "non-positive frequency in the PN solution",
                ));
            }
            // Euler angles in the J frame
            alphaaux.push(lnhat[1].atan2(lnhat[0]));
            cosbeta.push(lnhat[2]);
            fgw.push(fgw_mf);
            mfmax_pn = fgw_mf;
        }
        let i_max = len_pn - 1;
        if fgw
            .windows(2)
            .any(|w| w[0] >= w[1] || w[0].is_nan() || w[1].is_nan())
        {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "PN frequencies are not strictly increasing",
            ));
        }

        let mut fmax_inspiral = mfmax_pn - wf.delta_mf;
        if fmax_inspiral > wf.f_ring - wf.f_damp {
            fmax_inspiral = 1.020 * wf.f_meco;
        }
        let ftrans_mrd = 0.98 * fmax_inspiral;

        let alpha = unwrap(&alphaaux);
        let alpha_spline = CubicSpline::new(&fgw, &alpha);
        let cosbeta_spline = CubicSpline::new(&fgw, &cosbeta);
        if !cosbeta_spline.contains(fmax_inspiral) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "cos(beta) could not be evaluated at the end of the inspiral",
            ));
        }
        let cosbetamax = cosbeta_spline.eval(fmax_inspiral);

        // Estimate the final spin from the spins at the end of the PN
        // integration (PhenomXPFinalSpinMod 4).
        let m1 = wf.m1_si / wf.m_tot_si;
        let m2 = wf.m2_si / wf.m_tot_si;
        let at = |s: &[Vec<f64>; 3]| [s[0][i_max], s[1][i_max], s[2][i_max]];
        let dot = |a: [f64; 3], b: [f64; 3]| crate::spintaylor::fp::dot(&a, &b);
        let lnf = at(&arrays.lnhat);
        let lnorm = dot(lnf, lnf).sqrt();
        let s1f = at(&arrays.s1);
        let s2f = at(&arrays.s2);
        let dot_s1l = dot(s1f, lnf) / lnorm;
        let dot_s2l = dot(s2f, lnf) / lnorm;
        let perp = |s: [f64; 3], sl: f64, mass: f64| {
            let scale = mass * mass;
            std::array::from_fn::<f64, 3, _>(|k| scale * (s[k] - sl * lnf[k]))
        };
        let s1_perp = perp(s1f, dot_s1l, m1);
        let s2_perp = perp(s2f, dot_s2l, m2);
        let stot_perp: [f64; 3] = std::array::from_fn(|k| s1_perp[k] + s2_perp[k]);
        let s_perp_norm = dot(stot_perp, stot_perp).sqrt();
        let chi_perp_norm = s_perp_norm * (m1 + m2).powf(2.) / m1.powf(2.);
        let a_final = 1.0f64.copysign(cosbetamax)
            * precessing_final_spin_2017(wf.eta, dot_s1l, dot_s2l, chi_perp_norm);
        if a_final.is_nan() || a_final.abs() > 1.0 {
            // LAL evaluates the QNM fits outside of their domain.
            return Err(IMRPhenomXError::UnphysicalFinalSpin(a_final));
        }
        wf.set_final_spin(a_final);

        let (alpha_mrd, beta_mrd) = if matches!(prec.version, 320 | 321 | 330) {
            (
                AlphaMrd::new(&alpha_spline, fmax_inspiral)?,
                BetaMrd::new(&cosbeta_spline, fmax_inspiral, fmax_inspiral, wf),
            )
        } else {
            (AlphaMrd::default(), BetaMrd::default())
        };

        Ok(Self {
            alpha_spline,
            cosbeta_spline,
            gamma_spline: None,
            fmax_inspiral,
            ftrans_mrd,
            alpha_mrd,
            beta_mrd,
            alpha_ref: 0.0,
            gamma_ref: 0.0,
            version: prec.version,
        })
    }

    fn continues_analytically(&self) -> bool {
        self.version == 320 || self.version == 321
    }

    /// `gamma_from_alpha_cosbeta`: the change of `gamma` from `Mf` to
    /// `Mf + delta_mf` from the minimal rotation condition
    /// `dgamma/df = -cos(beta) dalpha/df`, with Boole's rule. `hint` caches
    /// the spline interval between calls (`gsl_interp_accel`).
    #[inline]
    pub(crate) fn gamma_step(
        &self,
        mf: f64,
        delta_mf: f64,
        hint: &mut usize,
    ) -> Result<f64, IMRPhenomXError> {
        let mf_high = mf + delta_mf;
        let step = (mf_high - mf) * 0.25;
        let mut alphadotcosbeta = [0.0; 5];
        if mf <= self.ftrans_mrd {
            for (j, value) in alphadotcosbeta.iter_mut().enumerate() {
                let mf_aux = fma(j as f64, step, mf);
                if !self.alpha_spline.contains(mf_aux) {
                    return Err(IMRPhenomXError::SpinTaylorAngles(
                        "could not evaluate splines for alpha and/or gamma angles",
                    ));
                }
                // alpha and cos(beta) share their knots.
                *hint = self.alpha_spline.locate(mf_aux, *hint);
                let cosbeta = self.cosbeta_spline.eval_in(*hint, mf_aux);
                let alphadot = self.alpha_spline.derivative_in(*hint, mf_aux);
                *value = cosbeta * alphadot;
            }
        } else {
            for (j, value) in alphadotcosbeta.iter_mut().enumerate() {
                let mf_aux = fma(j as f64, step, mf);
                let cosbeta = self.beta_mrd.beta(mf_aux).cos();
                let alphadot = self.alpha_mrd.dalpha(mf_aux);
                *value = cosbeta * alphadot;
            }
        }
        let a = alphadotcosbeta;
        // 7 a0 + 32 a1 + 12 a2 + 32 a3 + 7 a4 as contracted by clang
        let sum = fma(
            7.,
            a[4],
            fma(32., a[3], fma(12., a[2], fma(7., a[0], 32. * a[1]))),
        );
        Ok(-2. * step / 45. * sum)
    }

    /// `IMRPhenomX_InterpolateGamma_SpinTaylor`: integrate `gamma` on the
    /// uniform grid `Mf_min + i dMf` (with the spacing of
    /// `get_deltaF_from_wfstruct`) covering `[f_min, f_max]` (Hz), build its
    /// spline and evaluate `gamma` at the reference frequency.
    ///
    /// If the spline will only be evaluated up to `mf_needed`, the grid is
    /// truncated [`GAMMA_SPLINE_MARGIN`] knots beyond it. The coefficients of
    /// a natural cubic spline depend on a knot `k` positions away through a
    /// factor `~(2 - sqrt(3))^k`, `~1e-37` for 64 knots, so the spline below
    /// `mf_needed` is unchanged to the last bit while the (expensive)
    /// merger-ringdown part of the integration is skipped.
    pub(crate) fn interpolate_gamma(
        &mut self,
        wf: &WaveformStruct,
        f_min: f64,
        f_max: f64,
        mf_needed: Option<f64>,
    ) -> Result<(), IMRPhenomXError> {
        let mf_min = hz_to_mf(f_min, wf.m_tot);
        if f_min <= 0.0 || f_max <= 0.0 || f_max <= f_min || wf.f_ref < f_min {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "SpinTaylor angles need 0 < f_min <= f_ref and f_min < f_max",
            ));
        }
        let delta_mf = delta_mf_from_wfstruct(wf);
        let delta_f = mf_to_hz(delta_mf, wf.m_tot);
        // length of the frequency series expected by the user
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max / delta_f) as usize + 1;
        let output_length = i_stop - i_start;
        if output_length < 4 {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "too few points for the spline interpolation of gamma",
            ));
        }
        let mf_last = fma((output_length - 1) as f64, delta_mf, mf_min);
        if wf.mf_ref.is_nan() || wf.mf_ref < mf_min || wf.mf_ref > mf_last {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "gamma could not be evaluated at the reference frequency",
            ));
        }
        let n_knots = match mf_needed {
            Some(mf_needed) => {
                let last_needed = ((mf_needed - mf_min) / delta_mf).ceil().max(0.0) as usize;
                (last_needed + 1 + GAMMA_SPLINE_MARGIN).min(output_length)
            }
            None => output_length,
        };
        let mut frequencies = Vec::with_capacity(n_knots);
        let mut gamma = Vec::with_capacity(n_knots);
        frequencies.push(mf_min);
        gamma.push(0.);
        let mut hint = 0;
        for i in 1..n_knots {
            let mf = fma(i as f64, delta_mf, mf_min);
            frequencies.push(mf);
            let previous = gamma[i - 1];
            let deltagamma = if mf < self.ftrans_mrd || self.continues_analytically() {
                self.gamma_step(mf, delta_mf, &mut hint)?
            } else {
                0.0
            };
            gamma.push(previous + deltagamma);
        }
        let spline = CubicSpline::new(&frequencies, &gamma);
        self.gamma_ref = spline.eval_unchecked(wf.mf_ref);
        self.gamma_spline = Some(spline);
        Ok(())
    }

    /// `alpha` at the reference frequency: from the spline in the inspiral,
    /// from the MRD continuation for versions 32*, otherwise the value at
    /// the transition frequency.
    pub(crate) fn set_alpha_ref(&mut self, mf_ref: f64) -> Result<(), IMRPhenomXError> {
        let error = IMRPhenomXError::SpinTaylorAngles("could not evaluate alpha(f_ref)");
        self.alpha_ref = if mf_ref <= self.ftrans_mrd {
            if !self.alpha_spline.contains(mf_ref) {
                return Err(error);
            }
            self.alpha_spline.eval(mf_ref)
        } else if self.continues_analytically() {
            self.alpha_mrd.alpha(mf_ref)
        } else {
            if !self.alpha_spline.contains(self.ftrans_mrd) {
                return Err(error);
            }
            self.alpha_spline.eval(self.ftrans_mrd)
        };
        Ok(())
    }

    /// `IMRPhenomXPSpinTaylorAnglesIMR`: build the spline of `gamma` between
    /// the first and last frequencies and return `(alpha, cos(beta), gamma)`
    /// on the grid `frequencies` (Hz, increasing), with `alpha` offset to
    /// `alpha0` at the reference frequency and `gamma` not yet offset (the
    /// twisting up subtracts `gamma_ref`).
    ///
    /// In the merger-ringdown `gamma` is integrated from point to point,
    /// using (like LAL) the step to the previous frequency forward from the
    /// current one.
    ///
    /// Angles are only returned up to the last frequency with `Mf <= mf_cut`
    /// (the waveform vanishes above), which leaves those values unchanged.
    pub(crate) fn angles_imr(
        &mut self,
        wf: &WaveformStruct,
        prec: &PrecessionStruct,
        frequencies: &[f64],
        mf_cut: f64,
    ) -> Result<EulerAngles, IMRPhenomXError> {
        let n = frequencies.len();
        let f_min = frequencies[0];
        let f_max = frequencies[n - 1];
        let mf_min = hz_to_mf(f_min, wf.m_tot);
        if f_min <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies("fmin must be positive"));
        }
        if f_max <= f_min {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "fmax must be larger than fmin",
            ));
        }
        if wf.f_ref < f_min {
            return Err(IMRPhenomXError::InvalidFrequencies("fRef must be >= fmin"));
        }
        // The gamma spline is evaluated below ftrans_MRD and at f_ref.
        self.interpolate_gamma(wf, f_min, f_max, Some(self.ftrans_mrd.max(wf.mf_ref)))?;
        if !self.alpha_spline.contains(mf_min) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "could not evaluate alpha(f_min)",
            ));
        }
        let alphamin = self.alpha_spline.eval(mf_min);
        let cosbetamin = self.cosbeta_spline.eval(mf_min);
        self.set_alpha_ref(wf.mf_ref)?;
        let gamma_spline = self.gamma_spline.as_ref().expect("gamma spline");

        let alpha_offset = -self.alpha_ref + prec.alpha0;
        let mut angles = EulerAngles::with_capacity(n);
        angles.alpha.push(alphamin - self.alpha_ref + prec.alpha0);
        angles.cos_beta.push(cosbetamin);
        angles.gamma.push(0.);

        let mut hint = 0;
        let mut gamma_hint = 0;
        let mut step_hint = 0;
        for i in 1..n {
            let mf = hz_to_mf(frequencies[i], wf.m_tot);
            if mf > mf_cut {
                break;
            }
            if mf < self.ftrans_mrd {
                let gamma_ok = gamma_spline.contains(mf);
                if !(self.alpha_spline.contains(mf) && gamma_ok) {
                    return Err(IMRPhenomXError::SpinTaylorAngles(
                        "interpolation of the SpinTaylor angles out of range",
                    ));
                }
                hint = self.alpha_spline.locate(mf, hint);
                gamma_hint = gamma_spline.locate(mf, gamma_hint);
                angles
                    .alpha
                    .push(self.alpha_spline.eval_in(hint, mf) + alpha_offset);
                angles.cos_beta.push(self.cosbeta_spline.eval_in(hint, mf));
                angles.gamma.push(gamma_spline.eval_in(gamma_hint, mf));
            } else if self.continues_analytically() {
                let delta_mf = hz_to_mf(frequencies[i] - frequencies[i - 1], wf.m_tot);
                angles.alpha.push(self.alpha_mrd.alpha(mf) + alpha_offset);
                angles.cos_beta.push(self.beta_mrd.beta(mf).cos());
                let previous = angles.gamma[i - 1];
                // A failed step keeps the previous value.
                let deltagamma = self.gamma_step(mf, delta_mf, &mut step_hint).unwrap_or(0.0);
                angles.gamma.push(previous + deltagamma);
            } else {
                // continue the angles into the ringdown with constant values
                angles.alpha.push(angles.alpha[i - 1]);
                angles.cos_beta.push(angles.cos_beta[i - 1]);
                angles.gamma.push(angles.gamma[i - 1]);
            }
        }
        Ok(angles)
    }
}

/// Euler angles on a frequency grid.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct EulerAngles {
    /// Precession angle `alpha` of L about J.
    pub alpha: Vec<f64>,
    /// Cosine of the opening angle `beta` between L and J.
    pub cos_beta: Vec<f64>,
    /// Third Euler angle `gamma` (minimal rotation condition).
    pub gamma: Vec<f64>,
}

impl EulerAngles {
    fn with_capacity(n: usize) -> Self {
        Self {
            alpha: Vec::with_capacity(n),
            cos_beta: Vec::with_capacity(n),
            gamma: Vec::with_capacity(n),
        }
    }
}

/// `XLALSimIMRPhenomXUnwrapArray`.
pub(crate) fn unwrap(input: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(input.len());
    out.push(input[0]);
    for i in 1..input.len() {
        let mut d = input[i] - input[i - 1];
        d = if d > PI {
            d - 2. * PI
        } else if d < -PI {
            d + 2. * PI
        } else {
            d
        };
        out.push(out[i - 1] + d);
    }
    out
}
