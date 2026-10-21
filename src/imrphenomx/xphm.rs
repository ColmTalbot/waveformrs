// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! The precessing higher-mode model IMRPhenomXPHM, following
//! `XLALSimIMRPhenomXPHMFrequencySequence`, `XLALSimIMRPhenomXPHM`,
//! `IMRPhenomXPHM_hplushcross` and `IMRPhenomXPHMTwistUp` in
//! `LALSimIMRPhenomXPHM.c`, and `IMRPhenomX_Initialize_Euler_Angles` in
//! `LALSimIMRPhenomX_precession.c`, without multibanding.

use std::f64::consts::PI;

use num_complex::Complex;

use crate::cmath::LibmCbrt;
use crate::constants::{MSUN_SI, PC_SI};
use crate::grid::{FrequencyGrid, ScalePowers};
use crate::imrphenomx::internals::{hz_to_mf, IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::phase::PhaseCoefficients;
use crate::imrphenomx::precessing_grid::{
    add_twisted_mode, angles_with_reuse, coprecessing_mode, coprecessing_mode_generator,
    rotate_polarizations, AngleArrays, TwistArrays,
};
use crate::imrphenomx::precession::{
    mf_to_hz, PrecessingBinary, PrecessionStruct, PrecessionVersion,
};
use crate::imrphenomx::spintaylor_angles::{inspiral_angles_spintaylor, SpinTaylorAngles};
use crate::imrphenomx::version_fallback::{check_version, forbid_spintaylor_fallback};
use crate::imrphenomx::wigner::half_angles;
use crate::imrphenomx::xas::TidalGrid22;
use crate::imrphenomx::xhm::{ModeGenerator, IMRPHENOMXHM_MODES};
use crate::imrphenomx::xp::{frequency_bounds, Polarizations, TidalInput};
use crate::nrtidal::NRTidalVersion;

pub(super) mod jframe;

/// `pPrec->sqrt*` as hard-coded in `IMRPhenomXGetAndSetPrecessionVariables`
/// (`sqrt10` is only accurate to 15 digits there).
// LAL's literal 1.4142135623730951 is the double nearest to sqrt(2).
const SQRT2: f64 = std::f64::consts::SQRT_2;
const SQRT5: f64 = 2.236_067_977_499_79;
const SQRT6: f64 = 2.449_489_742_783_178;
const SQRT7: f64 = 2.645_751_311_064_590_7;
const SQRT10: f64 = 3.162_277_660_168_38;
const SQRT14: f64 = 3.741_657_386_773_941_3;
const SQRT15: f64 = 3.872_983_346_207_417;
const SQRT30: f64 = 5.477_225_575_051_661;
const SQRT70: f64 = 8.366_600_265_340_756;

/// IMRPhenomXPHM (Pratten et al. 2021, arXiv:2004.06503): the aligned-spin
/// modes of IMRPhenomXHM, `(l, |m'|)` = (2,2), (2,1), (3,3), (3,2), (4,4),
/// twisted up with precession Euler angles, for the default LALSimulation
/// settings (`PhenomXPConvention` 1, `PhenomXPFinalSpinMod` 4) without
/// multibanding (`PhenomXHMThresholdMband = PhenomXPHMThresholdMband = 0`),
/// with the SpinTaylor angles (`PhenomXPrecVersion` 320, arXiv:2412.16721),
/// the MSA angles (223) or the NNLO angles (102).
///
/// # Conventions
///
/// As for [`crate::imrphenomx::IMRPhenomXP`]: [`IMRPhenomXPHM::polarizations`]
/// reproduces `SimInspiralChooseFDWaveformSequence(phi_ref, m1, m2, *spin_1,
/// *spin_2, f_ref, distance, inclination, params, IMRPhenomXPHM,
/// frequencies)` and [`IMRPhenomXPHM::polarizations_uniform`]
/// `SimInspiralChooseFDWaveform` with `deltaF > 0`, where `params` sets
/// `PhenomXPrecVersion` and switches multibanding off. Masses are in solar
/// masses, the distance in Mpc and the spins are the dimensionless spin
/// vectors in the L frame at the reference frequency.
///
/// The co-precessing modes are twisted up at the frequency `2 Mf / m'` of
/// the angles; each mode `(l, m')` contributes both `m'` and `-m'` (LAL's
/// default mode array). [`IMRPhenomXPHM::with_mode_array`] selects a subset,
/// which (as in LAL) also changes the frequency range of the SpinTaylor
/// angles.
///
/// # Precomputation and failures
///
/// With a reference frequency, [`IMRPhenomXPHM::new`] precomputes the
/// precession quantities and, for the NNLO and MSA angles, all mode
/// coefficients. The SpinTaylor angles (and the remnant spin they set, and
/// hence the mode coefficients) depend on the grid in LALSimulation and are
/// computed per call.
///
/// Like LALSimulation, the SpinTaylor angles fall back to the NNLO angles for
/// in-plane spins below `1e-7` and silently to the MSA angles when the PN
/// integration does not provide them (which depends on the grid and, through
/// the range of `m'`, on the mode array); the MSA angles fall back to the NNLO
/// angles when the MSA system cannot be initialised.
#[derive(Clone, Debug)]
pub struct IMRPhenomXPHM {
    binary: PrecessingBinary,
    version: PrecessionVersion,
    /// Co-precessing modes `(l, m' > 0)` in the order of LAL's mode loop.
    modes: Vec<(u32, u32)>,
    /// Whether the mode array was set with [`IMRPhenomXPHM::with_mode_array`]
    /// (a `ModeArray` in the LALDict), which sets the frequency range of the
    /// SpinTaylor angles of the J-frame modes.
    custom_mode_array: bool,
    /// Whether the prescription may change as in LAL (see
    /// [`IMRPhenomXPHM::with_version_fallback`]).
    version_fallback: bool,
    /// Parameter-only set-up for a non-uniform grid (`f_ref > 0`).
    setup: Option<Setup>,
    /// The NRTidal corrections of the (2,2) mode (only for the uniform-grid
    /// IMRPhenomXP_NRTidalv2/v3, see [`IMRPhenomXPHM::with_nrtidal`]).
    tidal: Option<TidalInput>,
}

/// One co-precessing mode ready for evaluation.
#[derive(Clone, Debug)]
struct Mode {
    ell: u32,
    emm: u32,
    generator: ModeGenerator,
}

/// Waveform and precession structs and, when they do not depend on the
/// grid, the mode generators.
#[derive(Clone, Debug)]
struct Setup {
    wf: WaveformStruct,
    prec: PrecessionStruct,
    modes: Option<Vec<Mode>>,
}

/// The twisting-up factors shared by all modes with the same `m'` at one
/// frequency.
#[derive(Clone, Copy, Debug, Default)]
struct Twist {
    /// `cos(k alpha)`, `sin(k alpha)` for `k = 0..4`.
    cos_k: [f64; 5],
    sin_k: [f64; 5],
    /// Powers `0..=8` of `cos(beta/2)` and `sin(beta/2)`.
    c: [f64; 9],
    s: [f64; 9],
    /// `exp(-i m' epsilon) / 2`.
    eps: Complex<f64>,
}

impl Twist {
    fn new(alpha: f64, cos_beta: f64, epsilon: f64, mprime: u32) -> Self {
        // Powers of exp(i alpha) by multiplication, as in LAL.
        let (sin_a, cos_a) = alpha.sin_cos();
        let e1 = Complex::new(cos_a, sin_a);
        let e2 = e1 * e1;
        let e3 = e1 * e2;
        let e4 = e1 * e3;
        let powers = [Complex::new(1.0, 0.0), e1, e2, e3, e4];
        let (cb, sb) = half_angles(cos_beta);
        let mut c = [1.0; 9];
        let mut s = [1.0; 9];
        for k in 1..9 {
            c[k] = cb * c[k - 1];
            s[k] = sb * s[k - 1];
        }
        // cBetah2 = cBetah * cBetah, ... as in LAL (the first product is
        // cBetah * 1.0 = cBetah exactly).
        let (sin_e, cos_e) = (mprime as f64 * epsilon).sin_cos();
        Self {
            cos_k: powers.map(|z| z.re),
            sin_k: powers.map(|z| z.im),
            c,
            s,
            eps: Complex::new(cos_e, -sin_e) * 0.5,
        }
    }
}

/// `d^l_{k, m'}(beta)` for `k = -l..l` (Appendix A of arXiv:2004.06503),
/// written as in `IMRPhenomXPHMTwistUp`.
#[inline(always)]
pub(crate) fn wigner_d(ell: u32, mprime: u32, c: &[f64; 9], s: &[f64; 9]) -> [f64; 9] {
    let mut d = [0.0; 9];
    match (ell, mprime) {
        (2, 2) => {
            d[..5].copy_from_slice(&[
                s[4],
                2.0 * c[1] * s[3],
                SQRT6 * s[2] * c[2],
                2.0 * c[3] * s[1],
                c[4],
            ]);
        }
        (2, 1) => {
            d[..5].copy_from_slice(&[
                2.0 * c[1] * s[3],
                3.0 * c[2] * s[2] - s[4],
                SQRT6 * (c[3] * s[1] - c[1] * s[3]),
                c[2] * (c[2] - 3.0 * s[2]),
                -2.0 * c[3] * s[1],
            ]);
        }
        (3, 3) => {
            d[..7].copy_from_slice(&[
                s[6],
                SQRT6 * c[1] * s[5],
                SQRT15 * c[2] * s[4],
                2.0 * SQRT5 * c[3] * s[3],
                SQRT15 * c[4] * s[2],
                SQRT6 * c[5] * s[1],
                c[6],
            ]);
        }
        (3, 2) => {
            d[..7].copy_from_slice(&[
                SQRT6 * c[1] * s[5],
                s[4] * (5.0 * c[2] - s[2]),
                SQRT10 * s[3] * (2.0 * c[3] - c[1] * s[2]),
                SQRT30 * c[2] * (c[2] - s[2]) * s[2],
                SQRT10 * c[3] * (c[2] * s[1] - 2.0 * s[3]),
                c[4] * (c[2] - 5.0 * s[2]),
                -1.0 * SQRT6 * c[5] * s[1],
            ]);
        }
        (4, 4) => {
            d.copy_from_slice(&[
                s[8],
                2.0 * SQRT2 * c[1] * s[7],
                2.0 * SQRT7 * c[2] * s[6],
                2.0 * SQRT14 * c[3] * s[5],
                SQRT70 * c[4] * s[4],
                2.0 * SQRT14 * c[5] * s[3],
                2.0 * SQRT7 * c[6] * s[2],
                2.0 * SQRT2 * c[7] * s[1],
                c[8],
            ]);
        }
        _ => unreachable!("mode ({ell},{mprime}) is not twisted up"),
    }
    d
}

/// `IMRPhenomXPHMTwistUp` for the co-precessing mode `h = h_{l,-m'}`:
/// returns the contributions `(h_+, h_x)` of the modes `(l, +-m')` in the J
/// frame. With `A_k = e^{-i k alpha} d^l_{k,-m'} Y_{lk}`,
/// `B_k = e^{i k alpha} d^l_{k,m'} Y_{lk}` and `sigma = (-1)^l`,
///
/// `h_+ = e^{-i m' epsilon} h / 2 sum_k (A_k + sigma B_k)`,
/// `h_x = e^{-i m' epsilon} h / 2 sum_k i (A_k - sigma B_k)`,
///
/// evaluated with real arithmetic (`Y_{lk}` is real for the azimuth zero).
#[inline]
fn twist_up(ell: u32, mprime: u32, twist: &Twist, ylm: &[f64], h: Complex<f64>) -> Polarization {
    let l = ell as usize;
    let d = wigner_d(ell, mprime, &twist.c, &twist.s);
    let sigma = if ell % 2 == 0 { 1.0 } else { -1.0 };
    let mut hp_re = 0.0;
    let mut hp_im = 0.0;
    let mut hc_re = 0.0;
    let mut hc_im = 0.0;
    for (index, &y) in ylm.iter().enumerate().take(2 * l + 1) {
        // k = index - l; d^l_{k,-m'} = (-1)^(k + m') d^l_{-k,m'}
        let k = index as i64 - l as i64;
        let dm = if (k + mprime as i64) % 2 == 0 {
            d[2 * l - index]
        } else {
            -d[2 * l - index]
        };
        let a = dm * y;
        let b = sigma * d[index] * y;
        let k_abs = k.unsigned_abs() as usize;
        let cos = twist.cos_k[k_abs];
        let sin = if k < 0 {
            -twist.sin_k[k_abs]
        } else {
            twist.sin_k[k_abs]
        };
        let sum = a + b;
        hp_re += sum * cos;
        hp_im += (b - a) * sin;
        hc_re += sum * sin;
        hc_im += (a - b) * cos;
    }
    let eps_h = twist.eps * h;
    (
        eps_h * Complex::new(hp_re, hp_im),
        eps_h * Complex::new(hc_re, hc_im),
    )
}

type Polarization = (Complex<f64>, Complex<f64>);

/// The `s = -2` spherical harmonics at `(thetaJN, 0)` for mode `l`.
fn harmonics(prec: &PrecessionStruct, ell: u32) -> &[f64] {
    match ell {
        2 => &prec.ylm.y2,
        3 => &prec.ylm.y3,
        _ => &prec.ylm.y4,
    }
}

/// Per-`m'` state of the SpinTaylor angle evaluation.
#[derive(Clone, Copy, Debug, Default)]
struct AngleCursor {
    /// Interval hints of the alpha/cos(beta) and gamma splines.
    hint: usize,
    gamma_hint: usize,
    step_hint: usize,
    /// `pPrec->gamma_in`, the previous value of gamma for this `m'`.
    gamma_in: f64,
}

impl Setup {
    /// `IMRPhenomXSetWaveformVariables` and
    /// `IMRPhenomXGetAndSetPrecessionVariables` as called by
    /// `XLALSimIMRPhenomXPHMFrequencySequence` (`delta_f = 0`, which passes
    /// the masses and spins after `XLALIMRPhenomXPCheckMassesAndSpins`) and
    /// `XLALSimIMRPhenomXPHM` (`delta_f > 0`, which passes them in the
    /// original order).
    fn new(
        binary: &PrecessingBinary,
        version: PrecessionVersion,
        modes: &[(u32, u32)],
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<Self, IMRPhenomXError> {
        Self::with_tidal(binary, version, None, modes, delta_f, f_min, f_max)
    }

    /// [`Setup::new`] with the NRTidal corrections of the (2,2) mode (LAL's
    /// `PhenomXTidalFlag`, as the uniform-grid IMRPhenomXP_NRTidalv2/v3 of
    /// `SimInspiralChooseFDWaveform` use it).
    fn with_tidal(
        binary: &PrecessingBinary,
        version: PrecessionVersion,
        tidal: Option<TidalInput>,
        modes: &[(u32, u32)],
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<Self, IMRPhenomXError> {
        let (m1, m2, chi1, chi2) = binary.ordered();
        let (wf_m1, wf_m2, wf_chi1, wf_chi2) = if delta_f > 0.0 {
            (
                binary.mass_1,
                binary.mass_2,
                binary.spin_1[2],
                binary.spin_2[2],
            )
        } else {
            (m1, m2, chi1[2], chi2[2])
        };
        let mut wf = WaveformStruct::from_si(
            wf_m1 * MSUN_SI,
            wf_m2 * MSUN_SI,
            wf_chi1,
            wf_chi2,
            delta_f,
            binary.f_ref,
            binary.phi_ref,
            f_min,
            f_max,
            binary.luminosity_distance * 1e6 * PC_SI,
            binary.inclination,
        )?;
        if let Some((version, lambda_1, lambda_2)) = tidal {
            wf.set_nrtidal(version, binary.mass_1, binary.mass_2, lambda_1, lambda_2)?;
        }
        // IMRPhenomX_GetandSetModes: the range of |m| in the mode array.
        let m_min = modes.iter().map(|&(_, m)| m).min().unwrap_or(4);
        let m_max = modes.iter().map(|&(_, m)| m).max().unwrap_or(1);
        let prec = PrecessionStruct::new(
            &mut wf,
            m1 * MSUN_SI,
            m2 * MSUN_SI,
            chi1,
            chi2,
            version,
            m_min,
            m_max,
        )?;
        let modes = (prec.precessing_tag != 3).then(|| build_modes(&wf, modes));
        Ok(Self { wf, prec, modes })
    }

    /// Update the frequency bounds of the waveform struct for a
    /// non-uniform grid, as `IMRPhenomXSetWaveformVariables` sets them.
    fn with_bounds(&self, f_min: f64, f_max: f64) -> Self {
        let mut setup = self.clone();
        let wf = &mut setup.wf;
        wf.f_min = f_min;
        wf.f_max = f_max;
        wf.mf_max = hz_to_mf(f_max, wf.m_tot);
        wf.f_max_prime = wf.f_max_prime(f_max);
        setup
    }

    /// `IMRPhenomX_Initialize_Euler_Angles` (with the SpinTaylor integration
    /// of `IMRPhenomXGetAndSetPrecessionVariables`) for multibanding off:
    /// the splines of alpha and cos(beta) and the spline of gamma between
    /// `(f_min - buffer) 2 / M_MAX` and
    /// `(max(f_max, fRING + 4 fDAMP) + buffer) 2 / M_MIN`. Also sets the
    /// final spin from the end of the PN integration. `fast` selects the
    /// vectorised merger-ringdown integration of gamma of the grid path.
    fn initialize_euler_angles(&mut self, fast: bool) -> Result<SpinTaylorAngles, IMRPhenomXError> {
        let wf = &mut self.wf;
        let prec = &self.prec;
        let buffer = prec.integration_buffer;
        // start below fMin to avoid interpolation artefacts
        let f_min_angles = (wf.f_min - buffer) * 2. / prec.m_max as f64;
        if f_min_angles.is_nan() || f_min_angles <= 0. {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "fMin is too low and numerical angles could not be computed",
            ));
        }
        let mf_max_angles =
            (wf.mf_max.max(wf.f_ring + 4. * wf.f_damp) + hz_to_mf(buffer, wf.m_tot)) * 2.
                / prec.m_min as f64;
        let f_max_angles = mf_to_hz(mf_max_angles, wf.m_tot);

        // The PN integration starts at the same frequency
        // (IMRPhenomXGetAndSetPrecessionVariables).
        let arrays = inspiral_angles_spintaylor(wf, prec, f_min_angles)?;
        // IMRPhenomX_SpinTaylorAnglesSplinesAll
        if f_max_angles <= f_min_angles || wf.f_ref < f_min_angles {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "SpinTaylor angles need f_min <= f_ref and f_min < f_max",
            ));
        }
        let mut angles = SpinTaylorAngles::new(wf, prec, arrays)?;
        // gamma is evaluated at Mf' = 2 Mf / m' up to the cutoff, at f_ref
        // and at ftrans_MRD: the integration can stop beyond these (see
        // SpinTaylorAngles::interpolate_gamma), which leaves the spline
        // unchanged where it is used.
        let mf_needed = (wf.f_max_prime * wf.m_sec * (2.0 / prec.m_min as f64))
            .max(wf.mf_ref)
            .max(angles.ftrans_mrd);
        angles.interpolate_gamma(wf, f_min_angles, f_max_angles, Some(mf_needed), fast)?;
        angles.set_alpha_ref(wf.mf_ref)?;
        // alpha, cos(beta) and gamma are evaluated at ftrans_MRD (used by
        // the versions with constant angles in the ringdown).
        let gamma_spline = angles.gamma_spline.as_ref().expect("gamma spline");
        if !(angles.alpha_spline.contains(angles.ftrans_mrd)
            && gamma_spline.contains(angles.ftrans_mrd))
        {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "could not compute the angles at the end of the inspiral",
            ));
        }
        Ok(angles)
    }

    /// `IMRPhenomXPHM_hplushcross` on the grid `frequencies`, including the
    /// rotation of the polarisations by `2 zeta`.
    fn generate(
        &mut self,
        mode_list: &[(u32, u32)],
        frequencies: &[f64],
    ) -> Result<Polarizations, IMRPhenomXError> {
        if self.wf.f_max_prime <= self.wf.f_min {
            return Err(IMRPhenomXError::InvalidFrequencies("fCut <= f_min"));
        }
        let angles = if self.prec.precessing_tag == 3 {
            Some(self.initialize_euler_angles(false)?)
        } else {
            None
        };
        let modes = match self.modes.take() {
            Some(modes) => modes,
            None => build_modes(&self.wf, mode_list),
        };
        let wf = &self.wf;
        let prec = &self.prec;
        let tidal22 = tidal_22(&modes, wf, frequencies);

        let n = frequencies.len();
        let zero = Complex::new(0.0, 0.0);
        let mut hp = vec![zero; n];
        let mut hc = vec![zero; n];
        let mf_cut = wf.f_max_prime * wf.m_sec;

        // The distinct m' of the active modes: the angles are evaluated
        // once per frequency and m'.
        let mut mprimes: Vec<u32> = modes.iter().map(|mode| mode.emm).collect();
        mprimes.sort_unstable();
        mprimes.dedup();
        let mut twists = [Twist::default(); 5];
        let mut cursors = [AngleCursor::default(); 5];

        for (i, &f) in frequencies.iter().enumerate() {
            let mf = wf.m_sec * f;
            // Do not generate the waveform above Mf_max (default Mf = 0.3)
            if mf > mf_cut {
                continue;
            }
            for &mprime in &mprimes {
                let (alpha, cos_beta, epsilon) = match &angles {
                    Some(angles) => spintaylor_angles::<true>(
                        angles,
                        wf,
                        prec,
                        mf,
                        mprime,
                        &mut cursors[mprime as usize],
                    )?,
                    None if prec.msa.is_some() => {
                        // IMRPhenomX_Return_phi_zeta_costhetaL_MSA at
                        // v = (pi Mf 2 / m')^(1/3), with the offsets of
                        // Get_alpha_epsilon_offset (the same for every m')
                        let v = (PI * mf * (2.0 / mprime as f64)).libm_cbrt();
                        let (alpha, epsilon, cos_beta) = prec.msa_angles(v);
                        (alpha, cos_beta, epsilon)
                    }
                    None => {
                        // Get_alpha_beta_epsilon at omega = pi Mf 2 / m'
                        let (alpha, epsilon, cos_beta) =
                            prec.nnlo_angles(PI * mf * 2. / mprime as f64, wf.eta);
                        (alpha, cos_beta, epsilon)
                    }
                };
                twists[mprime as usize] = Twist::new(alpha, cos_beta, epsilon, mprime);
            }
            let powers = UsefulPowers::new(mf);
            let mut sum_p = zero;
            let mut sum_c = zero;
            for mode in &modes {
                let h = match &tidal22 {
                    Some(tidal22) if (mode.ell, mode.emm) == (2, 2) => tidal22.h22(i, mf, &powers),
                    _ => mode.generator.h(mf, &powers),
                };
                let (p, c) = twist_up(
                    mode.ell,
                    mode.emm,
                    &twists[mode.emm as usize],
                    harmonics(prec, mode.ell),
                    h,
                );
                sum_p += p;
                sum_c += c;
            }
            hp[i] = sum_p;
            hc[i] = sum_c;
        }

        let zeta = prec.zeta_polarization;
        if zeta.abs() > 0.0 {
            let cos_pol = (2.0 * zeta).cos();
            let sin_pol = (2.0 * zeta).sin();
            for (p, c) in hp.iter_mut().zip(hc.iter_mut()) {
                let (pp, pc) = (*p, *c);
                *p = pp * cos_pol + pc * sin_pol;
                *c = pc * cos_pol - pp * sin_pol;
            }
        }
        Ok((hp, hc))
    }

    /// [`Self::generate`] on a prepared grid, written into `hplus` and
    /// `hcross` (see [`crate::imrphenomx::precessing_grid`]).
    fn generate_grid(
        &mut self,
        mode_list: &[(u32, u32)],
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<(), IMRPhenomXError> {
        if self.wf.f_max_prime <= self.wf.f_min {
            return Err(IMRPhenomXError::InvalidFrequencies("fCut <= f_min"));
        }
        let spintaylor = if self.prec.precessing_tag == 3 {
            Some(self.initialize_euler_angles(true)?)
        } else {
            None
        };
        let modes = match self.modes.take() {
            Some(modes) => modes,
            None => build_modes(&self.wf, mode_list),
        };
        let wf = &self.wf;
        let prec = &self.prec;
        let frequencies = grid.frequencies();
        let tidal22 = tidal_22(&modes, wf, frequencies);
        let mf_cut = wf.f_max_prime * wf.m_sec;
        // the frequencies with Mf <= mf_cut
        let n = frequencies.partition_point(|&f| wf.m_sec * f <= mf_cut);

        let mut mprimes: Vec<u32> = modes.iter().map(|mode| mode.emm).collect();
        mprimes.sort_unstable();
        mprimes.dedup();
        let angles = match &spintaylor {
            Some(angles) => spintaylor_angle_arrays(angles, wf, prec, &frequencies[..n], &mprimes)?,
            None if prec.msa.is_some() => {
                // Evaluate m' = 2 first: the other m' reuse its angles.
                mprimes.sort_by_key(|&mprime| mprime != 2);
                // v = (pi Mf 2 / m')^(1/3) as in the per-frequency path
                let argument = |mprime: u32, i: usize| {
                    PI * (wf.m_sec * frequencies[i]) * (2.0 / mprime as f64)
                };
                angles_with_reuse(&mprimes, n, argument, |arguments, angles| {
                    prec.msa_angles_many(arguments, |alpha, epsilon, cos_beta| {
                        angles.push(alpha, cos_beta, epsilon)
                    });
                })
            }
            None => std::array::from_fn(|mprime| {
                mprimes
                    .contains(&(mprime as u32))
                    .then(|| AngleArrays::nnlo(prec, wf.eta, grid, wf.m_sec, mprime as u32, n))
            }),
        };
        let twists: [Option<TwistArrays>; 5] =
            std::array::from_fn(|mprime| angles[mprime].as_ref().map(TwistArrays::new));

        let zero = Complex::new(0.0, 0.0);
        hplus.fill(zero);
        hcross.fill(zero);
        let mut h = vec![zero; n];
        let mut extra = Vec::new();
        let scale = ScalePowers::new(wf.m_sec);
        for mode in &modes {
            if mode.generator.vanishes() {
                continue;
            }
            let emm = mode.emm as usize;
            let epsilon = &angles[emm].as_ref().expect("angles of m'").epsilon;
            match &tidal22 {
                Some(tidal22) if (mode.ell, mode.emm) == (2, 2) => {
                    coprecessing_mode(&mut h, 2, epsilon, |i| {
                        let powers = grid.powers(&scale, i);
                        let mf = powers.itself;
                        tidal22.amplitude_phase(i, mf, &UsefulPowers::from_grid(&powers))
                    })
                }
                _ => coprecessing_mode_generator(
                    &mut h,
                    &mut extra,
                    mode.emm,
                    epsilon,
                    &mode.generator,
                    grid,
                    &scale,
                ),
            }
            let twist = twists[emm].as_ref().expect("twist of m'");
            let ylm = harmonics(prec, mode.ell);
            add_twisted_mode(mode.ell, mode.emm, ylm, twist, &h, hplus, hcross);
        }
        rotate_polarizations(prec.zeta_polarization, &mut hplus[..n], &mut hcross[..n]);
        Ok(())
    }
}

/// The SpinTaylor Euler angles of [`spintaylor_angles`] for each `m'` in
/// `mprimes` (increasing) at the frequencies `frequencies`, returning the
/// error the per-frequency loop of [`Setup::generate`] meets first.
/// `cos(beta)` in the merger-ringdown uses the vectorised functions of
/// [`SpinTaylorAngles::mrd_cos_beta_fast`] (agreeing to rounding).
fn spintaylor_angle_arrays(
    angles: &SpinTaylorAngles,
    wf: &WaveformStruct,
    prec: &PrecessionStruct,
    frequencies: &[f64],
    mprimes: &[u32],
) -> Result<[Option<AngleArrays>; 5], IMRPhenomXError> {
    let mut out: [Option<AngleArrays>; 5] = Default::default();
    // (frequency index, error) of the first failure in the order of the
    // per-frequency loop (frequencies outside, increasing m' inside)
    let mut first_error: Option<(usize, IMRPhenomXError)> = None;
    for &mprime in mprimes {
        let end = first_error.as_ref().map_or(frequencies.len(), |(i, _)| *i);
        let mut cursor = AngleCursor::default();
        let mut arrays = AngleArrays::with_capacity(frequencies.len());
        for (i, &f) in frequencies[..end].iter().enumerate() {
            let mf = wf.m_sec * f;
            match spintaylor_angles::<false>(angles, wf, prec, mf, mprime, &mut cursor) {
                Ok((alpha, cos_beta, epsilon)) => arrays.push(alpha, cos_beta, epsilon),
                Err(error) => {
                    first_error = Some((i, error));
                    break;
                }
            }
        }
        // cos(beta) in the merger-ringdown, Mf' >= ftrans_MRD (a suffix)
        let factor = 2.0 / mprime as f64;
        let n = arrays.cos_beta.len();
        let start =
            frequencies[..n].partition_point(|&f| wf.m_sec * f * factor < angles.ftrans_mrd);
        let mfprime: Vec<f64> = frequencies[start..n]
            .iter()
            .map(|&f| wf.m_sec * f * factor)
            .collect();
        angles.mrd_cos_beta_fast(&mfprime, &mut arrays.cos_beta[start..]);
        out[mprime as usize] = Some(arrays);
    }
    match first_error {
        Some((_, error)) => Err(error),
        None => Ok(out),
    }
}

/// The SpinTaylor Euler angles `(alpha, cos(beta), epsilon)` for mode `m'`
/// at the frequency `Mf' = Mf 2 / m'`, as in `IMRPhenomXPHMTwistUp`: from the
/// splines in the inspiral and the MRD continuation above `ftrans_MRD`,
/// with `gamma` from its spline or, beyond it, integrated from the previous
/// frequency. Without `MRD_COS_BETA`, `cos(beta)` is left as NaN in the
/// merger-ringdown (for the grid path, which evaluates it separately).
#[inline]
fn spintaylor_angles<const MRD_COS_BETA: bool>(
    angles: &SpinTaylorAngles,
    wf: &WaveformStruct,
    prec: &PrecessionStruct,
    mf: f64,
    mprime: u32,
    cursor: &mut AngleCursor,
) -> Result<(f64, f64, f64), IMRPhenomXError> {
    let mfprime = mf * (2.0 / mprime as f64);
    let gamma_spline = angles.gamma_spline.as_ref().expect("gamma spline");
    let (alpha_i, cos_beta, gamma) = if mfprime < angles.ftrans_mrd {
        if !(angles.alpha_spline.contains(mfprime) && gamma_spline.contains(mfprime)) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "failed to evaluate the angles",
            ));
        }
        cursor.hint = angles.alpha_spline.locate(mfprime, cursor.hint);
        cursor.gamma_hint = gamma_spline.locate(mfprime, cursor.gamma_hint);
        (
            angles.alpha_spline.eval_in(cursor.hint, mfprime),
            angles.cosbeta_spline.eval_in(cursor.hint, mfprime),
            gamma_spline.eval_in(cursor.gamma_hint, mfprime),
        )
    } else {
        let alpha_i = angles.alpha_mrd.alpha(mfprime);
        let cos_beta = if MRD_COS_BETA {
            angles.beta_mrd.beta(mfprime).cos()
        } else {
            f64::NAN
        };
        let gamma = if gamma_spline.contains(mfprime) {
            cursor.gamma_hint = gamma_spline.locate(mfprime, cursor.gamma_hint);
            gamma_spline.eval_in(cursor.gamma_hint, mfprime)
        } else {
            let delta_mf = wf.delta_mf * 2. / mprime as f64;
            cursor.gamma_in + angles.gamma_step(mfprime, delta_mf, &mut cursor.step_hint)?
        };
        (alpha_i, cos_beta, gamma)
    };
    cursor.gamma_in = gamma;
    let alpha_offset = angles.alpha_ref - prec.alpha0;
    let epsilon_offset = -angles.gamma_ref - prec.epsilon0;
    Ok((alpha_i - alpha_offset, cos_beta, -gamma - epsilon_offset))
}

/// `IMRPhenomXGetAndSetPrecessionVariables` falls back to the MSA angles
/// (version 223) when the SpinTaylor integration fails: redo the set-up with
/// [`PrecessionVersion::Msa`] and generate the polarisations on the same
/// grid. None of the SpinTaylor set-up enters the MSA angles, so this is
/// what LAL computes.
/// The (2,2) mode with the NRTidal corrections on the grid `frequencies`,
/// if the (2,2) mode is active and the waveform struct has them.
fn tidal_22(modes: &[Mode], wf: &WaveformStruct, frequencies: &[f64]) -> Option<TidalGrid22> {
    modes.iter().find_map(|mode| match &mode.generator {
        ModeGenerator::Xas22 { xas, .. } if xas.tidal.is_some() => {
            Some(TidalGrid22::new(xas, 2.0 * wf.phi0, frequencies))
        }
        _ => None,
    })
}

#[allow(clippy::too_many_arguments)]
fn fallback_to_msa(
    result: Result<Polarizations, IMRPhenomXError>,
    binary: &PrecessingBinary,
    tidal: Option<TidalInput>,
    modes: &[(u32, u32)],
    delta_f: f64,
    frequencies: &[f64],
    f_min: f64,
    f_max: f64,
) -> Result<Polarizations, IMRPhenomXError> {
    match result {
        Err(IMRPhenomXError::SpinTaylorFallbackToMsa(_)) => Setup::with_tidal(
            binary,
            PrecessionVersion::Msa,
            tidal,
            modes,
            delta_f,
            f_min,
            f_max,
        )?
        .generate(modes, frequencies),
        result => result,
    }
}

/// The mode generators of the active modes, in the order of LAL's mode
/// loop, skipping odd `m'` when they vanish (equal masses and spins).
fn build_modes(wf: &WaveformStruct, modes: &[(u32, u32)]) -> Vec<Mode> {
    let phase22 = PhaseCoefficients::new(wf);
    modes
        .iter()
        .filter(|&&(_, emm)| !(wf.q == 1.0 && wf.chi1l == wf.chi2l && emm % 2 != 0))
        .map(|&(ell, emm)| Mode {
            ell,
            emm,
            generator: ModeGenerator::new(wf, &phase22, ell, emm),
        })
        .collect()
}

impl IMRPhenomXPHM {
    /// IMRPhenomXPHM with the default SpinTaylor angles (`PhenomXPrecVersion`
    /// 320) and the default modes.
    ///
    /// Returns an error where LALSimulation would (masses, spins beyond the
    /// Kerr bound, non-positive distance, negative reference frequency,
    /// unphysical remnant fits).
    pub fn new(binary: &PrecessingBinary) -> Result<Self, IMRPhenomXError> {
        Self::with_version(binary, PrecessionVersion::default())
    }

    /// IMRPhenomXPHM with the given Euler angle prescription.
    pub fn with_version(
        binary: &PrecessingBinary,
        version: PrecessionVersion,
    ) -> Result<Self, IMRPhenomXError> {
        Self::build(binary, version, IMRPHENOMXHM_MODES.to_vec())
    }

    /// Restrict the waveform to the co-precessing modes `(l, m')` with
    /// `m' > 0` among (2,2), (2,1), (3,3), (3,2) and (4,4) (LAL's mode array
    /// with both signs of `m'`).
    pub fn with_mode_array(self, modes: &[(u32, u32)]) -> Result<Self, IMRPhenomXError> {
        if modes.is_empty() {
            return Err(IMRPhenomXError::InvalidModeArray);
        }
        for mode in modes {
            if !IMRPHENOMXHM_MODES.contains(mode) {
                return Err(IMRPhenomXError::InvalidModeArray);
            }
        }
        // LAL's mode loop order
        let ordered: Vec<(u32, u32)> = (2..=4u32)
            .flat_map(|ell| (1..=ell).map(move |emm| (ell, emm)))
            .filter(|mode| modes.contains(mode))
            .collect();
        let mut model = Self::build(&self.binary, self.version, ordered)?;
        model.custom_mode_array = true;
        model.version_fallback = self.version_fallback;
        Ok(model)
    }

    /// Whether the Euler-angle prescription may change automatically as in
    /// LALSimulation (the default, `true`): SpinTaylor to NNLO for in-plane
    /// spins below `1e-7`, SpinTaylor to MSA when the SpinTaylor angles
    /// cannot be computed, MSA to NNLO when the MSA system cannot be
    /// initialised. With `false` the requested prescription is used or
    /// the evaluation (polarizations and modes) returns
    /// [`IMRPhenomXError::PrecessionVersionUnavailable`].
    pub fn with_version_fallback(mut self, allow: bool) -> Self {
        self.version_fallback = allow;
        self
    }

    /// Whether the Euler-angle prescription may change automatically (see
    /// [`IMRPhenomXPHM::with_version_fallback`]).
    pub fn version_fallback(&self) -> bool {
        self.version_fallback
    }

    /// Without version fallbacks, the failure of the SpinTaylor angles as an
    /// error (otherwise unchanged, for [`fallback_to_msa`]).
    fn strict<T>(&self, result: Result<T, IMRPhenomXError>) -> Result<T, IMRPhenomXError> {
        if self.version_fallback {
            result
        } else {
            forbid_spintaylor_fallback(result)
        }
    }

    /// Without version fallbacks, an error unless `setup` uses the requested
    /// prescription.
    fn check(&self, setup: &Setup) -> Result<(), IMRPhenomXError> {
        if self.version_fallback {
            Ok(())
        } else {
            check_version(&setup.prec, self.version)
        }
    }

    fn build(
        binary: &PrecessingBinary,
        version: PrecessionVersion,
        modes: Vec<(u32, u32)>,
    ) -> Result<Self, IMRPhenomXError> {
        if binary.f_ref < 0.0 {
            return Err(IMRPhenomXError::NegativeReferenceFrequency);
        }
        if binary.luminosity_distance.is_nan() || binary.luminosity_distance <= 0.0 {
            return Err(IMRPhenomXError::NegativeDistance);
        }
        // Check the parameters (and precompute) with placeholder bounds,
        // which do not enter the parameter-only quantities.
        let f_placeholder = if binary.f_ref > 0.0 {
            binary.f_ref
        } else {
            20.0
        };
        let setup = Setup::new(
            binary,
            version,
            &modes,
            0.0,
            f_placeholder,
            2.0 * f_placeholder,
        )?;
        Ok(Self {
            binary: *binary,
            version,
            modes,
            custom_mode_array: false,
            version_fallback: true,
            setup: (binary.f_ref > 0.0).then_some(setup),
            tidal: None,
        })
    }

    /// The model with the NRTidal corrections of the (2,2) mode for the
    /// tidal deformabilities of the bodies of `binary.mass_1` and
    /// `binary.mass_2`, with only the (2,2) co-precessing mode: what
    /// `SimInspiralChooseFDWaveform` (`XLALSimIMRPhenomXPHM`) computes for
    /// IMRPhenomXP_NRTidalv2/v3 on uniform grids. (On frequency sequences
    /// LAL uses IMRPhenomXP, [`crate::imrphenomx::IMRPhenomXP::with_nrtidal`].)
    pub(crate) fn with_nrtidal(
        binary: &PrecessingBinary,
        version: PrecessionVersion,
        nrtidal: NRTidalVersion,
        lambda_1: f64,
        lambda_2: f64,
    ) -> Result<Self, IMRPhenomXError> {
        let tidal = Some((nrtidal, lambda_1, lambda_2));
        let modes = vec![(2, 2)];
        let mut model = Self::build(binary, version, modes.clone())?;
        // the parameter-only set-up with the tidal parameters (placeholder
        // bounds as in Self::build)
        let f_placeholder = if binary.f_ref > 0.0 {
            binary.f_ref
        } else {
            20.0
        };
        let setup = Setup::with_tidal(
            binary,
            version,
            tidal,
            &modes,
            0.0,
            f_placeholder,
            2.0 * f_placeholder,
        )?;
        model.setup = (binary.f_ref > 0.0).then_some(setup);
        model.custom_mode_array = true;
        model.tidal = tidal;
        Ok(model)
    }

    /// The binary parameters.
    pub fn binary(&self) -> &PrecessingBinary {
        &self.binary
    }

    /// The requested Euler angle prescription.
    pub fn version(&self) -> PrecessionVersion {
        self.version
    }

    /// The co-precessing modes `(l, m')`.
    pub fn modes(&self) -> &[(u32, u32)] {
        &self.modes
    }

    /// `(h_+, h_x)` at the given frequencies (Hz, positive and strictly
    /// increasing), reproducing `SimInspiralChooseFDWaveformSequence` with
    /// multibanding off: the first and last frequencies are LAL's `f_min`
    /// and `f_max`, frequencies above `min(f_max, f_cut)` with
    /// `f_cut = 0.3 / M` give zero.
    pub fn polarizations(&self, frequencies: &[f64]) -> Result<Polarizations, IMRPhenomXError> {
        let Some((f_min, f_max)) = frequency_bounds(frequencies)? else {
            return Ok((Vec::new(), Vec::new()));
        };
        let mut setup = self.setup(f_min, f_max)?;
        fallback_to_msa(
            self.strict(setup.generate(&self.modes, frequencies)),
            &self.binary,
            self.tidal,
            &self.modes,
            0.0,
            frequencies,
            f_min,
            f_max,
        )
    }

    /// [`Self::polarizations`] at the frequencies of a prepared
    /// [`FrequencyGrid`], written into `hplus` and `hcross` (each with the
    /// length of the grid).
    ///
    /// For evaluating many waveforms on the same frequencies: the Euler
    /// angles are those of [`Self::polarizations`] (evaluated once where
    /// two `m'` need them at the same argument), while the co-precessing
    /// modes use the powers of `Mf` stored in the grid and the twisting up
    /// is done in vectorisable passes over the grid. The result agrees with
    /// [`Self::polarizations`] to rounding (~1e-13 relative), including the
    /// fallback of the SpinTaylor angles to the MSA angles.
    ///
    /// # Panics
    ///
    /// If `hplus` or `hcross` do not have the length of the grid.
    pub fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<(), IMRPhenomXError> {
        assert_eq!(hplus.len(), grid.len(), "output length must match the grid");
        assert_eq!(
            hcross.len(),
            grid.len(),
            "output length must match the grid"
        );
        let Some((f_min, f_max)) = frequency_bounds(grid.frequencies())? else {
            return Ok(());
        };
        let modes = &self.modes;
        match self.strict(
            self.setup(f_min, f_max)?
                .generate_grid(modes, grid, hplus, hcross),
        ) {
            Err(IMRPhenomXError::SpinTaylorFallbackToMsa(_)) => Setup::with_tidal(
                &self.binary,
                PrecessionVersion::Msa,
                self.tidal,
                modes,
                0.0,
                f_min,
                f_max,
            )?
            .generate_grid(modes, grid, hplus, hcross),
            result => result,
        }
    }

    /// The set-up for a non-uniform grid from `f_min` to `f_max`.
    fn setup(&self, f_min: f64, f_max: f64) -> Result<Setup, IMRPhenomXError> {
        let setup = match &self.setup {
            Some(setup) => setup.with_bounds(f_min, f_max),
            None => Setup::with_tidal(
                &self.binary,
                self.version,
                self.tidal,
                &self.modes,
                0.0,
                f_min,
                f_max,
            )?,
        };
        self.check(&setup)?;
        Ok(setup)
    }

    /// `(h_+, h_x)` on the uniform grid `k delta_f`, `k = 0 .. n - 1`,
    /// reproducing `SimInspiralChooseFDWaveform(..., delta_f, f_min, f_max,
    /// f_ref, ...)` with multibanding off (`XLALSimIMRPhenomXPHM`): zero below
    /// `f_min` and above `min(f_max, f_cut)`, and `n - 1` the next power of
    /// two of `f_max / delta_f` (`f_max = 0` selects `f_cut = 0.3 / M`).
    pub fn polarizations_uniform(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<Polarizations, IMRPhenomXError> {
        if delta_f.is_nan() || delta_f <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "deltaF must be positive",
            ));
        }
        if f_min.is_nan() || f_min <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "f_min must be positive",
            ));
        }
        if f_max < 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "f_max must be non-negative",
            ));
        }
        let mut setup = Setup::with_tidal(
            &self.binary,
            self.version,
            self.tidal,
            &self.modes,
            delta_f,
            f_min,
            f_max,
        )?;
        self.check(&setup)?;
        let f_max_prime = setup.wf.f_max_prime;
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max_prime / delta_f) as usize + 1;
        let frequencies: Vec<f64> = (i_start..i_stop).map(|i| i as f64 * delta_f).collect();
        let (hp_band, hc_band) = fallback_to_msa(
            self.strict(setup.generate(&self.modes, &frequencies)),
            &self.binary,
            self.tidal,
            &self.modes,
            delta_f,
            &frequencies,
            f_min,
            f_max,
        )?;

        let last_freq = if f_max_prime < f_max {
            f_max
        } else {
            f_max_prime
        };
        let n_full =
            2.0f64.powf(((last_freq / delta_f) as usize as f64).log2().ceil()) as usize + 1;
        let mut hp = vec![Complex::new(0.0, 0.0); n_full.max(i_stop)];
        let mut hc = hp.clone();
        hp[i_start..i_stop].copy_from_slice(&hp_band);
        hc[i_start..i_stop].copy_from_slice(&hc_band);
        hp.truncate(n_full);
        hc.truncate(n_full);
        Ok((hp, hc))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imrphenomx::xp::IMRPhenomXP;

    fn binary() -> PrecessingBinary {
        PrecessingBinary {
            mass_1: 36.0,
            mass_2: 29.0,
            spin_1: [0.4, 0.2, 0.3],
            spin_2: [-0.2, 0.3, -0.1],
            luminosity_distance: 500.0,
            inclination: 0.7,
            phi_ref: 1.2,
            f_ref: 0.0,
        }
    }

    fn grid() -> Vec<f64> {
        (0..400).map(|i| 20.0 * 1.012f64.powi(i)).collect()
    }

    const VERSIONS: [PrecessionVersion; 2] =
        [PrecessionVersion::SpinTaylor, PrecessionVersion::Nnlo];

    #[test]
    fn default_reference_frequency_is_first_frequency() {
        let frequencies = grid();
        for version in VERSIONS {
            let default = IMRPhenomXPHM::with_version(&binary(), version)
                .unwrap()
                .polarizations(&frequencies)
                .unwrap();
            let explicit = IMRPhenomXPHM::with_version(
                &PrecessingBinary {
                    f_ref: frequencies[0],
                    ..binary()
                },
                version,
            )
            .unwrap()
            .polarizations(&frequencies)
            .unwrap();
            assert_eq!(default, explicit);
        }
    }

    /// The set-up precomputed with placeholder frequency bounds gives the
    /// same waveform as one built for the actual grid.
    #[test]
    fn precomputed_setup_does_not_depend_on_bounds() {
        let b = PrecessingBinary {
            f_ref: 30.0,
            ..binary()
        };
        let frequencies = grid();
        let modes = IMRPHENOMXHM_MODES.to_vec();
        for version in VERSIONS {
            let cached = IMRPhenomXPHM::with_version(&b, version)
                .unwrap()
                .polarizations(&frequencies)
                .unwrap();
            let (f_min, f_max) = (frequencies[0], frequencies[frequencies.len() - 1]);
            let fresh = Setup::new(&b, version, &modes, 0.0, f_min, f_max)
                .unwrap()
                .generate(&modes, &frequencies)
                .unwrap();
            assert_eq!(cached, fresh);
        }
    }

    #[test]
    fn body_labels_do_not_matter() {
        let frequencies = grid();
        let b = binary();
        let swapped = PrecessingBinary {
            mass_1: b.mass_2,
            mass_2: b.mass_1,
            spin_1: b.spin_2,
            spin_2: b.spin_1,
            ..b
        };
        for version in VERSIONS {
            let a = IMRPhenomXPHM::with_version(&b, version).unwrap();
            let s = IMRPhenomXPHM::with_version(&swapped, version).unwrap();
            assert_eq!(
                a.polarizations(&frequencies).unwrap(),
                s.polarizations(&frequencies).unwrap()
            );
        }
    }

    /// With only the (2,2) mode and NNLO angles IMRPhenomXPHM is IMRPhenomXP
    /// (LALSimulation gives identical results).
    #[test]
    fn dominant_mode_with_nnlo_angles_is_imrphenomxp() {
        let frequencies = grid();
        let (hp, hc) = IMRPhenomXPHM::with_version(&binary(), PrecessionVersion::Nnlo)
            .unwrap()
            .with_mode_array(&[(2, 2)])
            .unwrap()
            .polarizations(&frequencies)
            .unwrap();
        let (hp_xp, hc_xp) = IMRPhenomXP::with_version(&binary(), PrecessionVersion::Nnlo)
            .unwrap()
            .polarizations(&frequencies)
            .unwrap();
        let scale = hp_xp.iter().map(|h| h.norm()).fold(0.0, f64::max);
        for (a, b) in hp.iter().zip(&hp_xp).chain(hc.iter().zip(&hc_xp)) {
            assert!((a - b).norm() <= 1e-15 * scale, "{a} != {b}");
        }
    }

    #[test]
    fn mode_array_is_validated() {
        let model = IMRPhenomXPHM::new(&binary()).unwrap();
        assert_eq!(model.modes(), &IMRPHENOMXHM_MODES);
        assert_eq!(
            model.clone().with_mode_array(&[(4, 3)]).unwrap_err(),
            IMRPhenomXError::InvalidModeArray
        );
        assert_eq!(
            model.clone().with_mode_array(&[]).unwrap_err(),
            IMRPhenomXError::InvalidModeArray
        );
        let subset = model.with_mode_array(&[(4, 4), (2, 2)]).unwrap();
        assert_eq!(subset.modes(), &[(2, 2), (4, 4)]);
    }

    /// The NNLO and MSA angles of the grid path (with the angles shared
    /// between m') agree with those of the per-frequency loop.
    #[test]
    fn grid_angles_match_per_frequency_angles() {
        let grid = FrequencyGrid::uniform(20.0, 1024.0, 0.125);
        let frequencies = grid.frequencies();
        let n = grid.len();
        let modes = IMRPHENOMXHM_MODES.to_vec();
        for version in [PrecessionVersion::Nnlo, PrecessionVersion::Msa] {
            let setup = Setup::new(&binary(), version, &modes, 0.0, 20.0, 1024.0).unwrap();
            let (wf, prec) = (&setup.wf, &setup.prec);
            let angles: [Option<AngleArrays>; 5] = if version == PrecessionVersion::Msa {
                let argument = |mprime: u32, i: usize| {
                    PI * (wf.m_sec * frequencies[i]) * (2.0 / mprime as f64)
                };
                angles_with_reuse(&[2, 1, 3, 4], n, argument, |arguments, angles| {
                    prec.msa_angles_many(arguments, |alpha, epsilon, cos_beta| {
                        angles.push(alpha, cos_beta, epsilon)
                    });
                })
            } else {
                std::array::from_fn(|mprime| {
                    (mprime > 0)
                        .then(|| AngleArrays::nnlo(prec, wf.eta, &grid, wf.m_sec, mprime as u32, n))
                })
            };
            for mprime in 1..=4u32 {
                let angles = angles[mprime as usize].as_ref().unwrap();
                for (i, &f) in frequencies.iter().enumerate() {
                    let mf = wf.m_sec * f;
                    let (alpha, epsilon, cos_beta) = if version == PrecessionVersion::Msa {
                        prec.msa_angles((PI * mf * (2.0 / mprime as f64)).libm_cbrt())
                    } else {
                        prec.nnlo_angles(PI * mf * 2. / mprime as f64, wf.eta)
                    };
                    let close = |a: f64, b: f64| (a - b).abs() <= 1e-13 * b.abs().max(1.0);
                    assert!(
                        close(angles.alpha[i], alpha)
                            && close(angles.epsilon[i], epsilon)
                            && close(angles.cos_beta[i], cos_beta),
                        "{version:?} m' = {mprime}, f = {f}"
                    );
                }
            }
        }
    }

    /// Where the SpinTaylor integration does not reach the first frequency
    /// LALSimulation silently uses the MSA angles. With the default modes the
    /// integration starts at `(f_min - 0.5) / 2` and this binary has angles;
    /// with only the (2,2) mode it starts at `f_min - 0.5` and falls back.
    #[test]
    fn spintaylor_failure_falls_back_to_msa() {
        let heavy = PrecessingBinary {
            mass_1: 180.0,
            mass_2: 65.0,
            spin_1: [0.2, 0.1, -0.1],
            spin_2: [-0.1, -0.2, 0.1],
            luminosity_distance: 500.0,
            inclination: 0.4,
            phi_ref: 0.3,
            f_ref: 15.5,
        };
        let frequencies: Vec<f64> = (0..200).map(|i| 15.5 + i as f64).collect();
        let spintaylor = IMRPhenomXPHM::new(&heavy).unwrap();
        let msa = IMRPhenomXPHM::with_version(&heavy, PrecessionVersion::Msa).unwrap();
        assert_ne!(
            spintaylor.polarizations(&frequencies).unwrap(),
            msa.polarizations(&frequencies).unwrap()
        );
        let spintaylor = spintaylor.with_mode_array(&[(2, 2)]).unwrap();
        let msa = msa.with_mode_array(&[(2, 2)]).unwrap();
        assert_eq!(
            spintaylor.polarizations(&frequencies).unwrap(),
            msa.polarizations(&frequencies).unwrap()
        );
        assert_eq!(
            spintaylor.polarizations_uniform(0.25, 15.5, 200.0).unwrap(),
            msa.polarizations_uniform(0.25, 15.5, 200.0).unwrap()
        );
    }

    /// The vectorised merger-ringdown SpinTaylor angles of the grid path
    /// agree with those of the per-frequency path: gamma (integrated over
    /// ~20000 knots) and cos(beta).
    #[test]
    fn fast_spintaylor_angles_match() {
        let modes = IMRPHENOMXHM_MODES.to_vec();
        for b in [
            binary(),
            PrecessingBinary {
                mass_1: 40.0,
                mass_2: 20.0,
                spin_1: [0.4, -0.3, 0.2],
                spin_2: [0.1, 0.5, -0.3],
                f_ref: 20.0,
                ..binary()
            },
        ] {
            let setup = |fast| {
                let mut setup =
                    Setup::new(&b, PrecessionVersion::SpinTaylor, &modes, 0.0, 20.0, 1024.0)
                        .unwrap();
                let angles = setup.initialize_euler_angles(fast).unwrap();
                (setup, angles)
            };
            let (_, exact) = setup(false);
            let (_, fast) = setup(true);
            let (exact_gamma, fast_gamma) = (
                exact.gamma_spline.as_ref().unwrap(),
                fast.gamma_spline.as_ref().unwrap(),
            );
            assert_eq!(exact_gamma.x(), fast_gamma.x());
            let scale = exact_gamma
                .x()
                .iter()
                .map(|&mf| exact_gamma.eval(mf).abs())
                .fold(1.0, f64::max);
            for &mf in exact_gamma.x() {
                let difference = (exact_gamma.eval(mf) - fast_gamma.eval(mf)).abs();
                assert!(difference < 1e-13 * scale, "{mf}: {difference:e}");
            }
            assert!((exact.gamma_ref - fast.gamma_ref).abs() < 1e-13 * scale);
            let mf: Vec<f64> = (0..1000)
                .map(|i| exact.ftrans_mrd * (1.0 + 0.002 * i as f64))
                .collect();
            let mut cos_beta = vec![0.0; mf.len()];
            fast.mrd_cos_beta_fast(&mf, &mut cos_beta);
            for (&mf, &cos_beta) in mf.iter().zip(&cos_beta) {
                assert!((exact.beta_mrd.beta(mf).cos() - cos_beta).abs() < 1e-14);
            }
        }
    }
}
