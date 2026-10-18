//! The precessing higher-mode model IMRPhenomXPNR of lalsimulation 6.2.0.
//!
//! In LALSimulation, `IMRPhenomXPNR` is IMRPhenomXPHM (`XLALSimIMRPhenomXPHM`,
//! `XLALSimIMRPhenomXPHMFrequencySequence`) with the LALDict flags
//! (unless set by the user)
//!
//! - `PhenomXHMReleaseVersion = 122019` (always overwritten),
//! - `PNRUseTunedAngles = 1`: the PNR precession angles (Hamilton et al.
//!   2021, arXiv:2107.08876) with the frequency map of the higher modes
//!   (arXiv:2312.10025),
//! - `PNRUseTunedCoprec = 1`: the PhenomXCP deviations of the (2,2)
//!   co-precessing mode and the effective ringdown frequency of the higher
//!   modes (`PNRUseTunedCoprec33` stays 0),
//! - `PNRForceXHMAlignment = 0`,
//! - `AntisymmetricWaveform = 1`: the antisymmetric part of the (2,|2|)
//!   co-precessing modes (Ghosh et al. 2024, arXiv:2310.16980),
//! - `PrecVersion = 330`: SpinTaylor angles connected to the PNR angles,
//! - `FinalSpinMod = 7`: the precessing final spin of the SpinTaylor angles
//!   with the sign of `cos(beta)` in the ringdown,
//! - `MBandPrecVersion = 2` (multibanding of the angles disabled for
//!   `q > 7` with large opening angles; not used here since multibanding is
//!   off).
//!
//! The co-precessing modes are those of IMRPhenomXPHM, (2,2), (2,1),
//! (3,3), (3,2) and (4,4), with both signs of `m'`.
//!
//! LAL switches the PNR machinery off by itself in parts of parameter space:
//!
//! - in-plane spins below `1e-6` use the NNLO angles (102) inside the PNR
//!   machinery, below `1e-7` the PNR angles, the co-precessing deviations
//!   and the antisymmetric waveform are off (plain XPHM with NNLO angles,
//!   the 122019 release and the final spin of `FinalSpinMod` 7);
//! - when the SpinTaylor integration fails the MSA angles (223) are used
//!   (and the NNLO angles if the MSA system fails);
//! - outside `q <= 8.5`, `chi <= 0.85` (effective single spin) the PNR
//!   angles are blended with the PN angles, and beyond `q = 12` or
//!   `chi = 1.2` the PN angles are used; the higher modes then use the
//!   inspiral scaling `2 f / m'` of the frequency map;
//! - the co-precessing deviations are switched off between `q = 10` and
//!   `q = 20`.
//!
//! With [`IMRPhenomXPNR::with_version_fallback`]`(false)` the model fails
//! instead of leaving the PNR angles with the SpinTaylor angles (330).

mod jframe;

use num_complex::Complex;

use crate::constants::{MSUN_SI, PC_SI};
use crate::imrphenomx::internals::{IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::libm::cexp_i;
use crate::imrphenomx::phase::PhaseCoefficients;
use crate::imrphenomx::pnr::{ringdown_frequency, ModeRange, PnrPrecession, PnrSplines};
use crate::imrphenomx::pnr_angles::{
    linear_frequency_map, linear_frequency_map_transition_frequencies,
};
use crate::imrphenomx::precession::{mf_to_hz, PrecessingBinary};
use crate::imrphenomx::wigner::half_angles;
use crate::imrphenomx::xhm::{ModeGenerator, IMRPHENOMXHM_MODES};
use crate::imrphenomx::xp::{frequency_bounds, Polarizations};
use crate::imrphenomx::xphm::wigner_d;

/// IMRPhenomXPNR (lalsimulation 6.2.0): IMRPhenomXPHM with the PNR angles,
/// the tuned co-precessing (2,2) mode, the antisymmetric (2,|2|) modes and
/// the 122019 release of the higher modes (see the module documentation),
/// without multibanding.
///
/// # Conventions
///
/// As for [`crate::imrphenomx::IMRPhenomXPHM`]:
/// [`IMRPhenomXPNR::polarizations`] reproduces
/// `SimInspiralChooseFDWaveformSequence(phi_ref, m1, m2, *spin_1, *spin_2,
/// f_ref, distance, inclination, None, IMRPhenomXPNR, frequencies)` and
/// [`IMRPhenomXPNR::polarizations_uniform`] `SimInspiralChooseFDWaveform`
/// with `deltaF > 0` and multibanding switched off
/// (`PhenomXHMThresholdMband = PhenomXPHMThresholdMband = 0`; LAL's default
/// multibanding of the uniform grid is not reproduced). Masses are in solar
/// masses, the distance in Mpc and the spins are the dimensionless spin
/// vectors in the L frame at the reference frequency.
///
/// The PNR angles require the reference frequency to lie between the first
/// and last frequencies (for a uniform grid: between `f_min` and the
/// cutoff), as in LAL.
///
/// The J-frame modes are available through
/// [`crate::FrequencyDomainModel::modes`] (the single-mode functions of
/// LAL with the flags of IMRPhenomXPNR) and [`crate::choose_fd_modes`]
/// (`SimInspiralChooseFDModes`).
///
/// ```
/// use waveformrs::imrphenomx::{IMRPhenomXPNR, PrecessingBinary};
///
/// let binary = PrecessingBinary {
///     mass_1: 36.0,
///     mass_2: 29.0,
///     spin_1: [0.4, 0.2, 0.3],
///     spin_2: [-0.2, 0.3, -0.1],
///     luminosity_distance: 500.0,
///     inclination: 0.7,
///     phi_ref: 1.2,
///     f_ref: 20.0,
/// };
/// let frequencies: Vec<f64> = (0..200).map(|i| 20.0 + 2.0 * i as f64).collect();
/// let (hplus, hcross) = IMRPhenomXPNR::new(&binary)?.polarizations(&frequencies)?;
/// # assert_eq!(hplus.len(), 200);
/// # assert_eq!(hcross.len(), 200);
/// # Ok::<(), waveformrs::imrphenomx::IMRPhenomXError>(())
/// ```
#[derive(Clone, Debug)]
pub struct IMRPhenomXPNR {
    binary: PrecessingBinary,
    /// Co-precessing modes `(l, m' > 0)` in the order of LAL's mode loop.
    modes: Vec<(u32, u32)>,
    /// Whether the modes were restricted with a mode array (which then
    /// sets the frequency range of the angles of the single-mode
    /// functions).
    custom_mode_array: bool,
    /// Whether the precession prescription may leave the PNR angles with
    /// the SpinTaylor angles as in LAL.
    version_fallback: bool,
}

/// The co-precessing modes `(l, m' > 0)` among `modes` in the order of
/// LAL's mode loop (`l = 2..4`, `m' = 1..l`).
fn lal_order(modes: &[(u32, u32)]) -> Vec<(u32, u32)> {
    (2..=4u32)
        .flat_map(|ell| (1..=ell).map(move |emm| (ell, emm)))
        .filter(|mode| modes.contains(mode))
        .collect()
}

/// An error unless the PNR angles are used with the SpinTaylor angles
/// (version 330), for models without version fallbacks.
fn check_prescription(prec: &PnrPrecession) -> Result<(), IMRPhenomXError> {
    if prec.use_tuned_angles && prec.version == 330 {
        return Ok(());
    }
    let reason = if !prec.use_tuned_angles {
        "the in-plane spins are below 1e-7, which switches the PNR angles off"
    } else if prec.version == 223 {
        "the SpinTaylor angles cannot be computed"
    } else {
        "the in-plane spins are below 1e-6 or the MSA system cannot be initialised"
    };
    Err(IMRPhenomXError::PrecessionVersionUnavailable {
        requested: 330,
        fallback: prec.version,
        reason,
    })
}

/// The Euler angles `(alpha, cos(beta), epsilon)` of one co-precessing mode
/// `(l, m')`: the PNR angles at the frequency of the (2,2) mode given by
/// the frequency map, or the NNLO angles without the PNR angles.
struct ModeAngles<'a> {
    wf: &'a WaveformStruct,
    prec: &'a PnrPrecession,
    splines: Option<&'a PnrSplines>,
    ell: u32,
    emm: u32,
    /// Transition frequencies of the map and the (l, m') ringdown frequency.
    mf_low: f64,
    mf_high: f64,
    mf_rd_lm: f64,
    /// The interval of the splines at the previous frequency.
    hint: std::cell::Cell<usize>,
}

impl<'a> ModeAngles<'a> {
    fn new(
        wf: &'a WaveformStruct,
        prec: &'a PnrPrecession,
        splines: Option<&'a PnrSplines>,
        ell: u32,
        emm: u32,
    ) -> Self {
        let (mf_low, mf_high, mf_rd_lm) = if splines.is_some() && (ell, emm) != (2, 2) {
            let mf_rd_lm = ringdown_frequency(ell, emm, wf);
            let (low, high) = linear_frequency_map_transition_frequencies(
                emm,
                wf.f_ring,
                mf_rd_lm,
                prec.pnr_hm_mflow,
                prec.pnr_hm_mfhigh,
            );
            (low, high, mf_rd_lm)
        } else {
            (0.0, 0.0, 0.0)
        };
        Self {
            wf,
            prec,
            splines,
            ell,
            emm,
            mf_low,
            mf_high,
            mf_rd_lm,
            hint: std::cell::Cell::new(0),
        }
    }

    /// The angles at `Mf`, with the offsets at the reference frequency.
    #[inline]
    fn at(&self, mf: f64) -> Result<(f64, f64, f64), IMRPhenomXError> {
        let (wf, prec) = (self.wf, self.prec);
        match self.splines {
            Some(splines) => {
                let mf_mapped = linear_frequency_map(
                    mf,
                    self.ell,
                    self.emm,
                    self.mf_low,
                    self.mf_high,
                    wf.f_ring,
                    self.mf_rd_lm,
                    prec.pnr_inspiral_scaling,
                )
                .ok_or(IMRPhenomXError::PnrAngles(
                    "invalid transition frequencies of the frequency map",
                ))?;
                let f_mapped = mf_to_hz(mf_mapped, wf.m_tot);
                let mut hint = self.hint.get();
                let (alpha, beta, gamma) = splines.angles_from(f_mapped, &mut hint)?;
                self.hint.set(hint);
                Ok((
                    alpha - prec.alpha_offset,
                    beta.cos(),
                    -gamma - prec.epsilon_offset,
                ))
            }
            None => {
                let omega = std::f64::consts::PI * mf * 2. / self.emm as f64;
                Ok(prec.nnlo_angles(omega))
            }
        }
    }
}

/// The twisting-up factors shared by all modes with the same `m'` at one
/// frequency (as in `IMRPhenomXPHMTwistUp`).
#[derive(Clone, Copy, Debug, Default)]
struct Twist {
    cos_k: [f64; 5],
    sin_k: [f64; 5],
    c: [f64; 9],
    s: [f64; 9],
    /// `exp(-i m' epsilon) / 2`.
    eps: Complex<f64>,
}

impl Twist {
    fn new(alpha: f64, cos_beta: f64, epsilon: f64, mprime: u32) -> Self {
        // powers of exp(i alpha) by multiplication, as in LAL
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

/// `IMRPhenomXPHMTwistUp` for the co-precessing mode `h = h_{l,-m'}` with
/// `PolarizationSymmetry` `symmetry` (+1 for the symmetric modes, -1 for
/// the antisymmetric (2,2) mode): the contributions `(h_+, h_x)` in the J
/// frame.
#[inline]
fn twist_up(
    ell: u32,
    mprime: u32,
    twist: &Twist,
    ylm: &[f64],
    h: Complex<f64>,
    symmetry: f64,
) -> (Complex<f64>, Complex<f64>) {
    let l = ell as usize;
    let d = wigner_d(ell, mprime, &twist.c, &twist.s);
    let sigma = if ell % 2 == 0 { symmetry } else { -symmetry };
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

/// The `s = -2` spherical harmonics at `(thetaJN, 0)` for mode `l`.
fn harmonics(prec: &PnrPrecession, ell: u32) -> &[f64] {
    match ell {
        2 => &prec.ylm.y2,
        3 => &prec.ylm.y3,
        _ => &prec.ylm.y4,
    }
}

/// The LAL function that sets up the waveform and precession structs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Caller {
    /// `XLALSimIMRPhenomXPHMFrequencySequence` (`delta_f = 0`) and
    /// `XLALSimIMRPhenomXPHM` (`delta_f > 0`, which keeps the input order
    /// of the bodies in the waveform struct), which check that the
    /// reference frequency is inside the grid.
    Polarizations,
    /// The single-mode functions `XLALSimIMRPhenomXPHMOneMode` and
    /// `XLALSimIMRPhenomXPHMFrequencySequenceOneMode`.
    SingleMode,
}

/// Waveform and precession structs of one evaluation.
struct Setup {
    wf: WaveformStruct,
    prec: PnrPrecession,
}

impl Setup {
    /// `IMRPhenomXSetWaveformVariables` and
    /// `IMRPhenomXGetAndSetPrecessionVariables` as called by `caller`, with
    /// the angles covering the modes of `range`.
    fn new(
        binary: &PrecessingBinary,
        range: ModeRange,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
        caller: Caller,
    ) -> Result<Self, IMRPhenomXError> {
        let (m1, m2, chi1, chi2) = binary.ordered();
        let (wf_m1, wf_m2, wf_chi1, wf_chi2) = if delta_f > 0.0 && caller == Caller::Polarizations {
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
        // XPNR uses the original release of the higher modes
        wf.hm_release = 122019;
        // tuned angles: f_min <= f_ref <= f_max
        let upper = if delta_f > 0.0 { wf.f_max_prime } else { f_max };
        if caller == Caller::Polarizations && !(wf.f_ref >= f_min && wf.f_ref <= upper) {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "f_min <= f_ref <= f_max is required when using tuned angles",
            ));
        }
        let prec = PnrPrecession::new(&mut wf, m1 * MSUN_SI, m2 * MSUN_SI, chi1, chi2, range)?;
        Ok(Self { wf, prec })
    }

    /// The SpinTaylor angles (`IMRPhenomX_Initialize_Euler_Angles`) and the
    /// interpolants of the PNR angles for the co-precessing modes `modes`,
    /// with the offsets at the reference frequency and the J frame of the
    /// line of sight remapped with the PNR `beta`; `None` without the PNR
    /// angles.
    fn prepare_angles(
        &mut self,
        modes: &[(u32, u32)],
    ) -> Result<Option<PnrSplines>, IMRPhenomXError> {
        if self.wf.f_max_prime <= self.wf.f_min {
            return Err(IMRPhenomXError::InvalidFrequencies("fCut <= f_min"));
        }
        if self.prec.precessing_tag == 3 {
            self.prec.initialize_euler_angles(&self.wf)?;
        }
        if !self.prec.use_tuned_angles {
            return Ok(None);
        }
        let splines = self.prec.angle_interpolants(&self.wf, modes)?;
        let (alpha_ref, beta_ref, gamma_ref) = splines.angles(self.wf.f_ref)?;
        // the contribution of alpha0 is assigned in the remapping
        self.prec.alpha_offset = alpha_ref;
        // (the sign is flipped between gamma and epsilon)
        self.prec.epsilon_offset = -gamma_ref - self.prec.epsilon0;
        if !beta_ref.is_finite() {
            return Err(IMRPhenomXError::PnrAngles("invalid beta at f_ref"));
        }
        self.prec.remap_theta_jsf(&self.wf, beta_ref);
        Ok(Some(splines))
    }

    /// `IMRPhenomXPHM_hplushcross` with the PNR flags on the grid
    /// `frequencies` (the band of a uniform grid), including the rotation of
    /// the polarisations by `2 zeta`.
    fn generate(
        &mut self,
        modes: &[(u32, u32)],
        frequencies: &[f64],
    ) -> Result<Polarizations, IMRPhenomXError> {
        let splines = self.prepare_angles(modes)?;
        let tuned = splines.is_some();
        let wf = &self.wf;
        let prec = &self.prec;
        let phase22 = PhaseCoefficients::new(wf);

        let n = frequencies.len();
        let zero = Complex::new(0.0, 0.0);
        let mut hp = vec![zero; n];
        let mut hc = vec![zero; n];
        let mf_cut = wf.f_max_prime * wf.m_sec;

        for &(ell, emm) in modes {
            // skip twisting up when the co-precessing mode vanishes
            if wf.q == 1.0 && wf.chi1l == wf.chi2l && emm % 2 != 0 {
                continue;
            }
            let generator = ModeGenerator::new(wf, &phase22, ell, emm);
            let antisymmetric = if (ell, emm) == (2, 2) && prec.antisymmetric && tuned {
                Some(prec.antisymmetric_waveform(wf, frequencies)?)
            } else {
                None
            };
            let angles = ModeAngles::new(wf, prec, splines.as_ref(), ell, emm);
            let ylm = harmonics(prec, ell);
            for (i, &f) in frequencies.iter().enumerate() {
                let mf = wf.m_sec * f;
                // no waveform above Mf_max (default Mf = 0.3)
                if mf > mf_cut {
                    continue;
                }
                let h = generator.h(mf, &UsefulPowers::new(mf));
                let (alpha, cos_beta, epsilon) = angles.at(mf)?;
                let twist = Twist::new(alpha, cos_beta, epsilon, emm);
                let (mut p, mut c) = twist_up(ell, emm, &twist, ylm, h, 1.0);
                if let Some((amp, phase)) = &antisymmetric {
                    let h_antisymmetric = amp[i] * cexp_i(phase[i]);
                    let (pa, ca) = twist_up(ell, emm, &twist, ylm, h_antisymmetric, -1.0);
                    p += pa;
                    c += ca;
                }
                hp[i] += p;
                hc[i] += c;
            }
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
}

impl IMRPhenomXPNR {
    /// IMRPhenomXPNR with the default modes.
    ///
    /// Returns an error where LALSimulation would (masses, spins beyond the
    /// Kerr bound, non-positive distance, negative reference frequency,
    /// unphysical remnant fits).
    pub fn new(binary: &PrecessingBinary) -> Result<Self, IMRPhenomXError> {
        Self::build(binary, lal_order(&IMRPHENOMXHM_MODES), false)
    }

    /// Restrict the waveform to the co-precessing modes `(l, m')` with
    /// `m' > 0` among (2,2), (2,1), (3,3), (3,2) and (4,4) (LAL's mode array
    /// with both signs of `m'`), which also changes the frequency range of
    /// the angles.
    pub fn with_mode_array(self, modes: &[(u32, u32)]) -> Result<Self, IMRPhenomXError> {
        if modes.is_empty() || modes.iter().any(|mode| !IMRPHENOMXHM_MODES.contains(mode)) {
            return Err(IMRPhenomXError::InvalidModeArray);
        }
        Ok(Self {
            modes: lal_order(modes),
            custom_mode_array: true,
            ..self
        })
    }

    /// Whether the Euler angles may leave the PNR angles with the
    /// SpinTaylor angles (version 330) as LALSimulation does (`true`, the
    /// default): the NNLO angles inside the PNR machinery for in-plane spins
    /// below `1e-6`, plain IMRPhenomXPHM with the NNLO angles below `1e-7`,
    /// the MSA angles when the SpinTaylor integration fails and the NNLO
    /// angles when the MSA system cannot be initialised. With `false` these
    /// are [`IMRPhenomXError::PrecessionVersionUnavailable`] errors.
    pub fn with_version_fallback(self, allow: bool) -> Self {
        Self {
            version_fallback: allow,
            ..self
        }
    }

    /// Whether the precession prescription may change as in LALSimulation
    /// (see [`IMRPhenomXPNR::with_version_fallback`]).
    pub fn version_fallback(&self) -> bool {
        self.version_fallback
    }

    /// The waveform and precession structs of the polarizations.
    fn polarization_setup(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<Setup, IMRPhenomXError> {
        let setup = Setup::new(
            &self.binary,
            ModeRange::of_mode_array(&self.modes),
            delta_f,
            f_min,
            f_max,
            Caller::Polarizations,
        )?;
        if !self.version_fallback {
            check_prescription(&setup.prec)?;
        }
        Ok(setup)
    }

    fn build(
        binary: &PrecessingBinary,
        modes: Vec<(u32, u32)>,
        custom_mode_array: bool,
    ) -> Result<Self, IMRPhenomXError> {
        if binary.f_ref < 0.0 {
            return Err(IMRPhenomXError::NegativeReferenceFrequency);
        }
        if binary.luminosity_distance.is_nan() || binary.luminosity_distance <= 0.0 {
            return Err(IMRPhenomXError::NegativeDistance);
        }
        Ok(Self {
            binary: *binary,
            modes,
            custom_mode_array,
            version_fallback: true,
        })
    }

    /// The binary parameters.
    pub fn binary(&self) -> &PrecessingBinary {
        &self.binary
    }

    /// The co-precessing modes `(l, m')`.
    pub fn modes(&self) -> &[(u32, u32)] {
        &self.modes
    }

    /// `(h_+, h_x)` at the given frequencies (Hz, positive and strictly
    /// increasing), reproducing `SimInspiralChooseFDWaveformSequence`: the
    /// first and last frequencies are LAL's `f_min` and `f_max`, frequencies
    /// above `min(f_max, f_cut)` with `f_cut = 0.3 / M` give zero.
    pub fn polarizations(&self, frequencies: &[f64]) -> Result<Polarizations, IMRPhenomXError> {
        let Some((f_min, f_max)) = frequency_bounds(frequencies)? else {
            return Ok((Vec::new(), Vec::new()));
        };
        self.polarization_setup(0.0, f_min, f_max)?
            .generate(&self.modes, frequencies)
    }

    /// `(h_+, h_x)` on the uniform grid `k delta_f`, `k = 0 .. n - 1`,
    /// reproducing `SimInspiralChooseFDWaveform(..., delta_f, f_min, f_max,
    /// f_ref, ...)` with multibanding off: zero below `f_min` and above
    /// `min(f_max, f_cut)`, and `n - 1` the next power of two of
    /// `f_max / delta_f` (`f_max = 0` selects `f_cut = 0.3 / M`).
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
        let mut setup = self.polarization_setup(delta_f, f_min, f_max)?;
        let f_max_prime = setup.wf.f_max_prime;
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max_prime / delta_f) as usize + 1;
        let frequencies: Vec<f64> = (i_start..i_stop).map(|i| i as f64 * delta_f).collect();
        let (hp_band, hc_band) = setup.generate(&self.modes, &frequencies)?;
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
