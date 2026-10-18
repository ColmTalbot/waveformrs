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

use num_complex::Complex;

use crate::constants::{MSUN_SI, PC_SI};
use crate::imrphenomx::internals::{IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::libm::cexp_i;
use crate::imrphenomx::phase::PhaseCoefficients;
use crate::imrphenomx::pnr::{ringdown_frequency, PnrPrecession, PnrSplines};
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
#[derive(Clone, Debug)]
pub struct IMRPhenomXPNR {
    binary: PrecessingBinary,
    /// Co-precessing modes `(l, m' > 0)` in the order of LAL's mode loop.
    modes: Vec<(u32, u32)>,
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

/// Waveform and precession structs of one evaluation.
struct Setup {
    wf: WaveformStruct,
    prec: PnrPrecession,
}

impl Setup {
    /// `IMRPhenomXSetWaveformVariables` and
    /// `IMRPhenomXGetAndSetPrecessionVariables` as called by
    /// `XLALSimIMRPhenomXPHMFrequencySequence` (`delta_f = 0`) and
    /// `XLALSimIMRPhenomXPHM` (`delta_f > 0`), with the check of the
    /// reference frequency for the PNR angles.
    fn new(
        binary: &PrecessingBinary,
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
        // XPNR uses the original release of the higher modes
        wf.hm_release = 122019;
        // tuned angles: f_min <= f_ref <= f_max
        let upper = if delta_f > 0.0 { wf.f_max_prime } else { f_max };
        if !(wf.f_ref >= f_min && wf.f_ref <= upper) {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "f_min <= f_ref <= f_max is required when using tuned angles",
            ));
        }
        let prec = PnrPrecession::new(&mut wf, m1 * MSUN_SI, m2 * MSUN_SI, chi1, chi2, modes)?;
        Ok(Self { wf, prec })
    }

    /// `IMRPhenomXPHM_hplushcross` with the PNR flags on the grid
    /// `frequencies` (the band of a uniform grid), including the rotation of
    /// the polarisations by `2 zeta`.
    fn generate(
        &mut self,
        modes: &[(u32, u32)],
        frequencies: &[f64],
    ) -> Result<Polarizations, IMRPhenomXError> {
        if self.wf.f_max_prime <= self.wf.f_min {
            return Err(IMRPhenomXError::InvalidFrequencies("fCut <= f_min"));
        }
        if self.prec.precessing_tag == 3 {
            self.prec.initialize_euler_angles(&self.wf)?;
        }
        let tuned = self.prec.use_tuned_angles;
        let splines = if tuned {
            let splines = self.prec.angle_interpolants(&self.wf, modes)?;
            let f_ref = self.wf.f_ref;
            let (alpha_ref, beta_ref, gamma_ref) = splines.angles(f_ref)?;
            // the contribution of alpha0 is assigned in the remapping
            self.prec.alpha_offset = alpha_ref;
            // (the sign is flipped between gamma and epsilon)
            self.prec.epsilon_offset = -gamma_ref - self.prec.epsilon0;
            if !beta_ref.is_finite() {
                return Err(IMRPhenomXError::PnrAngles("invalid beta at f_ref"));
            }
            self.prec.remap_theta_jsf(&self.wf, beta_ref);
            Some(splines)
        } else {
            None
        };
        let wf = &self.wf;
        let prec = &self.prec;
        let phase22 = PhaseCoefficients::new(wf);

        let n = frequencies.len();
        let zero = Complex::new(0.0, 0.0);
        let mut hp = vec![zero; n];
        let mut hc = vec![zero; n];
        let mf_cut = wf.f_max_prime * wf.m_sec;
        let mf_rd_22 = wf.f_ring;

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
            // transition frequencies of the frequency map
            let (mf_low, mf_high, mf_rd_lm) = if tuned && (ell, emm) != (2, 2) {
                let mf_rd_lm = ringdown_frequency(ell, emm, wf);
                let (low, high) = linear_frequency_map_transition_frequencies(
                    emm,
                    mf_rd_22,
                    mf_rd_lm,
                    prec.pnr_hm_mflow,
                    prec.pnr_hm_mfhigh,
                );
                (low, high, mf_rd_lm)
            } else {
                (0.0, 0.0, 0.0)
            };
            let ylm = harmonics(prec, ell);
            for (i, &f) in frequencies.iter().enumerate() {
                let mf = wf.m_sec * f;
                // no waveform above Mf_max (default Mf = 0.3)
                if mf > mf_cut {
                    continue;
                }
                let h = generator.h(mf, &UsefulPowers::new(mf));
                let twist = match &splines {
                    Some(splines) => {
                        let mf_mapped = linear_frequency_map(
                            mf,
                            ell,
                            emm,
                            mf_low,
                            mf_high,
                            mf_rd_22,
                            mf_rd_lm,
                            prec.pnr_inspiral_scaling,
                        )
                        .ok_or(IMRPhenomXError::PnrAngles(
                            "invalid transition frequencies of the frequency map",
                        ))?;
                        let f_mapped = mf_to_hz(mf_mapped, wf.m_tot);
                        let (alpha, beta, gamma) = splines_at(splines, f_mapped)?;
                        let alpha = alpha - prec.alpha_offset;
                        let epsilon = -gamma - prec.epsilon_offset;
                        Twist::new(alpha, beta.cos(), epsilon, emm)
                    }
                    None => {
                        let omega = std::f64::consts::PI * mf * 2. / emm as f64;
                        let (alpha, cos_beta, epsilon) = prec.nnlo_angles(omega);
                        Twist::new(alpha, cos_beta, epsilon, emm)
                    }
                };
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

/// The PNR angles at `f` (Hz).
#[inline]
fn splines_at(splines: &PnrSplines, f: f64) -> Result<(f64, f64, f64), IMRPhenomXError> {
    splines.angles(f)
}

impl IMRPhenomXPNR {
    /// IMRPhenomXPNR with the default modes.
    ///
    /// Returns an error where LALSimulation would (masses, spins beyond the
    /// Kerr bound, non-positive distance, negative reference frequency,
    /// unphysical remnant fits).
    pub fn new(binary: &PrecessingBinary) -> Result<Self, IMRPhenomXError> {
        Self::build(binary, IMRPHENOMXHM_MODES.to_vec())
    }

    /// Restrict the waveform to the co-precessing modes `(l, m')` with
    /// `m' > 0` among (2,2), (2,1), (3,3), (3,2) and (4,4) (LAL's mode array
    /// with both signs of `m'`), which also changes the frequency range of
    /// the angles.
    pub fn with_mode_array(self, modes: &[(u32, u32)]) -> Result<Self, IMRPhenomXError> {
        if modes.is_empty() || modes.iter().any(|mode| !IMRPHENOMXHM_MODES.contains(mode)) {
            return Err(IMRPhenomXError::InvalidModeArray);
        }
        // LAL's mode loop order
        let ordered: Vec<(u32, u32)> = (2..=4u32)
            .flat_map(|ell| (1..=ell).map(move |emm| (ell, emm)))
            .filter(|mode| modes.contains(mode))
            .collect();
        Self::build(&self.binary, ordered)
    }

    fn build(binary: &PrecessingBinary, modes: Vec<(u32, u32)>) -> Result<Self, IMRPhenomXError> {
        if binary.f_ref < 0.0 {
            return Err(IMRPhenomXError::NegativeReferenceFrequency);
        }
        if binary.luminosity_distance.is_nan() || binary.luminosity_distance <= 0.0 {
            return Err(IMRPhenomXError::NegativeDistance);
        }
        Ok(Self {
            binary: *binary,
            modes,
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
        Setup::new(&self.binary, &self.modes, 0.0, f_min, f_max)?.generate(&self.modes, frequencies)
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
        let mut setup = Setup::new(&self.binary, &self.modes, delta_f, f_min, f_max)?;
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
