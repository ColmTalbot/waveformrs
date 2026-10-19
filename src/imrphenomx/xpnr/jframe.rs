//! The modes `h_lm` of IMRPhenomXPNR in the inertial J frame, as
//! `SimInspiralChooseFDModes` returns them: `XLALSimIMRPhenomXPHMOneMode`
//! (uniform grids) and `XLALSimIMRPhenomXPHMFrequencySequenceOneMode`
//! (frequency sequences) with the flags of IMRPhenomXPNR, through
//! `IMRPhenomXPHM_OneMode` and `IMRPhenomXPHMTwistUpOneMode`.
//!
//! For J-frame mode `(l, m)`, each co-precessing mode `h = h_{l,-m'}(f)`
//! (`m' > 0`, the same `l`) with antisymmetric part `h_as` (the (2,2) mode
//! only) contributes
//!
//! `h_lm(f) += e^{-i m' epsilon} (h + h_as) e^{-i m alpha} d^l_{m,-m'}(beta)`,
//! `h_lm(-f) += e^{i m' epsilon} (-1)^l conj(h - h_as) e^{-i m alpha} d^l_{m,m'}(beta)`,
//!
//! with the PNR angles at the frequency of the (2,2) mode given by the
//! frequency map (the NNLO angles at `Mf' = 2 Mf / m'` when LAL switches the
//! PNR angles off). As in LAL, without a mode array each J-frame mode `(l,
//! m)` has its own set-up, whose SpinTaylor angles cover `M_MIN = max(1,
//! |m|)` to `M_MAX = L_MAX_PNR = l`.

use num_complex::Complex;

use super::{check_prescription, Caller, IMRPhenomXPNR, ModeAngles, Setup};
use crate::constants::{MSUN_SI, PC_SI};
use crate::imrphenomx::internals::{IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::libm::cexp_i;
use crate::imrphenomx::phase::PhaseCoefficients;
use crate::imrphenomx::pnr::ModeRange;
use crate::imrphenomx::wigner::half_angles;
use crate::imrphenomx::xhm::modes::next_pow2;
use crate::imrphenomx::xhm::ModeGenerator;
use crate::imrphenomx::xp::frequency_bounds;
use crate::imrphenomx::xphm::jframe::EpsilonFactors;
use crate::imrphenomx::xphm::wigner_d;
use crate::model::{FrequencyMode, FrequencyModes};

/// One J-frame mode being accumulated.
struct Accumulator {
    ell: u32,
    m: i32,
    positive: Vec<Complex<f64>>,
    negative: Vec<Complex<f64>>,
}

impl IMRPhenomXPNR {
    /// `IMRPhenomXSetWaveformVariables` as called by the single-mode
    /// functions, with the bodies ordered by
    /// `XLALIMRPhenomXPCheckMassesAndSpins`.
    fn single_mode_struct(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<WaveformStruct, IMRPhenomXError> {
        let binary = &self.binary;
        let (m1, m2, chi1, chi2) = binary.ordered();
        WaveformStruct::from_si(
            m1 * MSUN_SI,
            m2 * MSUN_SI,
            chi1[2],
            chi2[2],
            delta_f,
            binary.f_ref,
            binary.phi_ref,
            f_min,
            f_max,
            binary.luminosity_distance * 1e6 * PC_SI,
            binary.inclination,
        )
    }

    /// The J-frame modes `(l, m)`, `l` among the co-precessing modes and
    /// `|m| <= l`, grouped by the mode range of their set-up.
    fn jframe_groups(&self) -> Vec<(ModeRange, Vec<(u32, i32)>)> {
        let mut ells: Vec<u32> = self.modes.iter().map(|&(ell, _)| ell).collect();
        ells.sort_unstable();
        ells.dedup();
        let mut groups: Vec<(ModeRange, Vec<(u32, i32)>)> = Vec::new();
        for ell in ells {
            for m in -(ell as i32)..=ell as i32 {
                let key = if self.custom_mode_array {
                    ModeRange::of_mode_array(&self.modes)
                } else {
                    ModeRange::of_jframe_mode(ell, m)
                };
                match groups.iter_mut().find(|(k, _)| *k == key) {
                    Some((_, modes)) => modes.push((ell, m)),
                    None => groups.push((key, vec![(ell, m)])),
                }
            }
        }
        groups
    }

    /// The modes `h_lm` in the J frame at `frequencies` (Hz, positive and
    /// increasing), as `XLALSimIMRPhenomXPHMFrequencySequenceOneMode`
    /// computes each of them with the flags of IMRPhenomXPNR, for every `l`
    /// of the co-precessing modes and `|m| <= l`.
    pub(crate) fn jframe_modes(
        &self,
        frequencies: &[f64],
    ) -> Result<FrequencyModes, IMRPhenomXError> {
        let Some((f_min, f_max)) = frequency_bounds(frequencies)? else {
            return Ok(FrequencyModes::default());
        };
        self.evaluate_jframe(0.0, f_min, f_max, frequencies, 0, frequencies.len())
    }

    /// The modes `h_lm` in the J frame on the uniform grid `k delta_f`, `k
    /// = 0 .. n - 1`, as `SimInspiralChooseFDModes`
    /// (`XLALSimIMRPhenomXPHMModes` with the flags of IMRPhenomXPNR): zero
    /// below `f_min` and above `min(f_max, f_cut)`, `n - 1` the next power
    /// of two of `max(f_max, f_cut) / delta_f` (`f_max = 0` selects
    /// `f_cut`).
    pub(crate) fn jframe_modes_uniform(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<FrequencyModes, IMRPhenomXError> {
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
        if f_max.is_nan() || f_max < 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "f_max must be non-negative",
            ));
        }
        // fCut (and hence f_max_prime) does not depend on the precession
        let f_max_prime = self.single_mode_struct(delta_f, f_min, f_max)?.f_max_prime;
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max_prime / delta_f) as usize + 1;
        let last_freq = if f_max_prime < f_max {
            f_max
        } else {
            f_max_prime
        };
        let n_full = next_pow2((last_freq / delta_f) as usize) + 1;
        let n = n_full.max(i_stop);
        let frequencies: Vec<f64> = (0..n).map(|i| i as f64 * delta_f).collect();
        let mut modes =
            self.evaluate_jframe(delta_f, f_min, f_max, &frequencies, i_start, i_stop)?;
        modes.frequencies.truncate(n_full);
        for mode in &mut modes.modes {
            mode.positive.truncate(n_full);
            mode.negative.truncate(n_full);
        }
        Ok(modes)
    }

    /// The J-frame modes on `frequencies`, evaluated at the indices
    /// `start..stop` (zero elsewhere), for LAL's grid bounds `f_min`,
    /// `f_max`.
    fn evaluate_jframe(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
        frequencies: &[f64],
        start: usize,
        stop: usize,
    ) -> Result<FrequencyModes, IMRPhenomXError> {
        let n = frequencies.len();
        let zero = Complex::new(0.0, 0.0);
        let band = &frequencies[start..stop];
        let mut out: Vec<Accumulator> = Vec::new();
        let mut frame = None;
        for (range, jframe) in self.jframe_groups() {
            let mut accumulators: Vec<Accumulator> = jframe
                .iter()
                .map(|&(ell, m)| Accumulator {
                    ell,
                    m,
                    positive: vec![zero; n],
                    negative: vec![zero; n],
                })
                .collect();
            let mut setup = Setup::new(
                &self.binary,
                range,
                delta_f,
                f_min,
                f_max,
                Caller::SingleMode,
                self.parity_fix,
            )?;
            if !self.version_fallback {
                check_prescription(&setup.prec)?;
            }
            twist_up_modes(&mut setup, &self.modes, band, start, &mut accumulators)?;
            frame.get_or_insert((setup.prec.theta_jn, setup.prec.zeta_polarization));
            out.extend(accumulators);
        }
        out.sort_by_key(|a| (a.ell, a.m));
        let (theta, zeta) = frame.unwrap_or((0.0, 0.0));
        Ok(FrequencyModes {
            frequencies: frequencies.to_vec(),
            theta,
            phi: 0.0,
            polarization_angle: zeta,
            modes: out
                .into_iter()
                .map(|a| FrequencyMode {
                    ell: a.ell,
                    m: a.m,
                    positive: a.positive,
                    negative: a.negative,
                })
                .collect(),
        })
    }
}

/// `IMRPhenomXPHM_OneMode` with the PNR flags for the J-frame modes of
/// `accumulators` (all sharing `setup`): twist up the co-precessing modes
/// of the same `l` at `band`, the frequencies of indices `offset..` of the
/// accumulators.
fn twist_up_modes(
    setup: &mut Setup,
    coprecessing: &[(u32, u32)],
    band: &[f64],
    offset: usize,
    accumulators: &mut [Accumulator],
) -> Result<(), IMRPhenomXError> {
    let splines = setup.prepare_angles(coprecessing)?;
    let tuned = splines.is_some();
    let wf = &setup.wf;
    let prec = &setup.prec;
    let phase22 = PhaseCoefficients::new(wf);
    let mf_cut = wf.f_max_prime * wf.m_sec;
    let zero = Complex::new(0.0, 0.0);
    for &(ell, emm) in coprecessing {
        if !accumulators.iter().any(|a| a.ell == ell) {
            continue;
        }
        // skip twisting up when the co-precessing mode vanishes
        if wf.q == 1.0 && wf.chi1l == wf.chi2l && emm % 2 != 0 {
            continue;
        }
        let generator = ModeGenerator::new(wf, &phase22, ell, emm);
        let antisymmetric = if (ell, emm) == (2, 2) && prec.antisymmetric && tuned {
            Some(prec.antisymmetric_waveform(wf, band)?)
        } else {
            None
        };
        let angles = ModeAngles::new(wf, prec, splines.as_ref(), ell, emm);
        let l = ell as usize;
        let sigma = if ell % 2 == 0 { 1.0 } else { -1.0 };
        for (i, &f) in band.iter().enumerate() {
            let mf = wf.m_sec * f;
            if mf > mf_cut {
                continue;
            }
            let h = generator.h(mf, &UsefulPowers::new(mf));
            let h_as = match &antisymmetric {
                Some((amp, phase)) => amp[i] * cexp_i(phase[i]),
                None => zero,
            };
            let (alpha, cos_beta, epsilon) = angles.at(mf)?;
            let (cb, sb) = half_angles(cos_beta);
            let mut c = [1.0; 9];
            let mut s = [1.0; 9];
            for k in 1..9 {
                c[k] = cb * c[k - 1];
                s[k] = sb * s[k - 1];
            }
            let d = wigner_d(ell, emm, &c, &s);
            let e = cexp_i(emm as f64 * epsilon);
            let factors = EpsilonFactors::with_antisymmetric(h, h_as, e, sigma);
            for accumulator in accumulators.iter_mut().filter(|a| a.ell == ell) {
                let m = accumulator.m;
                let rotation = cexp_i(-(m as f64) * alpha);
                let (positive, negative) = factors.contribution(&d, l, emm, m, rotation);
                accumulator.positive[offset + i] += positive;
                accumulator.negative[offset + i] += negative;
            }
        }
    }
    Ok(())
}
