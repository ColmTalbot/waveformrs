// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The modes `h_lm` of IMRPhenomXPHM in the inertial J frame, as
//! `SimInspiralChooseFDModes` returns them: `XLALSimIMRPhenomXPHMOneMode`
//! (uniform grids) and `XLALSimIMRPhenomXPHMFrequencySequenceOneMode`
//! (frequency sequences), with `IMRPhenomXPHM_OneMode` and
//! `IMRPhenomXPHMTwistUpOneMode` of `LALSimIMRPhenomXPHM.c`.
//!
//! For J-frame mode `(l, m)`, each co-precessing mode `h = h_{l,-m'}(f)`
//! (`m' > 0`, the same `l`) contributes
//!
//! `h_lm(f) += e^{-i m' epsilon} h e^{-i m alpha} d^l_{m,-m'}(beta)`,
//! `h_lm(-f) += e^{i m' epsilon} (-1)^l conj(h) e^{-i m alpha} d^l_{m,m'}(beta)`,
//!
//! with the Euler angles at `Mf' = 2 Mf / m'`. As in LAL, without a mode
//! array in the LALDict the SpinTaylor angles of mode `(l, m)` are
//! integrated over the frequency range of `M_MIN = max(1, |m|)` and
//! `M_MAX = l`, so each `(l, |m|)` has its own set-up (and fallback to the
//! MSA angles).

use std::f64::consts::PI;

use num_complex::Complex;

use super::{build_modes, spintaylor_angles, wigner_d, AngleCursor, IMRPhenomXPHM, Setup};
use crate::complex::cexp_i;
use crate::constants::{MSUN_SI, PC_SI};
use crate::imrphenomx::internals::{IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::precession::{PrecessingBinary, PrecessionStruct, PrecessionVersion};
use crate::imrphenomx::version_fallback::{check_version, forbid_spintaylor_fallback};
use crate::imrphenomx::wigner::half_angles;
use crate::imrphenomx::xhm::modes::next_pow2;
use crate::imrphenomx::xp::frequency_bounds;
use crate::model::{FrequencyMode, FrequencyModes};

/// How the frequency range of the SpinTaylor angles (`M_MIN`, `M_MAX` of
/// LAL's precession struct) is chosen.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AngleRange {
    /// LAL's single-mode functions without a `ModeArray` in the LALDict:
    /// `M_MIN = max(1, |m|)` and `M_MAX = l` of the J-frame mode.
    PerMode,
    /// The range of `|m'|` of the co-precessing mode array (a `ModeArray`
    /// in the LALDict).
    ModeArray,
}

/// The J-frame modes of a precessing binary.
#[derive(Clone, Copy, Debug)]
pub(crate) struct JFrameModes<'a> {
    pub(crate) binary: &'a PrecessingBinary,
    pub(crate) version: PrecessionVersion,
    /// Co-precessing modes `(l, m' > 0)` in LAL's order.
    pub(crate) coprecessing: &'a [(u32, u32)],
    pub(crate) range: AngleRange,
    /// Whether the prescription may change as in LAL.
    pub(crate) version_fallback: bool,
}

/// The angle range `(M_MIN, M_MAX)` of a set-up and its J-frame modes
/// `(l, m)`.
type Group = ((u32, u32), Vec<(u32, i32)>);

/// One J-frame mode being accumulated.
struct Accumulator {
    ell: u32,
    m: i32,
    positive: Vec<Complex<f64>>,
    negative: Vec<Complex<f64>>,
}

impl JFrameModes<'_> {
    /// `IMRPhenomXSetWaveformVariables` as called by the single-mode
    /// functions, which pass the masses and spins ordered by
    /// `XLALIMRPhenomXPCheckMassesAndSpins` (for both grids).
    fn waveform_struct(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<WaveformStruct, IMRPhenomXError> {
        let binary = self.binary;
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

    /// `IMRPhenomXSetWaveformVariables` and
    /// `IMRPhenomXGetAndSetPrecessionVariables` as called by the
    /// single-mode functions, for the angle range `(m_min, m_max)`.
    fn setup(
        &self,
        version: PrecessionVersion,
        m_min: u32,
        m_max: u32,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<Setup, IMRPhenomXError> {
        let (m1, m2, chi1, chi2) = self.binary.ordered();
        let mut wf = self.waveform_struct(delta_f, f_min, f_max)?;
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
        let modes = (prec.precessing_tag != 3).then(|| build_modes(&wf, self.coprecessing));
        Ok(Setup { wf, prec, modes })
    }

    /// The J-frame modes `(l, m)`, `l` among the co-precessing modes and
    /// `|m| <= l`, grouped by the angle range of their set-up.
    fn groups(&self) -> Vec<Group> {
        let mut ells: Vec<u32> = self.coprecessing.iter().map(|&(ell, _)| ell).collect();
        ells.sort_unstable();
        ells.dedup();
        let jframe = ells
            .iter()
            .flat_map(|&ell| (-(ell as i32)..=ell as i32).map(move |m| (ell, m)));
        let fixed = match self.range {
            AngleRange::ModeArray => {
                let m_min = self.coprecessing.iter().map(|&(_, m)| m).min();
                let m_max = self.coprecessing.iter().map(|&(_, m)| m).max();
                Some((m_min.unwrap_or(4), m_max.unwrap_or(1)))
            }
            // only the SpinTaylor angles depend on the range
            AngleRange::PerMode if self.version != PrecessionVersion::SpinTaylor => Some((1, 4)),
            AngleRange::PerMode => None,
        };
        let mut groups: Vec<Group> = Vec::new();
        for (ell, m) in jframe {
            let key = fixed.unwrap_or((m.unsigned_abs().max(1), ell));
            match groups.iter_mut().find(|(k, _)| *k == key) {
                Some((_, modes)) => modes.push((ell, m)),
                None => groups.push((key, vec![(ell, m)])),
            }
        }
        groups
    }

    /// The modes on `frequencies` (positive, increasing): the single-mode
    /// functions for a frequency sequence.
    pub(crate) fn sequence(&self, frequencies: &[f64]) -> Result<FrequencyModes, IMRPhenomXError> {
        let Some((f_min, f_max)) = frequency_bounds(frequencies)? else {
            return Ok(FrequencyModes::default());
        };
        self.evaluate(0.0, f_min, f_max, frequencies, 0, frequencies.len())
    }

    /// The modes on the uniform grid `k delta_f`, `k = 0 .. n - 1`, as
    /// `SimInspiralChooseFDModes`: zero below `f_min` and above
    /// `min(f_max, f_cut)`, `n - 1` the next power of two of
    /// `max(f_max, f_cut) / delta_f` (`f_max = 0` selects `f_cut`).
    pub(crate) fn uniform(
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
        let f_max_prime = self.waveform_struct(delta_f, f_min, f_max)?.f_max_prime;
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max_prime / delta_f) as usize + 1;
        let last_freq = if f_max_prime < f_max {
            f_max
        } else {
            f_max_prime
        };
        let n = (next_pow2((last_freq / delta_f) as usize) + 1).max(i_stop);
        let frequencies: Vec<f64> = (0..n).map(|i| i as f64 * delta_f).collect();
        let mut modes = self.evaluate(delta_f, f_min, f_max, &frequencies, i_start, i_stop)?;
        let n_full = next_pow2((last_freq / delta_f) as usize) + 1;
        modes.frequencies.truncate(n_full);
        for mode in &mut modes.modes {
            mode.positive.truncate(n_full);
            mode.negative.truncate(n_full);
        }
        Ok(modes)
    }

    /// The modes on `frequencies`, evaluated at the indices `start..stop`
    /// (zero elsewhere), for the grid bounds `f_min`, `f_max` of LAL.
    fn evaluate(
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
        let mut out: Vec<Accumulator> = Vec::new();
        let mut frame = None;
        for ((m_min, m_max), jframe) in self.groups() {
            let mut accumulators: Vec<Accumulator> = jframe
                .iter()
                .map(|&(ell, m)| Accumulator {
                    ell,
                    m,
                    positive: vec![zero; n],
                    negative: vec![zero; n],
                })
                .collect();
            let band = &frequencies[start..stop];
            let mut setup = self.setup(self.version, m_min, m_max, delta_f, f_min, f_max)?;
            if !self.version_fallback {
                check_version(&setup.prec, self.version)?;
            }
            let mut result = twist_up_modes(
                &mut setup,
                self.coprecessing,
                band,
                start,
                &mut accumulators,
            );
            if !self.version_fallback {
                result = forbid_spintaylor_fallback(result);
            }
            let setup = match result {
                Err(IMRPhenomXError::SpinTaylorFallbackToMsa(_)) => {
                    // LAL switches this mode to the MSA angles (223)
                    let mut setup =
                        self.setup(PrecessionVersion::Msa, m_min, m_max, delta_f, f_min, f_max)?;
                    for accumulator in &mut accumulators {
                        accumulator.positive.fill(zero);
                        accumulator.negative.fill(zero);
                    }
                    twist_up_modes(
                        &mut setup,
                        self.coprecessing,
                        band,
                        start,
                        &mut accumulators,
                    )?;
                    setup
                }
                Err(error) => return Err(error),
                Ok(()) => setup,
            };
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

/// `IMRPhenomXPHM_OneMode` for the J-frame modes of `accumulators` (all
/// sharing `setup`): twist up the co-precessing modes of the same `l` at
/// `band`, the frequencies of indices `offset..` of the accumulators.
fn twist_up_modes(
    setup: &mut Setup,
    coprecessing: &[(u32, u32)],
    band: &[f64],
    offset: usize,
    accumulators: &mut [Accumulator],
) -> Result<(), IMRPhenomXError> {
    if setup.wf.f_max_prime <= setup.wf.f_min {
        return Err(IMRPhenomXError::InvalidFrequencies("fCut <= f_min"));
    }
    let angles = if setup.prec.precessing_tag == 3 {
        Some(setup.initialize_euler_angles(false)?)
    } else {
        None
    };
    let modes = match setup.modes.take() {
        Some(modes) => modes,
        None => build_modes(&setup.wf, coprecessing),
    };
    let wf = &setup.wf;
    let prec = &setup.prec;
    let mf_cut = wf.f_max_prime * wf.m_sec;
    // the co-precessing modes of the l of the J-frame modes
    let modes: Vec<_> = modes
        .iter()
        .filter(|mode| accumulators.iter().any(|a| a.ell == mode.ell))
        .collect();
    let mut mprimes: Vec<u32> = modes.iter().map(|mode| mode.emm).collect();
    mprimes.sort_unstable();
    mprimes.dedup();
    let mut cursors = [AngleCursor::default(); 5];
    // per m': e^{-i k alpha} for k = -4..=4, the powers of cos(beta/2),
    // sin(beta/2) and e^{i m' epsilon}
    let mut twists = [(
        [Complex::new(0.0, 0.0); 9],
        [0.0; 9],
        [0.0; 9],
        Complex::new(0.0, 0.0),
    ); 5];
    for (i, &f) in band.iter().enumerate() {
        let mf = wf.m_sec * f;
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
                    let v = (PI * mf * (2.0 / mprime as f64)).cbrt();
                    let (alpha, epsilon, cos_beta) = prec.msa_angles(v);
                    (alpha, cos_beta, epsilon)
                }
                None => {
                    let (alpha, epsilon, cos_beta) =
                        prec.nnlo_angles(PI * mf * 2. / mprime as f64, wf.eta);
                    (alpha, cos_beta, epsilon)
                }
            };
            let (cb, sb) = half_angles(cos_beta);
            let mut c = [1.0; 9];
            let mut s = [1.0; 9];
            for k in 1..9 {
                c[k] = cb * c[k - 1];
                s[k] = sb * s[k - 1];
            }
            let rotations: [Complex<f64>; 9] =
                std::array::from_fn(|index| cexp_i(-((index as f64) - 4.0) * alpha));
            twists[mprime as usize] = (rotations, c, s, cexp_i(mprime as f64 * epsilon));
        }
        let powers = UsefulPowers::new(mf);
        for mode in &modes {
            let (rotations, c, s, e) = &twists[mode.emm as usize];
            let h = mode.generator.h(mf, &powers);
            let d = wigner_d(mode.ell, mode.emm, c, s);
            let l = mode.ell as usize;
            let sigma = if mode.ell % 2 == 0 { 1.0 } else { -1.0 };
            let factors = EpsilonFactors::new(h, *e, sigma);
            for accumulator in accumulators.iter_mut().filter(|a| a.ell == mode.ell) {
                let m = accumulator.m;
                let (positive, negative) =
                    factors.contribution(&d, l, mode.emm, m, rotations[(m + 4) as usize]);
                accumulator.positive[offset + i] += positive;
                accumulator.negative[offset + i] += negative;
            }
        }
    }
    Ok(())
}

impl IMRPhenomXPHM {
    /// The J-frame modes of the single-mode functions.
    fn jframe(&self) -> JFrameModes<'_> {
        JFrameModes {
            binary: &self.binary,
            version: self.version,
            coprecessing: &self.modes,
            range: if self.custom_mode_array {
                AngleRange::ModeArray
            } else {
                AngleRange::PerMode
            },
            version_fallback: self.version_fallback,
        }
    }

    /// The modes `h_lm` in the J frame at `frequencies` (Hz, positive and
    /// increasing), as `XLALSimIMRPhenomXPHMFrequencySequenceOneMode`
    /// computes each of them, for every `l` of the co-precessing modes and
    /// `|m| <= l`.
    pub(crate) fn jframe_modes(
        &self,
        frequencies: &[f64],
    ) -> Result<FrequencyModes, IMRPhenomXError> {
        self.jframe().sequence(frequencies)
    }

    /// The modes `h_lm` in the J frame on the uniform grid `k delta_f`, as
    /// `SimInspiralChooseFDModes` (`XLALSimIMRPhenomXPHMModes`).
    pub(crate) fn jframe_modes_uniform(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<FrequencyModes, IMRPhenomXError> {
        self.jframe().uniform(delta_f, f_min, f_max)
    }
}

/// The factors of `IMRPhenomXPHMTwistUpOneMode` that multiply the Wigner-d
/// matrix for a co-precessing mode `h = h_{l,-m'}`: `e^{-i m' epsilon} h`
/// and `e^{i m' epsilon} (-1)^l conj(h)`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct EpsilonFactors {
    positive: Complex<f64>,
    negative: Complex<f64>,
}

impl EpsilonFactors {
    /// The factors for `e = e^{i m' epsilon}` and `sigma = (-1)^l`.
    #[inline]
    pub(crate) fn new(h: Complex<f64>, e: Complex<f64>, sigma: f64) -> Self {
        Self {
            positive: Complex::new(1.0, 0.0) / e * h,
            negative: e * sigma * h.conj(),
        }
    }

    /// The factors with the antisymmetric part `h_as` of the co-precessing
    /// mode (IMRPhenomXPNR): `e^{-i m' epsilon} (h + h_as)` and
    /// `e^{i m' epsilon} (-1)^l conj(h - h_as)`.
    #[inline]
    pub(crate) fn with_antisymmetric(
        h: Complex<f64>,
        h_as: Complex<f64>,
        e: Complex<f64>,
        sigma: f64,
    ) -> Self {
        Self {
            positive: Complex::new(1.0, 0.0) / e * (h + h_as),
            negative: e * sigma * (h - h_as).conj(),
        }
    }

    /// The contributions to `h_lm(f)` and `h_lm(-f)` of the J-frame mode
    /// `(l, m)`, with `d = d^l_{k,m'}(beta)` for `k = -l..l` and `rotation =
    /// e^{-i m alpha}`.
    #[inline]
    pub(crate) fn contribution(
        &self,
        d: &[f64; 9],
        l: usize,
        mprime: u32,
        m: i32,
        rotation: Complex<f64>,
    ) -> (Complex<f64>, Complex<f64>) {
        let index = (m + l as i32) as usize;
        // d^l_{m,-m'} = (-1)^(m + m') d^l_{-m,m'}
        let d_minus = if (m + mprime as i32) % 2 == 0 {
            d[2 * l - index]
        } else {
            -d[2 * l - index]
        };
        (
            self.positive * (rotation * d_minus),
            self.negative * (rotation * d[index]),
        )
    }
}
