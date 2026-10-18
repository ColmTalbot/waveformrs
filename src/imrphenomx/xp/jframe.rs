//! The J-frame modes `(2, m)` of IMRPhenomXP: its co-precessing (2,2) mode
//! twisted up mode by mode with the Euler angles of its polarizations, in
//! the conventions of the J-frame modes of IMRPhenomXPHM. LALSimulation has
//! no modes for IMRPhenomXP.

use std::f64::consts::PI;

use num_complex::Complex;

use super::{frequency_bounds, phif_ref_xp, IMRPhenomXP, Setup};
use crate::cmath::LibmCbrt;
use crate::imrphenomx::internals::{IMRPhenomXError, UsefulPowers};
use crate::imrphenomx::libm::cexp_i;
use crate::imrphenomx::precession::PrecessionVersion;
use crate::imrphenomx::wigner::half_angles;
use crate::imrphenomx::xphm::jframe::EpsilonFactors;
use crate::imrphenomx::xphm::wigner_d;
use crate::model::{FrequencyMode, FrequencyModes};

impl Setup {
    /// The J-frame modes `(2, m)`, `m = -2..=2`, on `frequencies` with the
    /// angles of [`Setup::generate`]: for each, `h_2m(f)` and `h_2m(-f)`.
    fn jframe_modes(&mut self, frequencies: &[f64]) -> Result<FrequencyModes, IMRPhenomXError> {
        let n = frequencies.len();
        let zero = Complex::new(0.0, 0.0);
        let mut positive = vec![vec![zero; n]; 5];
        let mut negative = vec![vec![zero; n]; 5];
        let spintaylor = if self.prec.precessing_tag == 3 {
            Some(self.spintaylor_angles(frequencies, false)?)
        } else {
            None
        };
        let wf = &self.wf;
        let prec = &self.prec;
        let xas = self.xas.as_ref().expect("XAS coefficients");
        let phif_ref = phif_ref_xp(xas, wf);
        let mf_cut = wf.f_max_prime * wf.m_sec;
        for (i, &f) in frequencies.iter().enumerate() {
            let mf = wf.m_sec * f;
            if mf > mf_cut {
                continue;
            }
            let (alpha, cos_beta, epsilon) = match &spintaylor {
                Some((angles, gamma_ref)) => (
                    angles.alpha[i],
                    angles.cos_beta[i],
                    -(angles.gamma[i] - gamma_ref) + prec.epsilon0,
                ),
                None if prec.msa.is_some() => {
                    let (alpha, epsilon, cos_beta) = prec.msa_angles((PI * mf).libm_cbrt());
                    (alpha, cos_beta, epsilon)
                }
                None => {
                    let (alpha, epsilon, cos_beta) = prec.nnlo_angles(PI * mf, wf.eta);
                    (alpha, cos_beta, epsilon)
                }
            };
            let h22 = xas.h22_with_powers(mf, &UsefulPowers::new(mf), phif_ref);
            let (cb, sb) = half_angles(cos_beta);
            let mut c = [1.0; 9];
            let mut s = [1.0; 9];
            for k in 1..9 {
                c[k] = cb * c[k - 1];
                s[k] = sb * s[k - 1];
            }
            let d = wigner_d(2, 2, &c, &s);
            let factors = EpsilonFactors::new(h22, cexp_i(2.0 * epsilon), 1.0);
            for m in -2..=2i32 {
                let rotation = cexp_i(-(m as f64) * alpha);
                let (pos, neg) = factors.contribution(&d, 2, 2, m, rotation);
                positive[(m + 2) as usize][i] = pos;
                negative[(m + 2) as usize][i] = neg;
            }
        }
        Ok(FrequencyModes {
            frequencies: frequencies.to_vec(),
            theta: prec.theta_jn,
            phi: 0.0,
            polarization_angle: prec.zeta_polarization,
            modes: positive
                .into_iter()
                .zip(negative)
                .zip(-2..=2)
                .map(|((positive, negative), m)| FrequencyMode {
                    ell: 2,
                    m,
                    positive,
                    negative,
                })
                .collect(),
        })
    }
}

impl IMRPhenomXP {
    /// The J-frame modes `(2, m)` at `frequencies` (Hz, positive and
    /// increasing), consistent with [`IMRPhenomXP::polarizations`]
    /// (including its fallbacks).
    pub(crate) fn jframe_modes(
        &self,
        frequencies: &[f64],
    ) -> Result<FrequencyModes, IMRPhenomXError> {
        let Some((f_min, f_max)) = frequency_bounds(frequencies)? else {
            return Ok(FrequencyModes::default());
        };
        match self.strict(self.setup(f_min, f_max)?.jframe_modes(frequencies)) {
            Err(IMRPhenomXError::SpinTaylorFallbackToMsa(_)) => {
                Setup::new(&self.binary, PrecessionVersion::Msa, 0.0, f_min, f_max)?
                    .jframe_modes(frequencies)
            }
            result => result,
        }
    }
}
