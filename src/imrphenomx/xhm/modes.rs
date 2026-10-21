// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The modes of IMRPhenomXHM as `SimInspiralChooseFDModes` returns them
//! (`XLALSimIMRPhenomXHMModes`, with `PhenomXHMThresholdMband = 0`).

use num_complex::Complex;

use super::{IMRPhenomXHM, IMRPHENOMXHM_MODES};
use crate::imrphenomx::internals::IMRPhenomXError;
use crate::model::{FrequencyMode, FrequencyModes};

/// `NextPow2` of LALSimulation, `2^ceil(log2(n))`.
pub(crate) fn next_pow2(n: usize) -> usize {
    2.0f64.powf((n as f64).log2().ceil()) as usize
}

impl IMRPhenomXHM {
    /// The active modes in the order of LAL's mode loop (`l` then `m`).
    fn modes_in_lal_order(&self) -> Vec<(u32, i32)> {
        (2..=4u32)
            .flat_map(|ell| (-(ell as i32)..=ell as i32).map(move |emm| (ell, emm)))
            .filter(|&(ell, emm)| {
                emm != 0
                    && IMRPHENOMXHM_MODES.contains(&(ell, emm.unsigned_abs()))
                    && self.is_active(ell, emm)
            })
            .collect()
    }

    /// The active modes `h_lm` at `frequencies` (Hz, increasing) for the
    /// orbital reference phase `phi_ref`, each from
    /// `XLALSimIMRPhenomXHMFrequencySequenceOneMode` (see
    /// [`IMRPhenomXHM::mode`]): `m < 0` at positive frequencies, `m > 0` at
    /// negative frequencies, as `SimInspiralChooseFDModes`.
    pub(crate) fn frequency_modes(
        &self,
        frequencies: &[f64],
        phi_ref: f64,
        inclination: f64,
    ) -> FrequencyModes {
        let zero = vec![Complex::new(0.0, 0.0); frequencies.len()];
        let modes = self
            .modes_in_lal_order()
            .into_iter()
            .map(|(ell, emm)| {
                let h = self.mode(ell, emm, frequencies, phi_ref);
                let (positive, negative) = if emm < 0 {
                    (h, zero.clone())
                } else {
                    (zero.clone(), h)
                };
                FrequencyMode {
                    ell,
                    m: emm,
                    positive,
                    negative,
                }
            })
            .collect();
        FrequencyModes {
            frequencies: frequencies.to_vec(),
            theta: inclination,
            phi: std::f64::consts::FRAC_PI_2,
            polarization_angle: 0.0,
            modes,
        }
    }

    /// `SimInspiralChooseFDModes(..., delta_f, f_min, f_max, ...)` for
    /// IMRPhenomXHM (`XLALSimIMRPhenomXHMGenerateFDOneMode` for each mode):
    /// the modes on the frequencies `k delta_f`, `k = 0 .. n - 1`, zero
    /// below `f_min` and above `min(f_max, f_cut)`, with `n - 1` the next
    /// power of two of `max(f_max, f_cut) / delta_f` (`f_max = 0` selects
    /// `f_cut`). `f_ref = 0` selects `f_min`.
    pub(crate) fn frequency_modes_uniform(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
        phi_ref: f64,
        inclination: f64,
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
        let wf = self.waveform_struct(f_min, f_max, phi_ref)?;
        let f_max_prime = wf.f_max_prime(f_max);
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max_prime / delta_f) as usize + 1;
        if i_start >= i_stop {
            return Err(IMRPhenomXError::InvalidFrequencies("fCut <= f_min"));
        }
        let last_freq = if f_max_prime < f_max {
            f_max
        } else {
            f_max_prime
        };
        let n = next_pow2((last_freq / delta_f) as usize) + 1;
        // the reference frequency of the uniform drivers is f_min
        let mut model = self.clone();
        if model.f_ref == 0.0 {
            model.f_ref = f_min;
        }
        let band: Vec<f64> = (i_start..i_stop).map(|i| i as f64 * delta_f).collect();
        let mut modes = model.frequency_modes(&band, phi_ref, inclination);
        let zero = Complex::new(0.0, 0.0);
        let pad = |h: &mut Vec<Complex<f64>>| {
            let mut full = vec![zero; n.max(i_stop)];
            full[i_start..i_stop].copy_from_slice(h);
            full.truncate(n);
            *h = full;
        };
        for mode in &mut modes.modes {
            pad(&mut mode.positive);
            pad(&mut mode.negative);
        }
        modes.frequencies = (0..n).map(|i| i as f64 * delta_f).collect();
        Ok(modes)
    }
}
