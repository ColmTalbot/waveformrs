// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! Frequency grids with precomputed, parameter-independent powers.
//!
//! Every model in the crate depends on the frequency only through the
//! dimensionless combination `s f` (for example `Mf`, or `pi M f`), so the
//! transcendental functions of the frequency factorise:
//! `cbrt(s f) = cbrt(s) cbrt(f)`, `ln(s f) = ln(s) + ln(f)`. When many
//! waveforms are evaluated on the same frequencies (as in a stochastic
//! sampler) the `cbrt`, `ln` and divisions of the frequency can be done once
//! per grid, leaving only multiplications per waveform.
//!
//! The results agree with the per-frequency evaluation to rounding (the
//! products are not correctly rounded, so the last bits differ).

/// Frequencies (in Hz) with their parameter-independent powers.
#[derive(Clone, Debug)]
pub struct FrequencyGrid {
    frequencies: Vec<f64>,
    /// `f^(1/3)`, `f^(1/6)`, `1 / f`, `ln(f)` and `f^(-7/6)` one after the
    /// other, each starting `stride` values after the previous one.
    powers: Vec<f64>,
    stride: usize,
}

/// Index of each power of the frequency in [`FrequencyGrid`].
const CBRT: usize = 0;
const SIXTH: usize = 1;
const INV: usize = 2;
const LN: usize = 3;
const M_SEVEN_SIXTHS: usize = 4;

impl FrequencyGrid {
    /// Prepare `frequencies` (in Hz, all positive) for repeated evaluation.
    pub fn new(frequencies: &[f64]) -> Self {
        let n = frequencies.len();
        // The powers are read together at each frequency. Stored a whole
        // number of pages apart (as separately allocated arrays often are)
        // they map to the same cache sets and the same 4 KiB offsets as each
        // other, which can make the loops over the grid twice as slow; the
        // stride staggers them by 640 bytes.
        let stride = n.div_ceil(512) * 512 + 80;
        let mut powers = vec![0.0; 5 * stride];
        for (i, &f) in frequencies.iter().enumerate() {
            let cbrt = f.cbrt();
            let sixth = cbrt.sqrt();
            let inv = 1.0 / f;
            powers[CBRT * stride + i] = cbrt;
            powers[SIXTH * stride + i] = sixth;
            powers[INV * stride + i] = inv;
            powers[LN * stride + i] = f.ln();
            powers[M_SEVEN_SIXTHS * stride + i] = inv / sixth;
        }
        Self {
            frequencies: frequencies.to_vec(),
            powers,
            stride,
        }
    }

    /// The power `which` of each frequency.
    #[inline(always)]
    fn column(&self, which: usize) -> &[f64] {
        let start = which * self.stride;
        &self.powers[start..start + self.frequencies.len()]
    }

    /// The uniform grid `f_min, f_min + delta_f, ...` up to and including
    /// `f_max` (to within half a step).
    pub fn uniform(f_min: f64, f_max: f64, delta_f: f64) -> Self {
        let n = ((f_max - f_min) / delta_f + 0.5).floor() as usize + 1;
        let frequencies: Vec<f64> = (0..n).map(|i| f_min + i as f64 * delta_f).collect();
        Self::new(&frequencies)
    }

    /// The frequencies of the grid in Hz.
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }

    /// Number of frequencies.
    pub fn len(&self) -> usize {
        self.frequencies.len()
    }

    /// Whether the grid is empty.
    pub fn is_empty(&self) -> bool {
        self.frequencies.is_empty()
    }
}

/// Powers of the scale `s` that turns a frequency into the model's
/// dimensionless frequency `s f`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ScalePowers {
    pub(crate) value: f64,
    pub(crate) cbrt: f64,
    pub(crate) sixth: f64,
    pub(crate) inv: f64,
    pub(crate) ln: f64,
    pub(crate) m_seven_sixths: f64,
}

impl ScalePowers {
    pub(crate) fn new(scale: f64) -> Self {
        let cbrt = scale.cbrt();
        let sixth = cbrt.sqrt();
        let inv = 1.0 / scale;
        Self {
            value: scale,
            cbrt,
            sixth,
            inv,
            ln: scale.ln(),
            m_seven_sixths: inv / sixth,
        }
    }
}

/// Powers of the dimensionless frequency `x = s f` at one grid point.
#[derive(Clone, Copy, Debug)]
pub(crate) struct GridPowers {
    /// x
    pub(crate) itself: f64,
    /// x^(1/3)
    pub(crate) cbrt: f64,
    /// x^(1/6)
    pub(crate) sixth: f64,
    /// 1 / x
    pub(crate) inv: f64,
    /// ln(x)
    pub(crate) ln: f64,
    /// x^(-7/6)
    pub(crate) m_seven_sixths: f64,
}

impl FrequencyGrid {
    /// Powers of `s f_i` for the scale `s`.
    #[inline(always)]
    pub(crate) fn powers(&self, scale: &ScalePowers, index: usize) -> GridPowers {
        GridPowers {
            itself: scale.value * self.frequencies[index],
            cbrt: scale.cbrt * self.column(CBRT)[index],
            sixth: scale.sixth * self.column(SIXTH)[index],
            inv: scale.inv * self.column(INV)[index],
            ln: scale.ln + self.column(LN)[index],
            m_seven_sixths: scale.m_seven_sixths * self.column(M_SEVEN_SIXTHS)[index],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn powers_match_direct_evaluation() {
        let grid = FrequencyGrid::uniform(20.0, 1024.0, 0.125);
        assert_eq!(grid.len(), 8033);
        let scale = ScalePowers::new(60.0 * 4.925490947641267e-6);
        for index in [0, 1, 1000, 8032] {
            let x = scale.value * grid.frequencies()[index];
            let p = grid.powers(&scale, index);
            let rel = |a: f64, b: f64| ((a - b) / b).abs();
            assert!(rel(p.itself, x) < 1e-15);
            assert!(rel(p.cbrt, x.cbrt()) < 4e-16);
            assert!(rel(p.sixth, x.powf(1.0 / 6.0)) < 4e-16);
            assert!(rel(p.inv, 1.0 / x) < 4e-16);
            assert!((p.ln - x.ln()).abs() < 1e-15);
            assert!(rel(p.m_seven_sixths, x.powf(-7.0 / 6.0)) < 1e-15);
        }
    }
}
