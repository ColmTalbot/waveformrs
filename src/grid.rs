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
    /// f^(1/3)
    pub(crate) cbrt: Vec<f64>,
    /// f^(1/6)
    pub(crate) sixth: Vec<f64>,
    /// 1 / f
    pub(crate) inv: Vec<f64>,
    /// ln(f)
    pub(crate) ln: Vec<f64>,
    /// f^(-7/6)
    pub(crate) m_seven_sixths: Vec<f64>,
}

impl FrequencyGrid {
    /// Prepare `frequencies` (in Hz, all positive) for repeated evaluation.
    pub fn new(frequencies: &[f64]) -> Self {
        let cbrt: Vec<f64> = frequencies.iter().map(|f| f.cbrt()).collect();
        let sixth: Vec<f64> = cbrt.iter().map(|c| c.sqrt()).collect();
        let inv: Vec<f64> = frequencies.iter().map(|f| 1.0 / f).collect();
        let m_seven_sixths = inv
            .iter()
            .zip(&sixth)
            .map(|(inv, sixth)| inv / sixth)
            .collect();
        Self {
            frequencies: frequencies.to_vec(),
            cbrt,
            sixth,
            inv,
            ln: frequencies.iter().map(|f| f.ln()).collect(),
            m_seven_sixths,
        }
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
            cbrt: scale.cbrt * self.cbrt[index],
            sixth: scale.sixth * self.sixth[index],
            inv: scale.inv * self.inv[index],
            ln: scale.ln + self.ln[index],
            m_seven_sixths: scale.m_seven_sixths * self.m_seven_sixths[index],
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
