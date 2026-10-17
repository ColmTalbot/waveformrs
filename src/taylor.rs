use std::f64::consts::PI;

use num_complex::Complex;

use crate::{
    constants::{MPC_SI, SOLAR_RADIUS_IN_M, SOLAR_RADIUS_IN_S},
    pn::{PNPhasing, Params},
    waveform::Waveform,
};

#[derive(Clone, Copy, Debug)]
pub struct TaylorF2 {
    params: Params,
    pub phasing: PNPhasing,
    pub total_mass: f64,
    pub luminosity_distance: f64,
    /// SPA phase at the reference frequency, zero for reference at coalescence.
    reference_phase: f64,
}

impl TaylorF2 {
    pub fn new(
        total_mass: f64,
        mass_ratio: f64,
        chi_1: f64,
        chi_2: f64,
        luminosity_distance: f64,
    ) -> Self {
        let params = Params::new(mass_ratio, chi_1, chi_2);
        let mut new = Self {
            params,
            phasing: PNPhasing {
                v: [0.0; 16],
                vlogv: [0.0; 16],
                vlogvlogv: [0.0; 16],
            },
            total_mass,
            luminosity_distance: luminosity_distance * MPC_SI,
            reference_phase: 0.0,
        };
        new.phasing = new.phasing_coefficients();
        new
    }

    fn phasing_coefficients(&self) -> PNPhasing {
        PNPhasing::new(&self.params)
    }

    /// Define `phi_c` as the orbital phase at the reference frequency (in Hz)
    /// rather than at coalescence, matching the `f_ref` argument of
    /// `XLALSimInspiralTaylorF2`.
    pub fn with_reference_frequency(mut self, reference_frequency: f64) -> Self {
        self.reference_phase = 0.0;
        if reference_frequency > 0.0 {
            let v_ref = self.orbital_speed(reference_frequency);
            self.reference_phase = self.phase(v_ref, -PI / 8.0);
        }
        self
    }

    /// Leading-order amplitude prefactor, A(v) = prefactor * v^(-7/2).
    fn amplitude_prefactor(&self) -> f64 {
        let mass_1 = self.total_mass / (1.0 + self.params.mass_ratio);
        let mass_2 = self.total_mass - mass_1;
        let amp_0 =
            -4.0 * mass_1 * mass_2 * SOLAR_RADIUS_IN_M * SOLAR_RADIUS_IN_S * (PI / 12.0).sqrt()
                / self.luminosity_distance;
        amp_0 * (5.0 / 32.0 / self.params.eta).sqrt()
    }

    /// Snapshot of everything needed to evaluate the waveform, taken from the
    /// current (public, mutable) fields so batch evaluation stays consistent
    /// with the single-frequency methods.
    fn evaluator(&self, phi_c: f64) -> Evaluator {
        Evaluator {
            pi_m: PI * self.total_mass * SOLAR_RADIUS_IN_S,
            amplitude: self.amplitude_prefactor(),
            phase_shift: 2.0 * phi_c + PI / 4.0 + self.reference_phase,
            power: SparsePolynomial::new(&self.phasing.v),
            log: SparsePolynomial::new(&self.phasing.vlogv),
        }
    }
}

/// sum_{i=lo}^{hi} c_i v^i evaluated with Horner's scheme.
#[derive(Clone, Copy, Debug)]
struct SparsePolynomial {
    coefficients: [f64; 16],
    lo: usize,
    hi: usize,
    is_zero: bool,
}

impl SparsePolynomial {
    fn new(coefficients: &[f64; 16]) -> Self {
        let lo = coefficients.iter().position(|&c| c != 0.0);
        let hi = coefficients.iter().rposition(|&c| c != 0.0);
        match (lo, hi) {
            (Some(lo), Some(hi)) => Self {
                coefficients: *coefficients,
                lo,
                hi,
                is_zero: false,
            },
            _ => Self {
                coefficients: *coefficients,
                lo: 0,
                hi: 0,
                is_zero: true,
            },
        }
    }

    /// sum_i c_i v^(i - 5), with `inv_v5` = v^-5.
    #[inline(always)]
    fn eval_over_v5(&self, v: f64, inv_v5: f64) -> f64 {
        if self.is_zero {
            return 0.0;
        }
        let mut result = self.coefficients[self.hi];
        for ii in (self.lo..self.hi).rev() {
            result = result * v + self.coefficients[ii];
        }
        if self.lo >= 5 {
            for _ in 5..self.lo {
                result *= v;
            }
            result
        } else {
            for _ in 0..self.lo {
                result *= v;
            }
            result * inv_v5
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Evaluator {
    pi_m: f64,
    amplitude: f64,
    phase_shift: f64,
    power: SparsePolynomial,
    log: SparsePolynomial,
}

impl Evaluator {
    #[inline(always)]
    fn evaluate(&self, frequency: f64) -> Complex<f64> {
        let v = (self.pi_m * frequency).cbrt();
        let v2 = v * v;
        let v3 = v2 * v;
        let inv_v5 = 1.0 / (v3 * v2);
        let mut phasing = self.power.eval_over_v5(v, inv_v5);
        if !self.log.is_zero {
            phasing += self.log.eval_over_v5(v, inv_v5) * v.ln();
        }
        phasing -= self.phase_shift;
        // v^(-7/2) = v^-5 * v^(3/2)
        let amplitude = self.amplitude * inv_v5 * v * v.sqrt();
        let (sin, cos) = phasing.sin_cos();
        Complex::new(amplitude * cos, -amplitude * sin)
    }
}

impl Waveform for TaylorF2 {
    fn orbital_speed(&self, frequency: f64) -> f64 {
        (PI * self.total_mass * SOLAR_RADIUS_IN_S * frequency).cbrt()
    }

    fn amplitude(&self, v: f64) -> f64 {
        self.amplitude_prefactor() / (v.powi(3) * v.sqrt())
    }

    fn phase(&self, v: f64, phi_c: f64) -> f64 {
        let inv_v5 = 1.0 / v.powi(5);
        let mut phasing = SparsePolynomial::new(&self.phasing.v).eval_over_v5(v, inv_v5);
        phasing += SparsePolynomial::new(&self.phasing.vlogv).eval_over_v5(v, inv_v5) * v.ln();
        phasing -= 2.0 * phi_c + PI / 4.0 + self.reference_phase;
        phasing
    }

    fn waveform_single_frequency(&self, frequency: f64, phi_c: f64) -> Complex<f64> {
        self.evaluator(phi_c).evaluate(frequency)
    }

    fn waveform(&self, frequencies: &[f64], phi_c: f64) -> Vec<Complex<f64>> {
        let evaluator = self.evaluator(phi_c);
        frequencies.iter().map(|&f| evaluator.evaluate(f)).collect()
    }
}
