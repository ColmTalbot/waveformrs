use std::f64::consts::PI;

use num_complex::Complex;

use crate::{
    constants::{MPC_SI, SOLAR_RADIUS_IN_M, SOLAR_RADIUS_IN_S},
    pn::{PNPhasing, Params},
    waveform::Waveform,
};

mod amplitude;
mod coefficients;
mod constants;
mod phase;
mod qnm_data;
mod utils;

use amplitude::{imrphenomd_amplitude, imrphenomd_amplitude_from_powers, AmplitudePrefactors};
use phase::{
    imrphenomd_phase, imrphenomd_phase_from_powers, inspiral_prefactors, mrd_phase_derivative,
    PhaseCoefficients,
};
use utils::MfPowers;

#[derive(Clone, Copy, Debug)]
pub struct IMRPhenomD {
    pub total_mass: f64,
    pub luminosity_distance: f64,
    coeffs: PhaseCoefficients,
    phase_prefactors: PNPhasing,
    amplitude_prefactors: AmplitudePrefactors,
    t0: f64,
    /// Dimensionless frequency at which the linear time shift vanishes.
    reference_mf: f64,
    /// Phase subtracted so that the phase at the reference frequency is `-phi_c`.
    reference_phase: f64,
}

impl IMRPhenomD {
    pub fn new(
        total_mass: f64,
        mass_ratio: f64,
        chi_1: f64,
        chi_2: f64,
        luminosity_distance: f64,
    ) -> Self {
        let params = Params::new(mass_ratio, chi_1, chi_2);
        let coeffs = PhaseCoefficients::new(&params);
        let phase_prefactors = inspiral_prefactors(&params, &coeffs);
        let amplitude_prefactors = AmplitudePrefactors::new(&params);

        Self {
            total_mass,
            luminosity_distance: luminosity_distance * MPC_SI,
            coeffs,
            phase_prefactors,
            amplitude_prefactors,
            t0: mrd_phase_derivative(amplitude_prefactors.f_peak, &coeffs),
            reference_mf: amplitude_prefactors.f_peak,
            reference_phase: 0.0,
        }
    }

    /// Define the phase relative to a reference frequency (in Hz) following
    /// the LALSimulation convention: with `phi_c = 2 * phi_ref` the waveform
    /// matches `XLALSimIMRPhenomDGenerateFD` with reference frequency `f_ref`
    /// and orbital phase `phi_ref`.
    ///
    /// Without calling this the phase is not shifted and the time shift
    /// vanishes at the peak frequency.
    pub fn with_reference_frequency(mut self, reference_frequency: f64) -> Self {
        let mf_ref = self.orbital_speed(reference_frequency);
        self.reference_mf = mf_ref;
        self.reference_phase = imrphenomd_phase(mf_ref, &self.coeffs, &self.phase_prefactors);
        self
    }

    /// Overall amplitude scaling, read from the current public fields.
    fn amplitude_scale(&self) -> f64 {
        2. * (5.0 / (64.0 * PI)).sqrt()
            * self.total_mass
            * SOLAR_RADIUS_IN_M
            * self.total_mass
            * SOLAR_RADIUS_IN_S
            / self.luminosity_distance
    }

    /// Evaluate the waveform at dimensionless frequency `mf`, computing the
    /// shared powers of `mf` once for both amplitude and phase.
    #[inline(always)]
    fn evaluate(&self, mf: f64, amplitude_scale: f64, phase_shift: f64) -> Complex<f64> {
        let powers = MfPowers::new(mf);
        let amplitude = amplitude_scale
            * imrphenomd_amplitude_from_powers(&powers, &self.amplitude_prefactors);
        let phase = imrphenomd_phase_from_powers(&powers, &self.coeffs, &self.phase_prefactors)
            - phase_shift
            - self.t0 * (mf - self.reference_mf);
        let (sin, cos) = phase.sin_cos();
        Complex::new(amplitude * cos, -amplitude * sin)
    }
}

impl Waveform for IMRPhenomD {
    fn phase(&self, v: f64, phi_c: f64) -> f64 {
        imrphenomd_phase(v, &self.coeffs, &self.phase_prefactors)
            - phi_c
            - self.reference_phase
            - self.t0 * (v - self.reference_mf)
    }

    fn amplitude(&self, v: f64) -> f64 {
        self.amplitude_scale() * imrphenomd_amplitude(v, &self.amplitude_prefactors)
    }

    fn orbital_speed(&self, frequency: f64) -> f64 {
        frequency * self.total_mass * SOLAR_RADIUS_IN_S
    }

    fn waveform_single_frequency(&self, frequency: f64, phi_c: f64) -> Complex<f64> {
        self.evaluate(
            self.orbital_speed(frequency),
            self.amplitude_scale(),
            phi_c + self.reference_phase,
        )
    }

    fn waveform(&self, frequencies: &[f64], phi_c: f64) -> Vec<Complex<f64>> {
        let mass_scale = self.total_mass * SOLAR_RADIUS_IN_S;
        let amplitude_scale = self.amplitude_scale();
        let phase_shift = phi_c + self.reference_phase;
        frequencies
            .iter()
            .map(|&f| self.evaluate(f * mass_scale, amplitude_scale, phase_shift))
            .collect()
    }
}
