use std::f64::consts::PI;

use crate::{
    constants::{MPC_SI, SOLAR_RADIUS_IN_M, SOLAR_RADIUS_IN_S},
    pn::{PNPhasing, Params},
    waveform::Waveform,
};

mod amplitude;
mod coefficients;
mod constants;
mod phase;
mod utils;

use amplitude::{imrphenomd_amplitude, AmplitudePrefactors};
use phase::{imrphenomd_phase, mrd_phase_derivative, PhaseCoefficients, inspiral_prefactors};

#[derive(Clone, Copy, Debug)]
pub struct IMRPhenomD {
    pub total_mass: f64,
    pub luminosity_distance: f64,
    coeffs: PhaseCoefficients,
    phase_prefactors: PNPhasing,
    amplitude_prefactors: AmplitudePrefactors,
    t0: f64,
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
        }
    }
}

impl Waveform for IMRPhenomD {
    fn phase(&self, v: f64, phi_c: f64) -> f64 {
        imrphenomd_phase(v, &self.coeffs, &self.phase_prefactors)
            - phi_c
            - self.t0 * (v - self.amplitude_prefactors.f_peak)
    }

    fn amplitude(&self, v: f64) -> f64 {
        let amp0 = 2.
            * (5.0 / (64.0 * PI)).sqrt()
            * self.total_mass
            * SOLAR_RADIUS_IN_M
            * self.total_mass
            * SOLAR_RADIUS_IN_S
            / self.luminosity_distance;
        amp0 * imrphenomd_amplitude(v, &self.amplitude_prefactors)
    }

    fn orbital_speed(&self, frequency: f64) -> f64 {
        frequency * self.total_mass * SOLAR_RADIUS_IN_S
    }
}
