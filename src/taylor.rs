use std::f64::consts::PI;

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
        };
        new.phasing = new.phasing_coefficients();
        new
    }

    fn phasing_coefficients(&self) -> PNPhasing {
        PNPhasing::new(&self.params)
    }
}

impl Waveform for TaylorF2 {
    fn orbital_speed(&self, frequency: f64) -> f64 {
        (PI * self.total_mass * SOLAR_RADIUS_IN_S * frequency).cbrt()
    }

    fn amplitude(&self, v: f64) -> f64 {
        let mass_1 = self.total_mass / (1.0 + self.params.mass_ratio);
        let mass_2 = self.total_mass - mass_1;
        let amp_0: f64 =
            -4.0 * mass_1 * mass_2 * SOLAR_RADIUS_IN_M * SOLAR_RADIUS_IN_S * (PI / 12.0).sqrt()
                / self.luminosity_distance;
        let d_energy_d_flux: f64 = 5.0 / 32.0 / self.params.eta / v.powi(9);
        amp_0 * d_energy_d_flux.sqrt() * v
    }

    fn phase(&self, v: f64, phi_c: f64) -> f64 {
        let mut phasing = 0.0;
        let mut cumulative_power_frequency = v.powi(-5);
        let log_orbital_speed = v.ln();
        for ii in 0..16 {
            phasing += self.phasing.v[ii] * cumulative_power_frequency;
            phasing += self.phasing.vlogv[ii] * cumulative_power_frequency * log_orbital_speed;
            cumulative_power_frequency *= v;
        }
        phasing -= 2.0 * phi_c + PI / 4.0;

        return phasing;
    }
}
