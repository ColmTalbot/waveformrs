use std::f64::consts::PI;

use crate::imrphenomd::coefficients::{ALPHA_COEFFS, BETA_COEFFS, SIGMA_COEFFS};
use crate::imrphenomd::constants::POWERS_OF_PI;
use crate::imrphenomd::constants::PHI_FJOIN_INS;
use crate::imrphenomd::utils::{fdamp, final_spin_0815, fring, phenomenological_function};
use crate::pn::{chi_pn, taylor_3pn_ss, PNPhasing, Params};

#[derive(Clone, Copy, Debug)]
pub(crate) struct PhaseCoefficients {
    eta_inv: f64,

    pub(crate) sigma: [f64; 4],
    beta: [f64; 3],
    alpha: [f64; 5],
    f_ringdown: f64,
    f_damping: f64,

    intermediate_connection: [f64; 2],
    mrd_connection: [f64; 2],
}

impl PhaseCoefficients {
    pub fn new(params: &Params) -> Self {
        let eta = params.eta;
        let eta_inv = 1.0 / eta;
        let chi = chi_pn(params);
        let xi = chi - 1.0;

        let sigma = [
            phenomenological_function(eta, xi, &SIGMA_COEFFS[0]),
            phenomenological_function(eta, xi, &SIGMA_COEFFS[1]),
            phenomenological_function(eta, xi, &SIGMA_COEFFS[2]),
            phenomenological_function(eta, xi, &SIGMA_COEFFS[3]),
        ];
        let beta = [
            phenomenological_function(eta, xi, &BETA_COEFFS[0]),
            phenomenological_function(eta, xi, &BETA_COEFFS[1]),
            phenomenological_function(eta, xi, &BETA_COEFFS[2]),
        ];
        let alpha = [
            phenomenological_function(eta, xi, &ALPHA_COEFFS[0]),
            phenomenological_function(eta, xi, &ALPHA_COEFFS[1]),
            phenomenological_function(eta, xi, &ALPHA_COEFFS[2]),
            phenomenological_function(eta, xi, &ALPHA_COEFFS[3]),
            phenomenological_function(eta, xi, &ALPHA_COEFFS[4]),
        ];
        let final_spin = final_spin_0815(params);
        let f_ringdown = fring(final_spin);
        let f_damping = fdamp(final_spin);

        let mut new = Self {
            eta_inv,
            sigma,
            beta,
            alpha,
            f_ringdown,
            f_damping,
            intermediate_connection: [0.0, 0.0],
            mrd_connection: [0.0, 0.0],
        };

        let prefactors = inspiral_prefactors(&params, &new);
        let (intermediate_connection, mrd_connection) =
            phase_connection_coefficients(&new, &prefactors);
        new.intermediate_connection = intermediate_connection;
        new.mrd_connection = mrd_connection;
        new
    }
}

pub(crate) fn inspiral_prefactors(params: &Params, coeffs: &PhaseCoefficients) -> PNPhasing {
    let mut pn = PNPhasing::new(&params);
    let _3pnss = taylor_3pn_ss(&params) * pn.v[0];
    pn.v[5] -= PI / 4.0;
    pn.v[6] -= _3pnss;
    pn.v[8] = coeffs.sigma[0] / params.eta / POWERS_OF_PI.one;
    pn.v[9] = coeffs.sigma[1] * 3.0 / 4.0 / params.eta / POWERS_OF_PI.four_thirds;
    pn.v[10] = coeffs.sigma[2] * 3.0 / 5.0 / params.eta / POWERS_OF_PI.five_thirds;
    pn.v[11] = coeffs.sigma[3] / 2.0 / params.eta / POWERS_OF_PI.two;
    pn
}

fn phase_connection_coefficients(
    coeffs: &PhaseCoefficients,
    prefactors: &PNPhasing,
) -> ([f64; 2], [f64; 2]) {
    let f2: f64 = coeffs.f_ringdown / 2.0;

    let d_phi_ins = inspiral_phase_derivative(PHI_FJOIN_INS, prefactors);
    let d_phi_int = intermediate_phase_derivative(PHI_FJOIN_INS, coeffs);
    let c2_int = d_phi_ins - d_phi_int;
    let c1_int = inspiral_phase(PHI_FJOIN_INS, prefactors)
        - intermediate_phase(PHI_FJOIN_INS, coeffs)
        - c2_int * PHI_FJOIN_INS;

    let intermediate_phase_connection = intermediate_phase(f2, coeffs) + c1_int + c2_int * f2;
    let d_phi_int = c2_int + intermediate_phase_derivative(f2, coeffs);
    let d_phi_mrd = mrd_phase_derivative(f2, coeffs);
    let c2_mrd = d_phi_int - d_phi_mrd;
    let c1_mrd = intermediate_phase_connection - mrd_phase(f2, coeffs) - c2_mrd * f2;

    ([c1_int, c2_int], [c1_mrd, c2_mrd])
}

fn inspiral_phase(frequency: f64, prefactors: &PNPhasing) -> f64 {
    let orbital_speed = frequency.cbrt() * POWERS_OF_PI.third;
    let logv = orbital_speed.ln();

    let mut phasing = 0.0;
    let mut cumulative_orbital_speed = orbital_speed.powi(-5);

    for ii in 0..13 {
        phasing += prefactors.v[ii] * cumulative_orbital_speed;
        phasing += prefactors.vlogv[ii] * cumulative_orbital_speed * logv;
        cumulative_orbital_speed *= orbital_speed;
    }

    phasing
}

fn inspiral_phase_derivative(frequency: f64, prefactors: &PNPhasing) -> f64 {
    let orbital_speed = frequency.cbrt() * POWERS_OF_PI.third;
    let logv = orbital_speed.ln();

    let mut phasing = 0.0;
    let mut cumulative_orbital_speed = orbital_speed.powi(-5);
    let mut power = -5.0 / 3.0;

    for ii in 0..13 {
        phasing += power * prefactors.v[ii] * cumulative_orbital_speed;
        phasing += prefactors.vlogv[ii] * cumulative_orbital_speed * (power * logv + 1.0 / 3.0);
        cumulative_orbital_speed *= orbital_speed;
        power += 1.0 / 3.0;
    }

    phasing / frequency
}

fn intermediate_phase(frequency: f64, prefactors: &PhaseCoefficients) -> f64 {
    (prefactors.beta[0] * frequency - prefactors.beta[2] / (3.0 * frequency.powi(3))
        + prefactors.beta[1] * frequency.ln())
        * prefactors.eta_inv
}

fn intermediate_phase_derivative(frequency: f64, prefactors: &PhaseCoefficients) -> f64 {
    (prefactors.beta[0] + prefactors.beta[2] / frequency.powi(4) + prefactors.beta[1] / frequency)
        * prefactors.eta_inv
}

fn mrd_phase(frequency: f64, prefactors: &PhaseCoefficients) -> f64 {
    let rho_lm = 1.0;
    let tau_lm = 1.0;
    (prefactors.alpha[0] * frequency - prefactors.alpha[1] / frequency
        + 4.0 / 3.0 * prefactors.alpha[2] * frequency.powf(0.75)
        + prefactors.alpha[3]
            * rho_lm
            * ((frequency - prefactors.alpha[4] * prefactors.f_ringdown)
                / (prefactors.f_damping * tau_lm * rho_lm))
                .atan())
        * prefactors.eta_inv
}

pub(crate) fn mrd_phase_derivative(frequency: f64, prefactors: &PhaseCoefficients) -> f64 {
    let rho_lm = 1.0;
    let tau_lm = 1.0;
    (prefactors.alpha[0]
        + prefactors.alpha[1] / frequency.powi(2)
        + prefactors.alpha[2] / frequency.powf(0.25)
        + prefactors.alpha[3]
            / (prefactors.f_damping
                * tau_lm
                * (1.0
                    + (frequency - prefactors.alpha[4] * prefactors.f_ringdown).powi(2)
                        / (prefactors.f_damping * tau_lm * rho_lm).powi(2))))
        * prefactors.eta_inv
}

pub fn imrphenomd_phase(
    frequency: f64,
    coeffs: &PhaseCoefficients,
    prefactors: &PNPhasing,
) -> f64 {
    let f1 = PHI_FJOIN_INS;
    let f2 = coeffs.f_ringdown / 2.0;
    match frequency {
        f if f < f1 => inspiral_phase(f, prefactors),
        f if f1 <= f && f < f2 => {
            intermediate_phase(f, coeffs)
                + coeffs.intermediate_connection[0]
                + coeffs.intermediate_connection[1] * f
        }
        f if f2 <= f => {
            mrd_phase(f, coeffs) + coeffs.mrd_connection[0] + coeffs.mrd_connection[1] * f
        }
        _ => 0.0,
    }
}
