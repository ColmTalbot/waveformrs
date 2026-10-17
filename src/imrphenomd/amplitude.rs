use std::f64::consts::PI;

use crate::imrphenomd::coefficients::{COLLOCATION_COEFFS, GAMMA_COEFFS, RHO_COEFFS};
use crate::imrphenomd::constants::{AMP_FJOIN_INS, POWERS_OF_PI};
use crate::imrphenomd::utils::{
    fdamp, final_spin_0815, fring, phenomenological_function, solve_linear_system, MfPowers,
};
use crate::pn::{chi_pn, PNAmplitude, Params};

#[derive(Clone, Copy, Debug)]
pub struct AmplitudePrefactors {
    amp0: f64,
    pn: PNAmplitude,
    f_ringdown: f64,
    f_damping: f64,
    pub f_peak: f64,
    gamma1: f64,
    gamma2: f64,
    gamma3: f64,
    deltas: [f64; 5],
}

impl AmplitudePrefactors {
    pub fn new(params: &Params) -> Self {
        let eta = params.eta;
        let xi = chi_pn(params) - 1.0;
        let final_spin = final_spin_0815(params);
        let mut pn = phenomd_inspiral_amplitude_pn(params);
        pn.v[7] = phenomenological_function(eta, xi, &RHO_COEFFS[0]) / POWERS_OF_PI.seven_thirds;
        pn.v[8] = phenomenological_function(eta, xi, &RHO_COEFFS[1]) / POWERS_OF_PI.eight_thirds;
        pn.v[9] = phenomenological_function(eta, xi, &RHO_COEFFS[2]) / POWERS_OF_PI.three;

        let mut new = Self {
            amp0: (2.0 / 3.0 * eta).sqrt() * POWERS_OF_PI.minus_one_sixth,
            pn,
            f_ringdown: fring(params, final_spin),
            f_damping: fdamp(params, final_spin),
            f_peak: 0.0,
            gamma1: phenomenological_function(eta, xi, &GAMMA_COEFFS[0]),
            gamma2: phenomenological_function(eta, xi, &GAMMA_COEFFS[1]),
            gamma3: phenomenological_function(eta, xi, &GAMMA_COEFFS[2]),
            deltas: [0.0, 0.0, 0.0, 0.0, 0.0],
        };
        new.f_peak = f_peak(&new);
        new.deltas = deltas(params, &new);
        new
    }
}

/// PN amplitude coefficients used in the IMRPhenomD inspiral ansatz
/// (`init_amp_ins_prefactors` in LALSimulation), expressed as coefficients of
/// powers of the orbital speed v = (pi M f)^(1/3).
fn phenomd_inspiral_amplitude_pn(params: &Params) -> PNAmplitude {
    let eta = params.eta;
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let chi1 = params.chi_1;
    let chi2 = params.chi_2;
    let chi12 = chi1 * chi1;
    let chi22 = chi2 * chi2;
    let seta = params.seta;
    let seta_plus_1 = 1.0 + seta;
    let pi = PI;
    let pi2 = PI * PI;
    let mut v = [0.0; 10];
    v[0] = 1.0;
    v[2] = (-969.0 + 1804.0 * eta) / 672.0;
    v[3] = (chi1 * (81.0 * seta_plus_1 - 44.0 * eta) + chi2 * (81.0 - 81.0 * seta - 44.0 * eta))
        / 48.0;
    v[4] = (-27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * seta_plus_1
        + 10287648.0 * chi22 * seta
        + 24.0 * (-1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 + 857304.0 * chi22) * eta
        + 35371056.0 * eta2)
        / 8.128512e6;
    v[5] = (chi2 * (-285197.0 * (-1.0 + seta) + 4.0 * (-91902.0 + 1579.0 * seta) * eta
        - 35632.0 * eta2)
        + chi1 * (285197.0 * seta_plus_1 - 4.0 * (91902.0 + 1579.0 * seta) * eta - 35632.0 * eta2)
        + 42840.0 * (-1.0 + 4.0 * eta) * pi)
        / 32256.0;
    v[6] = -(-336.0
        * (-3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2
            + 2943675504.0 * chi22)
        * eta2
        - 324322727232.0 * eta3
        - 7.0
            * (-177520268561.0 + 107414046432.0 * chi22 + 107414046432.0 * chi12 * seta_plus_1
                - 107414046432.0 * chi22 * seta
                + 11087290368.0 * (chi1 + chi2 + chi1 * seta - chi2 * seta) * pi)
        + 12.0
            * eta
            * (-545384828789.0 - 176491177632.0 * chi1 * chi2 + 202603761360.0 * chi22
                + 77616.0 * chi12 * (2610335.0 + 995766.0 * seta)
                - 77287373856.0 * chi22 * seta
                + 5841690624.0 * (chi1 + chi2) * pi
                + 21384760320.0 * pi2))
        / 6.0085960704e10;
    PNAmplitude { v }
}

pub fn inspiral_amplitude(frequency: f64, prefactors: &AmplitudePrefactors) -> f64 {
    let orbital_speed = (frequency * POWERS_OF_PI.one).cbrt();
    inspiral_amplitude_from_speed(orbital_speed, prefactors)
}

/// Horner evaluation of sum_i pn_i v^i.
#[inline(always)]
fn inspiral_amplitude_from_speed(orbital_speed: f64, prefactors: &AmplitudePrefactors) -> f64 {
    let pn = &prefactors.pn.v;
    let mut amplitude = pn[9];
    for i in (0..9).rev() {
        amplitude = amplitude * orbital_speed + pn[i];
    }
    amplitude
}

pub fn inspiral_amplitude_derivative(
    frequency: f64,
    prefactors: &AmplitudePrefactors,
) -> f64 {
    let orbital_speed = (frequency * POWERS_OF_PI.one).cbrt();
    let mut cumulative_power_frequency = 1.0;
    let mut amplitude = 0.0;
    let mut power = 0.0;
    for i in 0..10 {
        amplitude += power * cumulative_power_frequency * prefactors.pn.v[i];
        power += 1.0 / 3.0;
        cumulative_power_frequency *= orbital_speed;
    }
    amplitude / frequency
}

fn mrd_amplitude(frequency: f64, prefactors: &AmplitudePrefactors) -> f64 {
    let scaled_damping = prefactors.f_damping * prefactors.gamma3;
    let delta_f = frequency - prefactors.f_ringdown;
    prefactors.gamma1 * scaled_damping / (delta_f.powi(2) + scaled_damping.powi(2))
        * (-prefactors.gamma2 * delta_f / scaled_damping).exp()
}

fn mrd_amplitude_derivative(frequency: f64, prefactors: &AmplitudePrefactors) -> f64 {
    let scaled_damping = prefactors.f_damping * prefactors.gamma3;
    let delta_f = frequency - prefactors.f_ringdown;
    let expfactor = (delta_f * prefactors.gamma2 / scaled_damping).exp();

    ((-2.0 * delta_f * scaled_damping * prefactors.gamma1)
        / (delta_f.powi(2) + scaled_damping.powi(2))
        - (prefactors.gamma2 * prefactors.gamma1))
        / (expfactor * (delta_f.powi(2) + scaled_damping.powi(2)))
}

fn f_peak(prefactors: &AmplitudePrefactors) -> f64 {
    // NOTE: There's a problem with this expression from the paper becoming imaginary if gamma2>=1
    // Fix: if gamma2 >= 1 then set the square root term to zero.
    if !(prefactors.gamma2 > 1.0) {
        (prefactors.f_ringdown
            + (prefactors.f_damping
                * ((1.0 - prefactors.gamma2.powi(2)).sqrt() - 1.0)
                * prefactors.gamma3)
                / prefactors.gamma2)
            .abs()
    } else {
        (prefactors.f_ringdown + (-prefactors.f_damping * prefactors.gamma3) / prefactors.gamma2)
            .abs()
    }
}

fn intermediate_amplitude_collocation(params: &Params) -> f64 {
    let eta = params.eta;
    let xi = chi_pn(params) - 1.0;
    phenomenological_function(eta, xi, &COLLOCATION_COEFFS)
}

fn intermediate_amplitude(frequency: f64, prefactors: &AmplitudePrefactors) -> f64 {
    let d = &prefactors.deltas;
    d[0] + frequency * (d[1] + frequency * (d[2] + frequency * (d[3] + frequency * d[4])))
}

fn deltas(params: &Params, prefactors: &AmplitudePrefactors) -> [f64; 5] {
    let f1 = AMP_FJOIN_INS;
    let f3 = prefactors.f_peak;
    let f2 = (f1 + f3) / 2.0;

    let v1 = inspiral_amplitude(f1, prefactors);
    let d1 = inspiral_amplitude_derivative(f1, prefactors);
    let v2 = intermediate_amplitude_collocation(params);
    let v3 = mrd_amplitude(f3, prefactors);
    let d3 = mrd_amplitude_derivative(f3, prefactors);
    let f1s = [1.0, f1, f1 * f1, f1 * f1 * f1, f1 * f1 * f1 * f1];
    let f3s = [1.0, f3, f3 * f3, f3 * f3 * f3, f3 * f3 * f3 * f3];

    let f2s = [1.0, f2, f2 * f2, f2 * f2 * f2, f2 * f2 * f2 * f2];
    let matrix = [
        f1s,
        f2s,
        f3s,
        [0.0, 1.0, 2.0 * f1s[1], 3.0 * f1s[2], 4.0 * f1s[3]],
        [0.0, 1.0, 2.0 * f3s[1], 3.0 * f3s[2], 4.0 * f3s[3]],
    ];
    solve_linear_system(matrix, [v1, v2, v3, d1, d3])
}

pub fn imrphenomd_amplitude(frequency: f64, prefactors: &AmplitudePrefactors) -> f64 {
    imrphenomd_amplitude_from_powers(&MfPowers::new(frequency), prefactors)
}

/// The IMRPhenomD amplitude, including the leading `amp0 * Mf^(-7/6)` scaling.
#[inline(always)]
pub(crate) fn imrphenomd_amplitude_from_powers(
    powers: &MfPowers,
    prefactors: &AmplitudePrefactors,
) -> f64 {
    let f1 = AMP_FJOIN_INS;
    let f3 = prefactors.f_peak;
    let frequency = powers.mf;

    // Mf^(-7/6) = 1 / (Mf * Mf^(1/6))
    let prefactor = prefactors.amp0 / (frequency * powers.third.sqrt());

    let amplitude = if frequency < f1 {
        inspiral_amplitude_from_speed(powers.orbital_speed, prefactors)
    } else if frequency < f3 {
        intermediate_amplitude(frequency, prefactors)
    } else if f3 <= frequency {
        mrd_amplitude(frequency, prefactors)
    } else {
        0.0
    };
    amplitude * prefactor
}
