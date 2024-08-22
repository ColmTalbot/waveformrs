// use autodiff::*;
use ndarray::array;
use ndarray_linalg::Solve;

use crate::imrphenomd::coefficients::{COLLOCATION_COEFFS, GAMMA_COEFFS, RHO_COEFFS};
use crate::imrphenomd::constants::{AMP_FJOIN_INS, POWERS_OF_PI};
use crate::imrphenomd::utils::{fdamp, final_spin_0815, fring, phenomenological_function};
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
        let mut pn = PNAmplitude::new(&params);
        pn.v[7] = phenomenological_function(eta, xi, &RHO_COEFFS[0]) / POWERS_OF_PI.seven_thirds;
        pn.v[8] = phenomenological_function(eta, xi, &RHO_COEFFS[1]) / POWERS_OF_PI.eight_thirds;
        pn.v[9] = phenomenological_function(eta, xi, &RHO_COEFFS[2]) / POWERS_OF_PI.three;

        let mut new = Self {
            amp0: (2.0 / 3.0 * eta).sqrt() * POWERS_OF_PI.minus_one_sixth,
            pn,
            f_ringdown: fring(final_spin),
            f_damping: fdamp(final_spin),
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

pub fn inspiral_amplitude(
    frequency: f64,
    prefactors: &AmplitudePrefactors,
) -> f64 {
    let orbital_speed = (frequency * POWERS_OF_PI.one).cbrt();
    let mut cumulative_power_frequency = 1.0;
    let mut amplitude = 0.0;
    for i in 0..10 {
        amplitude += cumulative_power_frequency * prefactors.pn.v[i];
        cumulative_power_frequency *= orbital_speed;
    }
    amplitude
}

// fn inspiral_amplitude_2(frequency: FT<f64>, prefactors: &AmplitudePrefactors) -> FT<f64> {
//     1.0 + frequency.powf((2.0 / 3.0).into()) * prefactors.two_thirds
//         + frequency.powf((4.0 / 3.0).into()) * prefactors.four_thirds
//         + frequency.powf((5.0 / 3.0).into()) * prefactors.five_thirds
//         + frequency.powf((7.0 / 3.0).into()) * prefactors.seven_thirds
//         + frequency.powf((8.0 / 3.0).into()) * prefactors.eight_thirds
//         + frequency
//             * (prefactors.one
//                 + frequency * prefactors.two
//                 + frequency.powi(2.into()) * prefactors.three)
// }

// #[allow(dead_code)]
// pub fn inspiral_derivative(frequency: f64, prefactors: &AmplitudePrefactors) -> f64 {
//     let f = |x: &[FT<f64>]| inspiral_amplitude_2(x[0], prefactors);
//     let g = grad(f, &[frequency]);
//     g[0]
// }

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
    let power = [
        1.0,
        frequency,
        frequency.powi(2),
        frequency.powi(3),
        frequency.powi(4),
    ];
    power
        .iter()
        .zip(prefactors.deltas.iter())
        .map(|(p, d)| p * d)
        .sum()
}

fn deltas(params: &Params, prefactors: &AmplitudePrefactors) -> [f64; 5] {
    let f1 = AMP_FJOIN_INS;
    let f3 = prefactors.f_peak;
    let f2 = (f1 + f3) / 2.0;

    let v1 = inspiral_amplitude(f1, prefactors);
    let d1 = inspiral_amplitude_derivative(f1, &prefactors);
    let v2 = intermediate_amplitude_collocation(params);
    let v3 = mrd_amplitude(f3, prefactors);
    let d3 = mrd_amplitude_derivative(f3, prefactors);
    let f1s = [1.0, f1, f1.powi(2), f1.powi(3), f1.powi(4)];
    let f3s = [1.0, f3, f3.powi(2), f3.powi(3), f3.powi(4)];

    let matrix = array![
        f1s,
        [1.0, f2, f2.powi(2), f2.powi(3), f2.powi(4)],
        f3s,
        [0.0, 1.0, 2.0 * f1s[1], 3.0 * f1s[2], 4.0 * f1s[3]],
        [0.0, 1.0, 2.0 * f3s[1], 3.0 * f3s[2], 4.0 * f3s[3]],
    ];
    let result = array![v1, v2, v3, d1, d3];
    let solution = matrix.solve(&result).unwrap();

    [
        solution[0],
        solution[1],
        solution[2],
        solution[3],
        solution[4],
    ]
}

pub fn imrphenomd_amplitude(frequency: f64, prefactors: &AmplitudePrefactors) -> f64 {
    let f1 = AMP_FJOIN_INS;
    let f3 = prefactors.f_peak;

    let prefactor = prefactors.amp0 * frequency.powf(-7.0 / 6.0);

    let amplitude = match frequency {
        f if f < f1 => inspiral_amplitude(f, prefactors),
        f if f1 <= f && f < f3 => intermediate_amplitude(f, prefactors),
        f if f3 <= f => mrd_amplitude(f, prefactors),
        _ => 0.0,
    };
    amplitude * prefactor
}
