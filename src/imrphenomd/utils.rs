use crate::imrphenomd::coefficients::{
    DAMPING_COEFFICIENTS, FINAL_SPIN_COEFFS, RINGDOWN_COEFFICIENTS,
};
use crate::pn::Params;

#[derive(Debug)]
pub struct Powers {
    pub third: f64,
    pub one: f64,
    pub four_thirds: f64,
    pub five_thirds: f64,
    pub two: f64,
    pub seven_thirds: f64,
    pub eight_thirds: f64,
    pub three: f64,
    pub minus_one_sixth: f64,
}

impl Powers {
    pub fn new(value: f64) -> Self {
        let third = value.cbrt();
        let sixth = third.sqrt();
        Self {
            minus_one_sixth: sixth.powi(-1),
            third,
            one: value,
            four_thirds: third.powi(4),
            five_thirds: third.powi(5),
            two: value.powi(2),
            seven_thirds: third.powi(7),
            eight_thirds: third.powi(8),
            three: value.powi(3),
        }
    }
}

pub fn final_spin_0815(params: &Params) -> f64 {
    let spin = params.m1_on_m.powi(2) * params.chi_1 + params.m2_on_m.powi(2) * params.chi_2;
    let eta = params.eta;
    let eta_powers = [1.0, eta, eta.powi(2), eta.powi(3), eta.powi(4)];
    let spin_powers = [1.0, spin, spin.powi(2), spin.powi(3), spin.powi(4)];
    FINAL_SPIN_COEFFS
        .iter()
        .zip(spin_powers.iter())
        .map(|(coeffs, s)| {
            coeffs
                .iter()
                .zip(eta_powers)
                .map(|(c, e)| c * e)
                .sum::<f64>()
                * s
        })
        .sum()
}

fn _evaluate_pade(x: f64, num_coeffs: &[f64], den_coeffs: &[f64], order: usize) -> f64 {
    let powers: Vec<f64> = (0..order).map(|i| x.powi(i as i32)).collect();
    let num = num_coeffs
        .iter()
        .zip(&powers)
        .map(|(c, p)| c * p)
        .sum::<f64>();
    let den = den_coeffs
        .iter()
        .zip(&powers)
        .map(|(c, p)| c * p)
        .sum::<f64>();
    num / den
}

pub fn fring(final_spin: f64) -> f64 {
    _evaluate_pade(
        final_spin,
        &RINGDOWN_COEFFICIENTS[0],
        &RINGDOWN_COEFFICIENTS[1],
        8,
    )
}

pub fn fdamp(final_spin: f64) -> f64 {
    _evaluate_pade(
        final_spin,
        &DAMPING_COEFFICIENTS[0],
        &DAMPING_COEFFICIENTS[1],
        7,
    )
}

pub(crate) fn phenomenological_function(eta: f64, xi: f64, coeffs: &[[f64; 3]; 4]) -> f64 {
    let etas = [1.0, eta, eta * eta];
    let xis = [1.0, xi, xi * xi, xi * xi * xi];
    coeffs
        .iter()
        .zip(xis.iter())
        .map(|(coeffs, x)| {
            coeffs
                .iter()
                .zip(etas.iter())
                .map(|(c, e)| c * e)
                .sum::<f64>()
                * x
        })
        .sum()
}
