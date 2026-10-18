use once_cell::sync::Lazy;

use crate::cmath::LibmCbrt;
use crate::grid::GridPowers;
use crate::imrphenomd::coefficients::FINAL_SPIN_COEFFS;
use crate::imrphenomd::constants::POWERS_OF_PI;
use crate::imrphenomd::qnm_data::{QNM_FDAMP, QNM_FRING, QNM_SPIN};
use crate::interpolate::CubicSpline;
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
        let third = value.libm_cbrt();
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

/// Powers of the dimensionless frequency shared between the amplitude and
/// phase evaluation.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MfPowers {
    pub mf: f64,
    /// Mf^(1/3)
    pub third: f64,
    /// (pi Mf)^(1/3)
    pub orbital_speed: f64,
}

impl MfPowers {
    #[inline(always)]
    pub fn new(mf: f64) -> Self {
        let third = mf.libm_cbrt();
        Self {
            mf,
            third,
            orbital_speed: third * POWERS_OF_PI.third,
        }
    }
}

/// The powers of `Mf` used at each frequency: computed on demand for a single
/// frequency ([`MfPowers`]) or read from a [`crate::grid::FrequencyGrid`]
/// ([`GridPowers`]).
pub(crate) trait PowersOfMf {
    fn mf(&self) -> f64;
    /// v = (pi Mf)^(1/3)
    fn orbital_speed(&self) -> f64;
    /// ln(v)
    fn log_orbital_speed(&self) -> f64;
    /// ln(Mf)
    fn ln_mf(&self) -> f64;
    /// 1 / Mf
    fn inv(&self) -> f64;
    /// Mf^(1/2)
    fn sqrt(&self) -> f64;
    /// Mf^(-7/6)
    fn m_seven_sixths(&self) -> f64;
}

impl PowersOfMf for MfPowers {
    #[inline(always)]
    fn mf(&self) -> f64 {
        self.mf
    }
    #[inline(always)]
    fn orbital_speed(&self) -> f64 {
        self.orbital_speed
    }
    #[inline(always)]
    fn log_orbital_speed(&self) -> f64 {
        self.orbital_speed.ln()
    }
    #[inline(always)]
    fn ln_mf(&self) -> f64 {
        self.mf.ln()
    }
    #[inline(always)]
    fn inv(&self) -> f64 {
        1.0 / self.mf
    }
    #[inline(always)]
    fn sqrt(&self) -> f64 {
        self.mf.sqrt()
    }
    #[inline(always)]
    fn m_seven_sixths(&self) -> f64 {
        // Mf^(-7/6) = 1 / (Mf * Mf^(1/6))
        1.0 / (self.mf * self.third.sqrt())
    }
}

impl PowersOfMf for GridPowers {
    #[inline(always)]
    fn mf(&self) -> f64 {
        self.itself
    }
    #[inline(always)]
    fn orbital_speed(&self) -> f64 {
        self.cbrt * POWERS_OF_PI.third
    }
    #[inline(always)]
    fn log_orbital_speed(&self) -> f64 {
        (self.ln + std::f64::consts::PI.ln()) * (1.0 / 3.0)
    }
    #[inline(always)]
    fn ln_mf(&self) -> f64 {
        self.ln
    }
    #[inline(always)]
    fn inv(&self) -> f64 {
        self.inv
    }
    #[inline(always)]
    fn sqrt(&self) -> f64 {
        self.cbrt * self.sixth
    }
    #[inline(always)]
    fn m_seven_sixths(&self) -> f64 {
        self.m_seven_sixths
    }
}

/// Solve the dense linear system `matrix x = rhs` by Gaussian elimination with
/// partial pivoting.
pub(crate) fn solve_linear_system<const N: usize>(
    mut matrix: [[f64; N]; N],
    mut rhs: [f64; N],
) -> [f64; N] {
    for col in 0..N {
        let pivot = (col..N)
            .max_by(|&a, &b| matrix[a][col].abs().total_cmp(&matrix[b][col].abs()))
            .unwrap();
        matrix.swap(col, pivot);
        rhs.swap(col, pivot);
        for row in col + 1..N {
            let factor = matrix[row][col] / matrix[col][col];
            for k in col..N {
                matrix[row][k] -= factor * matrix[col][k];
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    let mut solution = [0.0; N];
    for row in (0..N).rev() {
        let mut value = rhs[row];
        for (k, known) in solution.iter().enumerate().skip(row + 1) {
            value -= matrix[row][k] * known;
        }
        solution[row] = value / matrix[row][row];
    }
    solution
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

/// Fraction of the total mass radiated, `EradRational0815` in LALSimulation.
pub fn erad_rational_0815(params: &Params) -> f64 {
    let eta = params.eta;
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let m1s = params.m1_on_m * params.m1_on_m;
    let m2s = params.m2_on_m * params.m2_on_m;
    let s = (m1s * params.chi_1 + m2s * params.chi_2) / (m1s + m2s);
    (eta * (0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2
        + 3.352411249771192 * eta3)
        * (1.
            + (-0.0030302335878845507 - 2.0066110851351073 * eta + 7.7050567802399215 * eta2) * s))
        / (1. + (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s)
}

static FRING_SPLINE: Lazy<CubicSpline> = Lazy::new(|| CubicSpline::new(&QNM_SPIN, &QNM_FRING));
static FDAMP_SPLINE: Lazy<CubicSpline> = Lazy::new(|| CubicSpline::new(&QNM_SPIN, &QNM_FDAMP));

/// Ringdown frequency (in units of the total mass) interpolated from the QNM
/// tables and rescaled by the remnant mass.
pub fn fring(params: &Params, final_spin: f64) -> f64 {
    FRING_SPLINE.eval(final_spin) / (1.0 - erad_rational_0815(params))
}

/// Damping frequency (in units of the total mass) interpolated from the QNM
/// tables and rescaled by the remnant mass.
pub fn fdamp(params: &Params, final_spin: f64) -> f64 {
    FDAMP_SPLINE.eval(final_spin) / (1.0 - erad_rational_0815(params))
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
