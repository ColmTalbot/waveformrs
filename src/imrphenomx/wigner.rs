//! Spin-weighted spherical harmonics and Wigner-d matrices used to twist up
//! the co-precessing modes of IMRPhenomXP/XPHM (Appendix A of
//! arXiv:2004.06503), following `IMRPhenomXPTwistUp22`,
//! `IMRPhenomXPTwistUp22_NumericalAngles` (`LALSimIMRPhenomX_precession.c`)
//! and `XLALSpinWeightedSphericalHarmonic` (`SphericalHarmonics.c` in LAL).

use std::f64::consts::PI;

use num_complex::Complex;

/// `pPrec->sqrt6` as hard-coded in LAL (`2.44948974278317788`).
const SQRT6: f64 = 2.449_489_742_783_178;

/// The `s = -2` spin-weighted spherical harmonics with `l = 2, 3, 4` at
/// polar angle `theta` and azimuth zero (`pPrec->Y2m2 ... Y44`), which are
/// real. Indexed by `m + l`.
#[derive(Clone, Copy, Debug, Default)]
#[allow(dead_code)] // l = 3, 4 are used by IMRPhenomXPHM
pub(crate) struct SphericalHarmonics {
    pub(crate) y2: [f64; 5],
    pub(crate) y3: [f64; 7],
    pub(crate) y4: [f64; 9],
}

impl SphericalHarmonics {
    /// `XLALSpinWeightedSphericalHarmonic(theta, 0, -2, l, m)`, with the
    /// expressions written as in LAL.
    pub(crate) fn new(theta: f64) -> Self {
        let c = theta.cos();
        let s = theta.sin();
        let ch = (theta / 2.0).cos();
        let sh = (theta / 2.0).sin();
        let s2t = (2.0 * theta).sin();
        let s3t = (3.0 * theta).sin();
        let s4t = (4.0 * theta).sin();
        let c2t = (2.0 * theta).cos();
        let y2 = [
            (5.0 / (64.0 * PI)).sqrt() * (1.0 - c) * (1.0 - c),
            (5.0 / (16.0 * PI)).sqrt() * s * (1.0 - c),
            (15.0 / (32.0 * PI)).sqrt() * s * s,
            (5.0 / (16.0 * PI)).sqrt() * s * (1.0 + c),
            (5.0 / (64.0 * PI)).sqrt() * (1.0 + c) * (1.0 + c),
        ];
        let y3 = [
            (21.0 / (2.0 * PI)).sqrt() * ch * sh.powf(5.0),
            (7.0 / (4.0 * PI)).sqrt() * (2.0 + 3.0 * c) * sh.powf(4.0),
            (35.0 / (2.0 * PI)).sqrt() * (s + 4.0 * s2t - 3.0 * s3t) / 32.0,
            ((105.0 / (2.0 * PI)).sqrt() * c * s.powf(2.0)) / 4.0,
            -(35.0 / (2.0 * PI)).sqrt() * (s - 4.0 * s2t - 3.0 * s3t) / 32.0,
            (7.0 / PI).sqrt() * ch.powf(4.0) * (-2.0 + 3.0 * c) / 2.0,
            -(21.0 / (2.0 * PI)).sqrt() * ch.powf(5.0) * sh,
        ];
        let y4 = [
            3.0 * (7.0 / PI).sqrt() * ch.powf(2.0) * sh.powf(6.0),
            3.0 * (7.0 / (2.0 * PI)).sqrt() * ch * (1.0 + 2.0 * c) * sh.powf(5.0),
            (3.0 * (9.0 + 14.0 * c + 7.0 * c2t) * sh.powf(4.0)) / (4.0 * PI.sqrt()),
            (3.0 * (3.0 * s + 2.0 * s2t + 7.0 * s3t - 7.0 * s4t)) / (32.0 * (2.0 * PI).sqrt()),
            (3.0 * (5.0 / (2.0 * PI)).sqrt() * (5.0 + 7.0 * c2t) * s.powf(2.0)) / 16.0,
            (3.0 * (3.0 * s - 2.0 * s2t + 7.0 * s3t + 7.0 * s4t)) / (32.0 * (2.0 * PI).sqrt()),
            (3.0 * ch.powf(4.0) * (9.0 - 14.0 * c + 7.0 * c2t)) / (4.0 * PI.sqrt()),
            -3.0 * (7.0 / (2.0 * PI)).sqrt() * ch.powf(5.0) * (-1.0 + 2.0 * c) * sh,
            3.0 * (7.0 / PI).sqrt() * ch.powf(6.0) * sh.powf(2.0),
        ];
        Self { y2, y3, y4 }
    }
}

/// `IMRPhenomXWignerdCoefficients_cosbeta`: `(cos(beta/2), sin(beta/2))`,
/// both non-negative.
#[inline]
pub(crate) fn half_angles(cos_beta: f64) -> (f64, f64) {
    (
        ((1.0 + cos_beta).abs() / 2.0).sqrt(),
        ((1.0 - cos_beta).abs() / 2.0).sqrt(),
    )
}

/// Twist up the co-precessing (2,2) mode `h_as` with the Euler angles
/// `alpha`, `cos(beta)` and `epsilon = -gamma`, giving `(h_+, h_x)` in the
/// J frame before the polarisation rotation by `zeta`
/// (`IMRPhenomXPTwistUp22`, `IMRPhenomXPTwistUp22_NumericalAngles`).
///
/// `y2` are the `l = 2` spherical harmonics at `(thetaJN, 0)`. The sums of
/// LAL,
///
/// `h_+ = e^{-2 i epsilon} h / 2 sum_m (A_m + B_m)`,
/// `h_x = e^{-2 i epsilon} h / 2 sum_m i (A_m - B_m)`,
///
/// with `A_m = e^{-i m alpha} d^2_{m,-2}(beta) Y_{2m}` and
/// `B_m = e^{i m alpha} d^2_{m,2}(beta) Y_{2m}`, are evaluated with real
/// arithmetic as `sum_m (a_m + b_m) cos(m alpha) + i (b_m - a_m) sin(m alpha)`
/// and `sum_m (a_m + b_m) sin(m alpha) + i (a_m - b_m) cos(m alpha)` where
/// `a_m = d^2_{m,-2} Y_{2m}` and `b_m = d^2_{m,2} Y_{2m}`.
#[inline]
pub(crate) fn twist_up_22(
    h_as: Complex<f64>,
    alpha: f64,
    cos_beta: f64,
    epsilon: f64,
    y2: &[f64; 5],
) -> (Complex<f64>, Complex<f64>) {
    let (c, s) = half_angles(cos_beta);
    let c2 = c * c;
    let c3 = c * c2;
    let c4 = c * c3;
    let s2 = s * s;
    let s3 = s * s2;
    let s4 = s * s3;
    // d22 = {d^2_{-2,2}, d^2_{-1,2}, d^2_{0,2}, d^2_{1,2}, d^2_{2,2}}
    let d22 = [s4, 2.0 * c * s3, SQRT6 * s2 * c2, 2.0 * c3 * s, c4];
    // d2m2 = {d^2_{-2,-2}, ..., d^2_{2,-2}}: d^2_{m,-2} = (-1)^m d^2_{-m,2}
    let d2m2 = [d22[4], -d22[3], d22[2], -d22[1], d22[0]];

    let (sin_a, cos_a) = alpha.sin_cos();
    let cos_2a = cos_a * cos_a - sin_a * sin_a;
    let sin_2a = cos_a * sin_a + sin_a * cos_a;
    // cos(m alpha), sin(m alpha) for m = -2..2
    let cos_m = [cos_2a, cos_a, 1.0, cos_a, cos_2a];
    let sin_m = [-sin_2a, -sin_a, 0.0, sin_a, sin_2a];

    let mut hp_sum = Complex::new(0.0, 0.0);
    let mut hc_sum = Complex::new(0.0, 0.0);
    for k in 0..5 {
        let a = d2m2[k] * y2[k];
        let b = d22[k] * y2[k];
        let sum = a + b;
        hp_sum += Complex::new(sum * cos_m[k], (b - a) * sin_m[k]);
        hc_sum += Complex::new(sum * sin_m[k], (a - b) * cos_m[k]);
    }

    let (sin_2e, cos_2e) = (2.0 * epsilon).sin_cos();
    let eps_phase = Complex::new(cos_2e, -sin_2e) * h_as * 0.5;
    (eps_phase * hp_sum, eps_phase * hc_sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The real-arithmetic sums agree with the complex expressions of LAL.
    #[test]
    fn twist_up_matches_complex_sum() {
        let ylm = SphericalHarmonics::new(0.83);
        let h = Complex::new(0.3, -1.2);
        for &(alpha, cos_beta, epsilon) in &[(0.4, 0.9, -1.3), (2.9, -0.2, 0.7), (-5.0, 0.99, 3.0)]
        {
            let (c, s) = half_angles(cos_beta);
            let d22 = [
                s.powi(4),
                2.0 * c * s.powi(3),
                SQRT6 * s * s * c * c,
                2.0 * c.powi(3) * s,
                c.powi(4),
            ];
            let d2m2 = [d22[4], -d22[3], d22[2], -d22[1], d22[0]];
            let i = Complex::new(0.0, 1.0);
            let mut hp = Complex::new(0.0, 0.0);
            let mut hc = Complex::new(0.0, 0.0);
            for m in -2i32..=2 {
                let k = (m + 2) as usize;
                let a = (i * (-m as f64) * alpha).exp() * d2m2[k] * ylm.y2[k];
                let b = (i * (m as f64) * alpha).exp() * d22[k] * ylm.y2[k];
                hp += a + b;
                hc += i * (a - b);
            }
            let phase = (-2.0 * i * epsilon).exp() * h / 2.0;
            let (hp_test, hc_test) = twist_up_22(h, alpha, cos_beta, epsilon, &ylm.y2);
            assert!((hp_test - phase * hp).norm() < 1e-15);
            assert!((hc_test - phase * hc).norm() < 1e-15);
        }
    }

    /// Spherical harmonics against `lal.SpinWeightedSphericalHarmonic`.
    #[test]
    fn spherical_harmonics_match_lal() {
        let ylm = SphericalHarmonics::new(1.1);
        // lal.SpinWeightedSphericalHarmonic(1.1, 0, -2, l, m).real
        let y2 = [
            0.04708121107949998,
            0.1535828111034953,
            0.3067985007403769,
            0.4085757574051893,
            0.3332019925372699,
        ];
        let y3 = [
            0.0608045442445269,
            0.187220244691151,
            0.33915927241411536,
            0.36818965784548674,
            0.1378888481383768,
            -0.25200913263391983,
            -0.43032442949239696,
        ];
        let y4 = [
            0.0663682599801361,
            0.20085892620442078,
            0.3547039630119613,
            0.36830885707528604,
            0.11697131486368266,
            -0.2509149381405521,
            -0.32853973589698415,
            0.06917368412004446,
            0.4696998220643001,
        ];
        let pairs = ylm
            .y2
            .iter()
            .zip(&y2)
            .chain(ylm.y3.iter().zip(&y3))
            .chain(ylm.y4.iter().zip(&y4));
        for (a, b) in pairs {
            assert!((a - b).abs() <= 2.0 * f64::EPSILON * b.abs(), "{a} != {b}");
        }
    }
}
