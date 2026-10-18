//! Frequency-domain spherical-harmonic modes `h_lm(f)`, as returned by
//! LALSimulation's `SimInspiralChooseFDModes`.

use std::f64::consts::FRAC_PI_2;

use num_complex::Complex;

use crate::imrphenomx::{spin_weighted_spherical_harmonic, Polarizations};

/// One mode `h_lm` at the positive and negative frequencies `+-f_k` of a
/// [`FrequencyModes`].
///
/// In LAL's convention a mode of an aligned-spin binary with `m < 0` lives
/// at positive frequencies and one with `m > 0` at negative frequencies,
/// `h_lm(-f) = (-1)^l conj(h_{l,-m}(f))`; the other half is zero. The modes
/// of precessing binaries have both halves.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FrequencyMode {
    /// The degree `l`.
    pub ell: u32,
    /// The order `m`.
    pub m: i32,
    /// `h_lm(f_k)`.
    pub positive: Vec<Complex<f64>>,
    /// `h_lm(-f_k)`.
    pub negative: Vec<Complex<f64>>,
}

/// The modes `h_lm(f)` of a frequency-domain model without the spherical
/// harmonics, in the conventions of `SimInspiralChooseFDModes`.
///
/// The polarizations follow from the modes as
///
/// `h_+(f) = 1/2 sum_lm [h_lm(f) Y_lm + conj(h_lm(-f)) conj(Y_lm)]`,
/// `h_x(f) = i/2 sum_lm [h_lm(f) Y_lm - conj(h_lm(-f)) conj(Y_lm)]`,
///
/// with the `s = -2` spin-weighted spherical harmonics
/// `Y_lm(theta, phi)` at the direction `(theta, phi)` of the line of sight
/// in the frame of the modes, followed by a rotation of the polarization
/// basis by `polarization_angle` (`zeta`):
/// `h_+ <- h_+ cos 2 zeta + h_x sin 2 zeta`,
/// `h_x <- h_x cos 2 zeta - h_+ sin 2 zeta` ([`FrequencyModes::polarizations`]).
///
/// * Aligned-spin models (IMRPhenomXHM, and IMRPhenomXAS, IMRPhenomD and
///   TaylorF2 with their (2,+-2) modes): as in LAL, the modes include the
///   orbital reference phase `phi_ref`, `theta` is the inclination,
///   `phi = pi / 2` and `zeta = 0`.
/// * Precessing models (IMRPhenomXPHM, IMRPhenomXP): the modes are in LAL's
///   inertial J frame (`PhenomXPHMModesL0Frame = 0`), `theta = theta_JN`,
///   `phi = 0` and `zeta` is the angle of LAL's polarization basis in the J
///   frame.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FrequencyModes {
    /// The (non-negative) frequencies `f_k` in Hz.
    pub frequencies: Vec<f64>,
    /// Polar angle of the line of sight in the frame of the modes.
    pub theta: f64,
    /// Azimuth of the line of sight in the frame of the modes.
    pub phi: f64,
    /// Rotation `zeta` of the polarization basis (zero except for the
    /// precessing models).
    pub polarization_angle: f64,
    /// The modes, ordered by `l` then `m` (LAL's order).
    pub modes: Vec<FrequencyMode>,
}

impl FrequencyModes {
    /// The mode `(ell, m)`, if present.
    pub fn mode(&self, ell: u32, m: i32) -> Option<&FrequencyMode> {
        self.modes
            .iter()
            .find(|mode| mode.ell == ell && mode.m == m)
    }

    /// The polarizations `(h_+, h_x)` at [`FrequencyModes::frequencies`]
    /// from the modes (see [`FrequencyModes`]), which agree with the
    /// model's polarizations to rounding where the model evaluates every
    /// mode (IMRPhenomXHM evaluates its higher modes in the polarizations
    /// also above the cutoff `Mf = 0.3` where the modes vanish).
    pub fn polarizations(&self) -> Polarizations {
        let n = self.frequencies.len();
        let zero = Complex::new(0.0, 0.0);
        let i = Complex::new(0.0, 1.0);
        let mut hplus = vec![zero; n];
        let mut hcross = vec![zero; n];
        for mode in &self.modes {
            let y = spin_weighted_spherical_harmonic(self.theta, self.phi, mode.ell as i32, mode.m);
            let y_conj = y.conj();
            for (((hp, hc), h_pos), h_neg) in hplus
                .iter_mut()
                .zip(hcross.iter_mut())
                .zip(&mode.positive)
                .zip(&mode.negative)
            {
                let a = h_pos * y;
                let b = h_neg.conj() * y_conj;
                *hp += 0.5 * (a + b);
                *hc += 0.5 * i * (a - b);
            }
        }
        let zeta = self.polarization_angle;
        if zeta != 0.0 {
            let (sin, cos) = (2.0 * zeta).sin_cos();
            for (p, c) in hplus.iter_mut().zip(hcross.iter_mut()) {
                let (pp, pc) = (*p, *c);
                *p = pp * cos + pc * sin;
                *c = pc * cos - pp * sin;
            }
        }
        (hplus, hcross)
    }
}

/// The modes `(2, -2)` (at positive frequencies) and `(2, 2)` (at negative
/// frequencies) of an aligned-spin (2,2)-mode model from `h22 =
/// h_{2,-2}(f)`, in the conventions of IMRPhenomXHM.
pub(crate) fn aligned_22(
    frequencies: &[f64],
    h22: Vec<Complex<f64>>,
    inclination: f64,
) -> FrequencyModes {
    let zero = vec![Complex::new(0.0, 0.0); h22.len()];
    let conjugate: Vec<Complex<f64>> = h22.iter().map(|h| h.conj()).collect();
    FrequencyModes {
        frequencies: frequencies.to_vec(),
        theta: inclination,
        phi: FRAC_PI_2,
        polarization_angle: 0.0,
        modes: vec![
            FrequencyMode {
                ell: 2,
                m: -2,
                positive: h22,
                negative: zero.clone(),
            },
            FrequencyMode {
                ell: 2,
                m: 2,
                positive: zero,
                negative: conjugate,
            },
        ],
    }
}
