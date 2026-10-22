// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
use std::f64::consts::PI;

use num_complex::Complex;

use crate::{
    constants::{MPC_SI, MSUN_SI, MTSUN_SI, SOLAR_RADIUS_IN_M, SOLAR_RADIUS_IN_S},
    fastmath,
    grid::{FrequencyGrid, ScalePowers},
    model::Extrinsic,
    nrtidal,
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
use utils::{MfPowers, PowersOfMf};

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
    /// Phase subtracted so that the phase at the reference frequency is
    /// `-2 phi_c`.
    reference_phase: f64,
    /// Inclination and orbital reference phase used by
    /// [`crate::model::FrequencyDomainModel`].
    pub(crate) extrinsic: Extrinsic,
    /// Whether [`crate::model::FrequencyDomainModel`] references the phase to
    /// the first frequency of each grid (LAL's `f_ref = 0`).
    pub(crate) reference_at_first_frequency: bool,
    /// The NRTidalv2 corrections of IMRPhenomD_NRTidalv2.
    tidal: Option<Tidal>,
}

/// The NRTidalv2 corrections of `IMRPhenomD_NRTidal_Core` and
/// `IMRPhenomDGenerateFD` (with `NRTidalv2_V`), for the heavier body first.
#[derive(Clone, Copy, Debug)]
struct Tidal {
    /// Total mass in solar masses as LAL computes it from the masses in kg.
    mtot: f64,
    xa: f64,
    xb: f64,
    kappa2t: f64,
    /// Merger frequency (Hz): the window ends at 1.2 times it.
    f_merger: f64,
    /// `pn_fac (SS_3p5PN + SSS_3p5PN)`, the coefficient of `f^(2/3)` (Hz) of
    /// the 3.5PN spin-induced multipole phase.
    ho_spin: f64,
}

impl Tidal {
    fn new(total_mass: f64, mass_ratio: f64, params: &Params, lambda: [f64; 2]) -> Self {
        let m1_si = total_mass / (1.0 + mass_ratio) * MSUN_SI;
        let m2_si = total_mass * mass_ratio / (1.0 + mass_ratio) * MSUN_SI;
        let m1 = m1_si / MSUN_SI;
        let m2 = m2_si / MSUN_SI;
        let mtot = m1 + m2;
        let eta = m1 * m2 / (mtot * mtot);
        let (xa, xb) = (m1 / mtot, m2 / mtot);
        let kappa2t = nrtidal::kappa2t(m1_si, m2_si, lambda[0], lambda[1]);
        let pi_m = PI * (mtot * MTSUN_SI);
        let pn_fac = 3. * pi_m.powf(2. / 3.) / (128. * eta);
        let (ss, sss) = nrtidal::ho_spin_terms(
            xa,
            xb,
            params.chi_1,
            params.chi_2,
            params.qm_def_1,
            params.qm_def_2,
        );
        Self {
            mtot,
            xa,
            xb,
            kappa2t,
            f_merger: nrtidal::merger_frequency(mtot, kappa2t, m1 / m2),
            ho_spin: pn_fac * (ss + sss),
        }
    }

    /// The amplitude correction to add to the IMRPhenomD amplitude (before
    /// `amp0`), the phase to add and the window at `f` (Hz).
    #[inline]
    fn correction(&self, f: f64) -> (f64, f64, f64) {
        let amplitude =
            2. * (PI / 5.).sqrt() * nrtidal::tidal_amplitude(f, self.mtot, self.kappa2t);
        let phase = nrtidal::tidal_phase_v2(f, self.xa, self.xb, self.mtot, self.kappa2t)
            + self.ho_spin * f.powf(2. / 3.);
        let window = 1.0 - nrtidal::planck_taper(f, self.f_merger, 1.2 * self.f_merger);
        (amplitude, phase, window)
    }
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
            extrinsic: Extrinsic::default(),
            reference_at_first_frequency: false,
            tidal: None,
        }
    }

    /// IMRPhenomD_NRTidalv2 (Dietrich et al. 2019, arXiv:1905.06011):
    /// IMRPhenomD with the quadrupole-monopole terms of the neutron stars in
    /// its post-Newtonian phase and the NRTidalv2 phase, amplitude and
    /// merger window, as `XLALSimIMRPhenomDNRTidalFrequencySequence`. The
    /// arguments are those of [`IMRPhenomD::new`] (body 1 the heavier) and
    /// the tidal deformabilities `Lambda_1`, `Lambda_2` of the two bodies.
    /// The quadrupole-monopole parameters follow from them as in LAL
    /// (`XLALSimInspiralSetQuadMonParamsFromLambdas`).
    ///
    /// # Panics
    ///
    /// If a tidal deformability is negative.
    pub fn with_nrtidal_v2(
        total_mass: f64,
        mass_ratio: f64,
        chi_1: f64,
        chi_2: f64,
        luminosity_distance: f64,
        lambda_1: f64,
        lambda_2: f64,
    ) -> Self {
        assert!(
            lambda_1 >= 0.0 && lambda_2 >= 0.0,
            "IMRPhenomD_NRTidalv2: the tidal deformabilities must be non-negative"
        );
        let mut params = Params::new(mass_ratio, chi_1, chi_2);
        params.qm_def_1 = nrtidal::default_quadparam(lambda_1);
        params.qm_def_2 = nrtidal::default_quadparam(lambda_2);
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
            extrinsic: Extrinsic::default(),
            reference_at_first_frequency: false,
            tidal: Some(Tidal::new(
                total_mass,
                mass_ratio,
                &params,
                [lambda_1, lambda_2],
            )),
        }
    }

    /// Whether the model has the NRTidalv2 corrections
    /// ([`IMRPhenomD::with_nrtidal_v2`]).
    pub fn is_tidal(&self) -> bool {
        self.tidal.is_some()
    }

    /// [`Waveform::waveform`] written into `out` (of the length of
    /// `frequencies`).
    pub(crate) fn waveform_into(&self, frequencies: &[f64], phi_c: f64, out: &mut [Complex<f64>]) {
        let mass_scale = self.total_mass * SOLAR_RADIUS_IN_S;
        let amplitude_scale = self.amplitude_scale();
        let phase_shift = 2.0 * phi_c + self.reference_phase;
        for (h, &f) in out.iter_mut().zip(frequencies) {
            *h = self.evaluate(f, f * mass_scale, amplitude_scale, phase_shift);
        }
    }

    /// Define the phase relative to a reference frequency (in Hz) following
    /// the LALSimulation convention: with `phi_c = phi_ref` the waveform
    /// matches `XLALSimIMRPhenomDGenerateFD` with reference frequency `f_ref`
    /// and orbital phase `phi_ref`.
    ///
    /// Without calling this the phase is not shifted and the time shift
    /// vanishes at the peak frequency.
    pub fn with_reference_frequency(mut self, reference_frequency: f64) -> Self {
        let mf_ref = self.orbital_speed(reference_frequency);
        self.reference_at_first_frequency = false;
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

    /// Evaluate the waveform at the frequency `f` (Hz) and dimensionless
    /// frequency `mf`, computing the shared powers of `mf` once for both
    /// amplitude and phase.
    #[inline(always)]
    fn evaluate(&self, f: f64, mf: f64, amplitude_scale: f64, phase_shift: f64) -> Complex<f64> {
        let (amplitude, phase) =
            self.amplitude_phase(f, &MfPowers::new(mf), amplitude_scale, phase_shift);
        let (sin, cos) = phase.sin_cos();
        Complex::new(amplitude * cos, -amplitude * sin)
    }

    /// Amplitude and phase at the frequency `f` (Hz), `h = A exp(-i
    /// phase)`, for powers of `Mf` computed on demand or read from a
    /// [`FrequencyGrid`].
    #[inline(always)]
    fn amplitude_phase<P: PowersOfMf>(
        &self,
        f: f64,
        powers: &P,
        amplitude_scale: f64,
        phase_shift: f64,
    ) -> (f64, f64) {
        let mf = powers.mf();
        let amplitude = imrphenomd_amplitude_from_powers(powers, &self.amplitude_prefactors);
        let phase = imrphenomd_phase_from_powers(powers, &self.coeffs, &self.phase_prefactors)
            - phase_shift
            - self.t0 * (mf - self.reference_mf);
        match &self.tidal {
            None => (amplitude_scale * amplitude, phase),
            Some(tidal) => {
                // amp0 (A + A_T) exp(-i phi) window exp(-i phi_T)
                let (tidal_amplitude, tidal_phase, window) = tidal.correction(f);
                (
                    amplitude_scale * (amplitude + tidal_amplitude) * window,
                    phase + tidal_phase,
                )
            }
        }
    }
}

impl Waveform for IMRPhenomD {
    /// The phase at `Mf = v` for the orbital reference phase `phi_c`, which
    /// shifts the (2,2) phase by `2 phi_c`.
    fn phase(&self, v: f64, phi_c: f64) -> f64 {
        let phase = imrphenomd_phase(v, &self.coeffs, &self.phase_prefactors)
            - 2.0 * phi_c
            - self.reference_phase
            - self.t0 * (v - self.reference_mf);
        match &self.tidal {
            None => phase,
            Some(tidal) => {
                phase
                    + tidal
                        .correction(v / (self.total_mass * SOLAR_RADIUS_IN_S))
                        .1
            }
        }
    }

    fn amplitude(&self, v: f64) -> f64 {
        let amplitude = imrphenomd_amplitude(v, &self.amplitude_prefactors);
        match &self.tidal {
            None => self.amplitude_scale() * amplitude,
            Some(tidal) => {
                let (tidal_amplitude, _, window) =
                    tidal.correction(v / (self.total_mass * SOLAR_RADIUS_IN_S));
                self.amplitude_scale() * (amplitude + tidal_amplitude) * window
            }
        }
    }

    fn orbital_speed(&self, frequency: f64) -> f64 {
        frequency * self.total_mass * SOLAR_RADIUS_IN_S
    }

    fn waveform_single_frequency(&self, frequency: f64, phi_c: f64) -> Complex<f64> {
        self.evaluate(
            frequency,
            self.orbital_speed(frequency),
            self.amplitude_scale(),
            2.0 * phi_c + self.reference_phase,
        )
    }

    fn waveform(&self, frequencies: &[f64], phi_c: f64) -> Vec<Complex<f64>> {
        let mass_scale = self.total_mass * SOLAR_RADIUS_IN_S;
        let amplitude_scale = self.amplitude_scale();
        let phase_shift = 2.0 * phi_c + self.reference_phase;
        frequencies
            .iter()
            .map(|&f| self.evaluate(f, f * mass_scale, amplitude_scale, phase_shift))
            .collect()
    }

    fn waveform_on_grid(&self, grid: &FrequencyGrid, phi_c: f64, out: &mut [Complex<f64>]) {
        assert_eq!(out.len(), grid.len(), "output length must match the grid");
        let scale = ScalePowers::new(self.total_mass * SOLAR_RADIUS_IN_S);
        let amplitude_scale = self.amplitude_scale();
        let phase_shift = 2.0 * phi_c + self.reference_phase;
        let frequencies = grid.frequencies();
        for (index, h) in out.iter_mut().enumerate() {
            let (amplitude, phase) = self.amplitude_phase(
                frequencies[index],
                &grid.powers(&scale, index),
                amplitude_scale,
                phase_shift,
            );
            *h = Complex::new(amplitude, phase);
        }
        fastmath::polar_to_complex_conj(out);
    }
}
