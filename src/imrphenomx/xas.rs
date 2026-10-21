// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! The aligned-spin (2,2)-mode model IMRPhenomXAS, following
//! `XLALSimIMRPhenomXASFrequencySequence` and `IMRPhenomXASGenerateFD` in
//! `LALSimIMRPhenomX.c`, and the polarisations built by
//! `XLALSimInspiralChooseFDWaveformSequence`.

use std::f64::consts::{FRAC_PI_4, PI};

use num_complex::Complex;

use crate::fastmath;
use crate::grid::{FrequencyGrid, ScalePowers};
use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::{fma, IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::nrtidal::Tidal;
use crate::imrphenomx::phase::{
    dphase_22, phase_22, phase_22_connection_coefficients, time_shift_22, PhaseCoefficients,
};
use crate::model::Extrinsic;
use crate::nrtidal::{NRTidalVersion, TidalSeries};
use crate::waveform::Waveform;

/// IMRPhenomXAS (Pratten et al. 2020, arXiv:2001.11412), the (2,2) mode of
/// non-precessing binary black holes, with the default LALSimulation model
/// versions (InsPhaseVersion 104, IntPhaseVersion 105, RDPhaseVersion 105,
/// InsAmpVersion 103, IntAmpVersion 104, RDAmpVersion 103).
///
/// # Conventions
///
/// With `phi_c = phi_ref`, the orbital reference phase passed to
/// LALSimulation, [`Waveform::waveform_modes`] reproduces the
/// `hplus`, `hcross` of `SimInspiralChooseFDWaveformSequence(phi_ref, m1, m2,
/// 0, 0, chi_1, 0, 0, chi_2, f_ref, distance, inclination, ..., IMRPhenomXAS,
/// frequencies)`. LALSimulation builds the polarisations from the (2,2) mode
/// as `hplus = (1 + cos^2 i) / 2 * Y h22` and `hcross = -i cos(i) Y h22` with
/// `Y = 2 sqrt(5 / 64 pi) exp(-i pi)`, the inclination dependence used by the
/// default `waveform_modes`; the constant `Y` is absorbed into
/// [`Waveform::waveform`], which returns `Y h22`, i.e. `hplus` for a face-on
/// binary.
///
/// The phase follows the [`Waveform`] convention `h = A exp(-i Psi)`, so
/// [`Waveform::phase`] returns `Psi = pi - phi_22` where `phi_22` is the LAL
/// phase of the (2,2) mode (`h22 = A exp(i phi_22)`). The frequency variable
/// used by [`Waveform::phase`] and [`Waveform::amplitude`] is the
/// dimensionless frequency `Mf` returned by [`Waveform::orbital_speed`].
///
/// Like LALSimulation, frequencies above the cutoff `Mf = 0.3` (0.33 for
/// `chi_eff > 0.99`) give zero.
///
/// # Reference frequency
///
/// [`IMRPhenomXAS::with_reference_frequency`] sets the frequency at which
/// the phase is referenced (LAL's `f_ref`). If it is not set (or set to
/// zero), [`Waveform::waveform`] and [`Waveform::waveform_modes`] use the
/// first frequency of the grid as LALSimulation does. The single-frequency
/// methods ([`Waveform::phase`], [`Waveform::waveform_single_frequency`])
/// cannot see the grid and, without a reference frequency, apply no
/// reference-phase alignment (the linear time shift `linb` is still applied).
#[derive(Clone, Copy, Debug)]
pub struct IMRPhenomXAS {
    /// Waveform variables (`pWF`).
    pub(crate) wf: WaveformStruct,
    /// Amplitude coefficients of the (2,2) mode (`pAmp22`).
    pub(crate) amp22: AmpCoefficients,
    /// Phase coefficients of the (2,2) mode, including the connection
    /// coefficients (`pPhase22`).
    pub(crate) phase22: PhaseCoefficients,
    /// Linear-in-frequency phase coefficient from `IMRPhenomX_TimeShift_22`.
    pub(crate) linb: f64,
    /// Constant phase shift `lina` (always zero for XAS).
    pub(crate) lina: f64,
    /// `-(phi22(MfRef) / eta + linb MfRef + lina)`: `phifRef` without the
    /// `2 phi0 + pi / 4` terms, zero if no reference frequency is set.
    reference_phase_shift: f64,
    /// Whether a reference frequency was supplied.
    has_reference_frequency: bool,
    /// 1 / eta
    inveta: f64,
    /// `amp0 * ampNorm`
    amp_prefactor: f64,
    /// Dimensionless cutoff `f_cut * m_sec` for the single-frequency methods.
    mf_cut: f64,
    /// Inclination and orbital reference phase used by
    /// [`crate::model::FrequencyDomainModel`].
    pub(crate) extrinsic: Extrinsic,
    /// The NRTidal corrections (IMRPhenomXAS_NRTidalv2/v3), if the waveform
    /// struct has them.
    pub(crate) tidal: Option<Tidal>,
}

/// `2 sqrt(5 / (64 pi)) exp(-2 i pi / 2)`, the factor from the spin-weighted
/// spherical harmonics applied by `ChooseFDWaveformSequence`.
fn ylm_factor() -> Complex<f64> {
    2.0 * (5.0 / (64.0 * PI)).sqrt() * Complex::new(0.0, -2.0 * (PI / 2.0)).exp()
}

impl IMRPhenomXAS {
    /// Construct the model from the total mass (solar masses), the mass
    /// ratio `m2 / m1 <= 1`, the aligned dimensionless spins and the
    /// luminosity distance (Mpc).
    ///
    /// # Panics
    ///
    /// If LALSimulation would refuse the parameters (see
    /// [`IMRPhenomXAS::from_components`]).
    pub fn new(
        total_mass: f64,
        mass_ratio: f64,
        chi_1: f64,
        chi_2: f64,
        luminosity_distance: f64,
    ) -> Self {
        let mass_1 = total_mass / (1.0 + mass_ratio);
        let mass_2 = total_mass - mass_1;
        Self::from_components(mass_1, mass_2, chi_1, chi_2, luminosity_distance)
            .unwrap_or_else(|e| panic!("IMRPhenomXAS: {e}"))
    }

    /// Construct the model from the component masses (solar masses, in any
    /// order: like LAL the heavier body is labelled 1), the aligned spins and
    /// the luminosity distance (Mpc).
    ///
    /// Returns an error where LALSimulation would: non-positive masses, mass
    /// ratio above 1000, spins outside [-1, 1], negative distance, or
    /// unphysical fits for the remnant.
    pub fn from_components(
        mass_1: f64,
        mass_2: f64,
        chi_1: f64,
        chi_2: f64,
        luminosity_distance: f64,
    ) -> Result<Self, IMRPhenomXError> {
        let wf = WaveformStruct::new(
            mass_1,
            mass_2,
            chi_1,
            chi_2,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            luminosity_distance,
            0.0,
        )?;
        Ok(Self::from_waveform_struct(wf, false))
    }

    /// IMRPhenomXAS_NRTidalv2 (`version` [`NRTidalVersion::V2`], Colleoni
    /// et al. 2023, arXiv:2311.15978) or IMRPhenomXAS_NRTidalv3
    /// ([`NRTidalVersion::V3`], Abac et al. 2024, arXiv:2311.07456): the
    /// arguments of [`IMRPhenomXAS::from_components`] and the tidal
    /// deformabilities `Lambda_1`, `Lambda_2` of the two bodies, as
    /// `XLALSimIMRPhenomXASFrequencySequence` with `PhenomXTidalFlag` 2 or
    /// 3. The quadrupole-monopole parameters follow from the
    /// deformabilities as in LAL.
    ///
    /// The tidal corrections depend on the whole frequency grid, as in LAL:
    /// the time shift is set at the merger frequency or at the last
    /// frequency of the grid if it ends below it, and the NRTidalv3 phase is
    /// held at its first minimum on the grid after merger. The
    /// single-frequency methods of [`Waveform`] use the time shift at the
    /// merger frequency and no minimum.
    #[allow(clippy::too_many_arguments)]
    pub fn with_nrtidal(
        version: NRTidalVersion,
        mass_1: f64,
        mass_2: f64,
        chi_1: f64,
        chi_2: f64,
        luminosity_distance: f64,
        lambda_1: f64,
        lambda_2: f64,
    ) -> Result<Self, IMRPhenomXError> {
        let mut wf = WaveformStruct::new(
            mass_1,
            mass_2,
            chi_1,
            chi_2,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            luminosity_distance,
            0.0,
        )?;
        wf.set_nrtidal(version, mass_1, mass_2, lambda_1, lambda_2)?;
        Ok(Self::from_waveform_struct(wf, false))
    }

    /// The version of the NRTidal corrections, if the model has them.
    pub fn nrtidal_version(&self) -> Option<NRTidalVersion> {
        self.tidal.map(|tidal| tidal.version)
    }

    /// Build the coefficients for given waveform variables, as
    /// `IMRPhenomXASGenerateFD` does. If `use_reference_frequency`, the phase
    /// is referenced to `wf.mf_ref`.
    pub(crate) fn from_waveform_struct(wf: WaveformStruct, use_reference_frequency: bool) -> Self {
        let amp22 = AmpCoefficients::new(&wf);
        let mut phase22 = PhaseCoefficients::new(&wf);
        phase_22_connection_coefficients(&mut phase22);
        let linb = time_shift_22(&phase22, &wf);
        let mut new = Self {
            wf,
            amp22,
            phase22,
            linb,
            lina: 0.0,
            reference_phase_shift: 0.0,
            has_reference_frequency: false,
            inveta: 1.0 / wf.eta,
            amp_prefactor: wf.amp0 * wf.amp_norm,
            mf_cut: wf.f_cut * wf.m_sec,
            extrinsic: Extrinsic::default(),
            tidal: wf.nrtidal.map(|version| Tidal::new(&wf, version)),
        };
        if use_reference_frequency {
            new.has_reference_frequency = true;
            new.reference_phase_shift = new.reference_phase_shift(wf.mf_ref);
        }
        new
    }

    /// Reference the phase to the frequency `reference_frequency` (Hz), LAL's
    /// `f_ref`: for the orbital reference phase `phi_c` the (2,2) phase at
    /// `f_ref` is `2 phi_c + pi / 4` before the spherical-harmonic factor. Zero
    /// selects the first frequency of the grid passed to
    /// [`Waveform::waveform`], as in LALSimulation.
    pub fn with_reference_frequency(mut self, reference_frequency: f64) -> Self {
        assert!(
            reference_frequency >= 0.0,
            "IMRPhenomXAS: {}",
            IMRPhenomXError::NegativeReferenceFrequency
        );
        self.wf.set_reference_frequency(reference_frequency);
        if reference_frequency > 0.0 {
            self.has_reference_frequency = true;
            self.reference_phase_shift = self.reference_phase_shift(self.wf.mf_ref);
        } else {
            self.has_reference_frequency = false;
            self.reference_phase_shift = 0.0;
        }
        self
    }

    /// Total mass in solar masses.
    pub fn total_mass(&self) -> f64 {
        self.wf.m_tot
    }

    /// `-(phi22(MfRef) / eta + phiTfRef + linb MfRef + lina)`, evaluated
    /// as in `IMRPhenomXASGenerateFD` (`phiTfRef = -phi_T(MfRef)` with the
    /// NRTidal corrections, zero otherwise).
    fn reference_phase_shift(&self, mf_ref: f64) -> f64 {
        let powers = UsefulPowers::new(mf_ref);
        let mut phi = self.inveta * phase_22(mf_ref, &powers, &self.phase22);
        if let Some(tidal) = &self.tidal {
            phi += -tidal.phase(mf_ref, &powers);
        }
        -(fma(self.linb, mf_ref, phi) + self.lina)
    }

    /// `phiTfRef`, the tidal contribution to the reference phase (zero
    /// without the NRTidal corrections).
    pub(crate) fn tidal_reference_phase(&self) -> f64 {
        match &self.tidal {
            Some(tidal) => -tidal.phase(self.wf.mf_ref, &UsefulPowers::new(self.wf.mf_ref)),
            None => 0.0,
        }
    }

    /// The model for a grid whose last frequency is `f_last` (Hz): with the
    /// NRTidal corrections, `linb` includes the time shift `tshift` of
    /// `IMRPhenomXASGenerateFD`, set by the phase derivative at the merger
    /// frequency or at `f_last` if it is lower. Without them, `self`.
    pub(crate) fn for_grid(&self, f_last: f64) -> Self {
        let Some(tidal) = &self.tidal else {
            return *self;
        };
        let mut model = *self;
        let f_final = f_last.min(tidal.f_merger());
        let mf_final = f_final * self.wf.m_sec;
        let powers = UsefulPowers::new(mf_final);
        let dphi_fmerger = 1. / self.wf.eta * dphase_22(mf_final, &powers, &self.phase22)
            + self.linb
            - tidal.phase_derivative(mf_final, &powers);
        model.linb += -dphi_fmerger;
        if model.has_reference_frequency {
            model.reference_phase_shift = model.reference_phase_shift(model.wf.mf_ref);
        }
        model
    }

    /// The (2,2) mode on `frequencies` (which set the grid-dependent tidal
    /// terms) times `factor`, written into `out`, with the NRTidal
    /// corrections, for the orbital reference phase `phi_c`.
    fn tidal_h22_into(
        &self,
        tidal: &Tidal,
        frequencies: &[f64],
        phi_c: f64,
        factor: Complex<f64>,
        out: &mut [Complex<f64>],
    ) {
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return;
        };
        let model = self.for_grid(f_max);
        let (phif_ref, mf_cut) = model.grid_phase_and_cutoff(f_min, f_max, phi_c);
        let series = tidal.series(frequencies);
        for (i, (h, &f)) in out.iter_mut().zip(frequencies).enumerate() {
            let mf = model.wf.m_sec * f;
            *h = if mf <= mf_cut {
                model.h22_tidal(tidal, mf, &UsefulPowers::new(mf), phif_ref, &series, i) * factor
            } else {
                Complex::new(0.0, 0.0)
            };
        }
    }

    /// The (2,2) mode with the NRTidal corrections at point `index` of
    /// `series`, for precomputed powers of `mf`, without the cutoff.
    #[inline]
    pub(crate) fn h22_tidal(
        &self,
        tidal: &Tidal,
        mf: f64,
        powers: &UsefulPowers,
        phif_ref: f64,
        series: &TidalSeries,
        index: usize,
    ) -> Complex<f64> {
        let phi = self.phase_22_full(mf, powers, phif_ref);
        let amp = self.amp22.ansatz(mf, powers);
        tidal.h22(&self.wf, powers, amp, phi, series, index)
    }

    /// The single-frequency (2,2) mode `h22(Mf)` with the NRTidal
    /// corrections for the time shift at the merger frequency, zero above
    /// the cutoff.
    fn tidal_h22_single(&self, tidal: &Tidal, mf: f64, two_phi0: f64) -> Complex<f64> {
        if mf > self.mf_cut {
            return Complex::new(0.0, 0.0);
        }
        let model = self.for_grid(f64::INFINITY);
        let series = tidal.series(&[mf / self.wf.m_sec]);
        let powers = UsefulPowers::fast(mf);
        model.h22_tidal(tidal, mf, &powers, model.phif_ref(two_phi0), &series, 0)
    }

    /// `phifRef` for the gravitational-wave phase `2 phi0` (twice LAL's
    /// orbital phase `phi0`): the phase offset added to the (2,2) phase in
    /// `IMRPhenomXASGenerateFD`.
    #[inline]
    pub(crate) fn phif_ref(&self, two_phi0: f64) -> f64 {
        self.reference_phase_shift + two_phi0 + FRAC_PI_4
    }

    /// The (2,2) phase `phi_22` at dimensionless frequency `mf` for a
    /// given `phifRef` (LAL convention, `h22 = A exp(i phi_22)`).
    #[inline]
    pub(crate) fn phase_22_full(&self, mf: f64, powers: &UsefulPowers, phif_ref: f64) -> f64 {
        let phi = self.phase22.phase_xas(mf, powers) * self.inveta;
        phi + (fma(self.linb, mf, self.lina) + phif_ref)
    }

    /// The (2,2) amplitude `amp0 ampNorm f^(-7/6) A(f)` (strain units).
    #[inline]
    pub(crate) fn amplitude_22_full(&self, mf: f64, powers: &UsefulPowers) -> f64 {
        self.amp_prefactor * powers.m_seven_sixths * self.amp22.ansatz(mf, powers)
    }

    /// The (2,2) mode `h22(Mf)` exactly as in the main loop of
    /// `IMRPhenomXASGenerateFD`, zero above `mf_cut`.
    #[inline]
    pub(crate) fn h22(&self, mf: f64, phif_ref: f64, mf_cut: f64) -> Complex<f64> {
        if mf <= mf_cut {
            self.h22_with_powers(mf, &UsefulPowers::fast(mf), phif_ref)
        } else {
            Complex::new(0.0, 0.0)
        }
    }

    /// The (2,2) mode at `mf` for precomputed powers of `mf`, without the
    /// cutoff.
    #[inline]
    pub(crate) fn h22_with_powers(
        &self,
        mf: f64,
        powers: &UsefulPowers,
        phif_ref: f64,
    ) -> Complex<f64> {
        let phi = self.phase_22_full(mf, powers, phif_ref);
        let amp = self.amplitude_22_full(mf, powers);
        let (sin, cos) = phi.sin_cos();
        Complex::new(amp * cos, amp * sin)
    }
}

/// The (2,2) mode of `IMRPhenomXASGenerateFD` with the NRTidal corrections
/// on one grid (as IMRPhenomXPHM calls it): the model with the time shift
/// of the grid, its `phifRef` and the tidal series on the grid.
pub(crate) struct TidalGrid22 {
    xas: IMRPhenomXAS,
    tidal: Tidal,
    phif_ref: f64,
    series: TidalSeries,
}

impl TidalGrid22 {
    /// For a model with the NRTidal corrections, the gravitational-wave
    /// phase `2 phi0` and the frequencies (Hz) of the grid.
    pub(crate) fn new(xas: &IMRPhenomXAS, two_phi0: f64, frequencies: &[f64]) -> Self {
        let tidal = xas.tidal.expect("NRTidal corrections");
        let xas = xas.for_grid(frequencies.last().copied().unwrap_or(f64::INFINITY));
        Self {
            phif_ref: xas.phif_ref(two_phi0),
            series: tidal.series(frequencies),
            xas,
            tidal,
        }
    }

    /// `h_{2,-2}` at the point `index` of the grid, at `mf`.
    #[inline]
    pub(crate) fn h22(&self, index: usize, mf: f64, powers: &UsefulPowers) -> Complex<f64> {
        self.xas
            .h22_tidal(&self.tidal, mf, powers, self.phif_ref, &self.series, index)
    }

    /// Amplitude and phase of [`TidalGrid22::h22`], `h = A exp(i phi)`.
    #[inline]
    pub(crate) fn amplitude_phase(
        &self,
        index: usize,
        mf: f64,
        powers: &UsefulPowers,
    ) -> (f64, f64) {
        let phi = self.xas.phase_22_full(mf, powers, self.phif_ref);
        let amp = self.xas.amp22.ansatz(mf, powers);
        self.tidal
            .amplitude_phase(&self.xas.wf, powers, amp, phi, &self.series, index)
    }
}

impl IMRPhenomXAS {
    /// `phifRef` and the cutoff `Mf` for a grid from `f_min` to `f_max`, as
    /// in [`Waveform::waveform`], for the orbital reference phase `phi_c`.
    fn grid_phase_and_cutoff(&self, f_min: f64, f_max: f64, phi_c: f64) -> (f64, f64) {
        let two_phi0 = 2.0 * phi_c;
        let phif_ref = if self.has_reference_frequency {
            self.phif_ref(two_phi0)
        } else {
            let mut wf = self.wf;
            wf.f_min = f_min;
            wf.set_reference_frequency(0.0);
            self.reference_phase_shift(wf.mf_ref) + two_phi0 + FRAC_PI_4
        };
        (phif_ref, self.wf.f_max_prime(f_max) * self.wf.m_sec)
    }
}

impl IMRPhenomXAS {
    /// [`Waveform::waveform`] written into `out` (of the length of
    /// `frequencies`).
    pub(crate) fn waveform_into(&self, frequencies: &[f64], phi_c: f64, out: &mut [Complex<f64>]) {
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return;
        };
        if let Some(tidal) = &self.tidal {
            return self.tidal_h22_into(tidal, frequencies, phi_c, ylm_factor(), out);
        }
        let (phif_ref, mf_cut) = self.grid_phase_and_cutoff(f_min, f_max, phi_c);
        let ylm = ylm_factor();
        for (h, &f) in out.iter_mut().zip(frequencies) {
            *h = self.h22(self.wf.m_sec * f, phif_ref, mf_cut) * ylm;
        }
    }

    /// The (2,2) mode `h_{2,-2}(f)` of `XLALSimIMRPhenomXASFrequencySequence`
    /// (without the spherical-harmonic factor) at `frequencies`, written into
    /// `out`, for the orbital reference phase `phi_ref`.
    pub(crate) fn h22_into(&self, frequencies: &[f64], phi_ref: f64, out: &mut [Complex<f64>]) {
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return;
        };
        if let Some(tidal) = &self.tidal {
            let one = Complex::new(1.0, 0.0);
            return self.tidal_h22_into(tidal, frequencies, phi_ref, one, out);
        }
        let (phif_ref, mf_cut) = self.grid_phase_and_cutoff(f_min, f_max, phi_ref);
        for (h, &f) in out.iter_mut().zip(frequencies) {
            *h = self.h22(self.wf.m_sec * f, phif_ref, mf_cut);
        }
    }
}

impl Waveform for IMRPhenomXAS {
    /// The dimensionless frequency `Mf`.
    fn orbital_speed(&self, frequency: f64) -> f64 {
        self.wf.m_sec * frequency
    }

    /// `pi - phi_22(Mf)` so that `amplitude * exp(-i phase)` is `Y h22`,
    /// for the orbital reference phase `phi_c`.
    fn phase(&self, v: f64, phi_c: f64) -> f64 {
        let powers = UsefulPowers::fast(v);
        if let Some(tidal) = &self.tidal {
            let model = self.for_grid(f64::INFINITY);
            let series = tidal.series(&[v / self.wf.m_sec]);
            let phi = model.phase_22_full(v, &powers, model.phif_ref(2.0 * phi_c));
            return PI - (phi - (series.phase[0] + tidal.multipole_phase(&powers)));
        }
        PI - self.phase_22_full(v, &powers, self.phif_ref(2.0 * phi_c))
    }

    /// `|Y| |h22(Mf)|`, zero above the cutoff frequency (with the NRTidal
    /// corrections, `|Y| (A + A_T) window`, which can be negative).
    fn amplitude(&self, v: f64) -> f64 {
        if let Some(tidal) = &self.tidal {
            if v > self.mf_cut {
                return 0.0;
            }
            let powers = UsefulPowers::fast(v);
            let series = tidal.series(&[v / self.wf.m_sec]);
            let amplitude = self.wf.amp0
                * (self.wf.amp_norm * powers.m_seven_sixths * self.amp22.ansatz(v, &powers)
                    + 2. * (PI / 5.).sqrt() * series.amplitude[0]);
            return 2.0 * (5.0 / (64.0 * PI)).sqrt() * amplitude * series.taper[0];
        }
        if v <= self.mf_cut {
            let powers = UsefulPowers::fast(v);
            2.0 * (5.0 / (64.0 * PI)).sqrt() * self.amplitude_22_full(v, &powers)
        } else {
            0.0
        }
    }

    fn waveform_single_frequency(&self, frequency: f64, phi_c: f64) -> Complex<f64> {
        let mf = self.orbital_speed(frequency);
        if let Some(tidal) = &self.tidal {
            return self.tidal_h22_single(tidal, mf, 2.0 * phi_c) * ylm_factor();
        }
        self.h22(mf, self.phif_ref(2.0 * phi_c), self.mf_cut) * ylm_factor()
    }

    /// `Y h22` on the frequency grid, reproducing
    /// `XLALSimIMRPhenomXASFrequencySequence`: the first frequency is the
    /// default reference frequency and frequencies above
    /// `min(last frequency, f_cut)` are zero.
    fn waveform(&self, frequencies: &[f64], phi_c: f64) -> Vec<Complex<f64>> {
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return Vec::new();
        };
        if let Some(tidal) = &self.tidal {
            let mut out = vec![Complex::new(0.0, 0.0); frequencies.len()];
            self.tidal_h22_into(tidal, frequencies, phi_c, ylm_factor(), &mut out);
            return out;
        }
        let (phif_ref, mf_cut) = self.grid_phase_and_cutoff(f_min, f_max, phi_c);
        let ylm = ylm_factor();
        frequencies
            .iter()
            .map(|&f| self.h22(self.wf.m_sec * f, phif_ref, mf_cut) * ylm)
            .collect()
    }

    fn waveform_on_grid(&self, grid: &FrequencyGrid, phi_c: f64, out: &mut [Complex<f64>]) {
        assert_eq!(out.len(), grid.len(), "output length must match the grid");
        let frequencies = grid.frequencies();
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return;
        };
        if let Some(tidal) = &self.tidal {
            return self.tidal_h22_into(tidal, frequencies, phi_c, ylm_factor(), out);
        }
        let (phif_ref, mf_cut) = self.grid_phase_and_cutoff(f_min, f_max, phi_c);
        // Y h22 = |Y| A exp(i (phi + arg Y)), stored as (amplitude, -phase)
        // for the vectorised conversion.
        let ylm = ylm_factor();
        let (ylm_abs, ylm_arg) = (ylm.norm(), ylm.arg());
        let scale = ScalePowers::new(self.wf.m_sec);
        for (index, h) in out.iter_mut().enumerate() {
            let powers = grid.powers(&scale, index);
            let mf = powers.itself;
            *h = if mf <= mf_cut {
                let powers = UsefulPowers::from_grid(&powers);
                let phi = self.phase_22_full(mf, &powers, phif_ref);
                let amp = self.amplitude_22_full(mf, &powers);
                Complex::new(amp * ylm_abs, -(phi + ylm_arg))
            } else {
                Complex::new(0.0, 0.0)
            };
        }
        fastmath::polar_to_complex_conj(out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{MSUN_SI, MTSUN_SI};
    use crate::imrphenomx::phase::dphase_22;

    /// `XLALSimIMRPhenomXASDuration`: SPA duration from `f_start` to the peak
    /// of the ringdown amplitude, a check of the phase derivative and of the
    /// amplitude coefficients.
    fn duration(mass_1: f64, mass_2: f64, chi_1: f64, chi_2: f64, f_start: f64) -> f64 {
        let xas = IMRPhenomXAS::from_components(mass_1, mass_2, chi_1, chi_2, 1.0).unwrap();
        let m_sec = MTSUN_SI * (mass_1 * MSUN_SI + mass_2 * MSUN_SI) / MSUN_SI;
        let dphi = |mf: f64| xas.inveta * dphase_22(mf, &UsefulPowers::new(mf), &xas.phase22);
        let dphi_start = dphi(f_start * m_sec);
        let dphi_end = dphi(xas.amp22.f_amp_rd_min);
        (dphi_start - dphi_end).abs() / 2.0 / PI * m_sec
    }

    #[test]
    fn duration_matches_lalsimulation() {
        // SimIMRPhenomXASDuration from the lalsimulation 6.2.0 python bindings.
        let cases = [
            (36.0, 29.0, 0.3, -0.4, 20.0, 0.7433398221006442),
            (60.0, 15.0, 0.8, 0.5, 15.0, 2.409389779639988),
            (1.6, 1.3, 0.05, 0.02, 30.0, 52.10480847992271),
            (8.0, 1.4, 0.9, 0.0, 20.0, 43.219328755376466),
            (29.0, 36.0, 0.3, -0.4, 10.0, 5.128074308176688),
        ];
        for (m1, m2, chi1, chi2, f_start, expected) in cases {
            let value = duration(m1, m2, chi1, chi2, f_start);
            assert!(
                ((value - expected) / expected).abs() < 1e-13,
                "{value} != {expected}"
            );
        }
    }

    #[test]
    fn single_frequency_methods_match_grid() {
        let xas =
            IMRPhenomXAS::new(65.0, 29.0 / 36.0, 0.3, -0.4, 500.0).with_reference_frequency(30.0);
        let frequencies: Vec<f64> = (0..200).map(|i| 20.0 + 5.0 * i as f64).collect();
        let grid = xas.waveform(&frequencies, 2.4);
        for (&f, h) in frequencies.iter().zip(grid) {
            assert_eq!(xas.waveform_single_frequency(f, 2.4), h);
            let mf = xas.orbital_speed(f);
            let polar = xas.amplitude(mf) * Complex::new(0.0, -xas.phase(mf, 2.4)).exp();
            assert!((polar - h).norm() <= 1e-12 * h.norm(), "{polar} != {h}");
        }
    }

    #[test]
    fn default_reference_frequency_is_first_frequency() {
        let xas = IMRPhenomXAS::new(40.0, 0.25, -0.2, 0.6, 300.0);
        let frequencies: Vec<f64> = (0..100).map(|i| 17.0 + 3.0 * i as f64).collect();
        let default = xas.waveform(&frequencies, 0.8);
        let explicit = xas
            .with_reference_frequency(17.0)
            .waveform(&frequencies, 0.8);
        assert_eq!(default, explicit);
    }

    #[test]
    fn zero_above_cutoff() {
        let xas = IMRPhenomXAS::new(90.0, 0.5, 0.0, 0.0, 100.0).with_reference_frequency(20.0);
        let f_cut = 0.3 / (90.0 * MTSUN_SI);
        let h = xas.waveform(&[20.0, 0.99 * f_cut, 1.01 * f_cut], 0.0);
        assert!(h[1].norm() > 0.0);
        assert_eq!(h[2], Complex::new(0.0, 0.0));
        assert_eq!(xas.amplitude(xas.orbital_speed(1.01 * f_cut)), 0.0);
    }

    #[test]
    fn rejects_parameters_lal_rejects() {
        assert_eq!(
            IMRPhenomXAS::from_components(30.0, 20.0, 1.2, 0.0, 100.0).unwrap_err(),
            IMRPhenomXError::UnphysicalSpin(1.2)
        );
        assert!(matches!(
            IMRPhenomXAS::from_components(2000.0, 1.0, 0.0, 0.0, 100.0),
            Err(IMRPhenomXError::MassRatioTooLarge(_))
        ));
        assert_eq!(
            IMRPhenomXAS::from_components(-1.0, 1.0, 0.0, 0.0, 100.0).unwrap_err(),
            IMRPhenomXError::NonPositiveMass
        );
    }

    #[test]
    fn masses_are_ordered_like_lal() {
        let a = IMRPhenomXAS::from_components(36.0, 29.0, 0.3, -0.4, 500.0).unwrap();
        let b = IMRPhenomXAS::from_components(29.0, 36.0, -0.4, 0.3, 500.0).unwrap();
        let frequencies = [20.0, 50.0, 150.0, 400.0];
        assert_eq!(a.waveform(&frequencies, 1.0), b.waveform(&frequencies, 1.0));
    }
}
