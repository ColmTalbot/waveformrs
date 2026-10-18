//! The aligned-spin higher-mode model IMRPhenomXHM, following
//! `IMRPhenomXHMGenerateFDOneMode`, `IMRPhenomXHM_MultiMode2` and
//! `IMRPhenomXHMFDAddMode` in `LALSimIMRPhenomXHM.c` for non-uniform
//! frequency grids (no multibanding).

use std::f64::consts::FRAC_PI_2;

use num_complex::Complex;

use crate::constants::{MPC_SI, MSUN_SI};
use crate::fastmath;
use crate::grid::{FrequencyGrid, ScalePowers};
use crate::imrphenomx::hm_internals::HMMode;
use crate::imrphenomx::internals::{IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::phase::PhaseCoefficients;
use crate::imrphenomx::wigner::SphericalHarmonics;
use crate::imrphenomx::xas::IMRPhenomXAS;

/// The modes of IMRPhenomXHM, `(l, m)` with `m > 0`.
pub const IMRPHENOMXHM_MODES: [(u32, u32); 5] = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4)];

/// The default mode array of IMRPhenomXHM: all modes with both signs of m.
const DEFAULT_MODE_ARRAY: [(u32, i32); 10] = [
    (2, 2),
    (2, -2),
    (2, 1),
    (2, -1),
    (3, 3),
    (3, -3),
    (3, 2),
    (3, -2),
    (4, 4),
    (4, -4),
];

/// IMRPhenomXHM (García-Quirós et al. 2020, arXiv:2001.10914, with the
/// amplitude recalibration of the 122022 release), the aligned-spin model
/// with the (2,|2|), (2,|1|), (3,|3|), (3,|2|) (including spheroidal-spherical
/// mode mixing) and (4,|4|) modes, for the default LALSimulation settings
/// with multibanding switched off (`PhenomXHMThresholdMband = 0`).
///
/// # Conventions
///
/// Masses are in solar masses with `mass_ratio = m2 / m1 <= 1`, spins are
/// the dimensionless components along the orbital angular momentum and the
/// distance is in Mpc, as for [`IMRPhenomXAS`].
///
/// [`IMRPhenomXHM::polarizations`] reproduces the `hplus`, `hcross` of
/// `SimInspiralChooseFDWaveformSequence(phi_ref, m1, m2, 0, 0, chi_1, 0, 0,
/// chi_2, f_ref, distance, inclination, params, IMRPhenomXHM, frequencies)`
/// with `PhenomXHMThresholdMband = 0`: `phi_ref` is LAL's orbital reference
/// phase and `inclination` the angle between the orbital angular momentum
/// and the line of sight. Following LAL, each mode enters with the
/// spin-weighted spherical harmonic `Y_{-2}^{lm}(inclination, pi / 2)`,
///
/// `hplus = sum_{l, m > 0} [Y^{l,-m} + (-1)^l conj(Y^{lm})] / 2 * h_{l,-m}`,
/// `hcross = sum_{l, m > 0} i [Y^{l,-m} - (-1)^l conj(Y^{lm})] / 2 * h_{l,-m}`,
///
/// where `h_{l,-m}(f)` is the mode returned by [`IMRPhenomXHM::mode`] for
/// negative `m`. The (3,2) mode reuses the (2,2) mode in its ringdown when
/// the latter is in the mode array, as LAL does.
///
/// # Reference frequency
///
/// [`IMRPhenomXHM::with_reference_frequency`] sets LAL's `f_ref`; if it is
/// not set (or zero) the first frequency of the grid is used.
#[derive(Clone, Debug)]
pub struct IMRPhenomXHM {
    /// Component masses in kg, as received by LAL.
    m1_si: f64,
    m2_si: f64,
    chi_1: f64,
    chi_2: f64,
    /// Luminosity distance in m.
    distance: f64,
    /// Reference frequency in Hz, zero for the first frequency of the grid.
    f_ref: f64,
    /// The active modes `(l, m)`.
    mode_array: Vec<(u32, i32)>,
}

impl IMRPhenomXHM {
    /// Construct the model from the total mass (solar masses), the mass
    /// ratio `m2 / m1 <= 1`, the aligned dimensionless spins and the
    /// luminosity distance (Mpc).
    ///
    /// # Panics
    ///
    /// If LALSimulation would refuse the parameters (see
    /// [`IMRPhenomXHM::from_components`]).
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
            .unwrap_or_else(|e| panic!("IMRPhenomXHM: {e}"))
    }

    /// Construct the model from the component masses (solar masses, in any
    /// order), the aligned spins and the luminosity distance (Mpc).
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
        let new = Self {
            m1_si: mass_1 * MSUN_SI,
            m2_si: mass_2 * MSUN_SI,
            chi_1,
            chi_2,
            distance: luminosity_distance * MPC_SI,
            f_ref: 0.0,
            mode_array: DEFAULT_MODE_ARRAY.to_vec(),
        };
        // Check the parameters as LAL does when the waveform is generated.
        new.waveform_struct(10.0, 0.0, 0.0)?;
        Ok(new)
    }

    /// Reference the phase to `reference_frequency` (Hz), LAL's `f_ref`.
    /// Zero selects the first frequency of the grid, as in LALSimulation.
    pub fn with_reference_frequency(mut self, reference_frequency: f64) -> Self {
        assert!(
            reference_frequency >= 0.0,
            "IMRPhenomXHM: {}",
            IMRPhenomXError::NegativeReferenceFrequency
        );
        self.f_ref = reference_frequency;
        self
    }

    /// Restrict the waveform to the modes `(l, m)` (LAL's mode array). Both
    /// signs of `m` can be given; a mode with only one sign active
    /// contributes only that sign.
    ///
    /// # Panics
    ///
    /// If a mode is not part of IMRPhenomXHM.
    pub fn with_mode_array(mut self, modes: &[(u32, i32)]) -> Self {
        for &(l, m) in modes {
            assert!(
                IMRPHENOMXHM_MODES.contains(&(l, m.unsigned_abs())),
                "IMRPhenomXHM: mode ({l},{m}) is not available"
            );
        }
        self.mode_array = modes.to_vec();
        self
    }

    /// Total mass in solar masses.
    pub fn total_mass(&self) -> f64 {
        self.m1_si / MSUN_SI + self.m2_si / MSUN_SI
    }

    fn is_active(&self, l: u32, m: i32) -> bool {
        self.mode_array.contains(&(l, m))
    }

    /// `IMRPhenomXSetWaveformVariables` for a grid from `f_min` to `f_max`.
    fn waveform_struct(
        &self,
        f_min: f64,
        f_max: f64,
        phi_ref: f64,
    ) -> Result<WaveformStruct, IMRPhenomXError> {
        WaveformStruct::from_si(
            self.m1_si,
            self.m2_si,
            self.chi_1,
            self.chi_2,
            0.0,
            self.f_ref,
            phi_ref,
            f_min,
            f_max,
            self.distance,
            0.0,
        )
    }

    /// The polarizations `(hplus, hcross)` at `frequencies` (Hz, increasing)
    /// for the orbital reference phase `phi_ref` and the inclination
    /// `inclination` (radians), reproducing
    /// `XLALSimIMRPhenomXHMFrequencySequence`.
    ///
    /// Like LALSimulation, the (2,2) mode vanishes above
    /// `min(last frequency, f_cut)` with `f_cut` at `Mf = 0.3` (0.33 for
    /// `chi_eff > 0.99`), while the higher modes are evaluated at every
    /// frequency of the grid.
    pub fn polarizations(
        &self,
        frequencies: &[f64],
        phi_ref: f64,
        inclination: f64,
    ) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
        let zero = Complex::new(0.0, 0.0);
        let n = frequencies.len();
        let mut hplus = vec![zero; n];
        let mut hcross = vec![zero; n];
        // The powers as LAL computes them: with these the phases of long
        // inspirals (~1e6 rad) agree with LAL to the last bit.
        self.fill_polarizations(
            frequencies,
            phi_ref,
            inclination,
            |_, mf| UsefulPowers::new(mf),
            &mut hplus,
            &mut hcross,
        );
        (hplus, hcross)
    }

    /// [`Self::polarizations`] on a prepared [`FrequencyGrid`], written into
    /// `hplus` and `hcross` (each with the length of the grid).
    ///
    /// For evaluating many waveforms on the same frequencies: the powers of
    /// `Mf` are products of factors stored in the grid rather than LAL's
    /// `pow` calls, and nothing is allocated per frequency. The result agrees
    /// with [`Self::polarizations`] to rounding (~1e-13 relative, more for
    /// inspirals with phases of order 1e6 rad).
    pub fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        phi_ref: f64,
        inclination: f64,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) {
        assert_eq!(hplus.len(), grid.len(), "output length must match the grid");
        assert_eq!(
            hcross.len(),
            grid.len(),
            "output length must match the grid"
        );
        self.fill_polarizations_grid(grid, phi_ref, inclination, hplus, hcross);
    }

    /// The waveform variables, the (2,2) and higher modes and their
    /// coefficients in the polarizations, for a grid from `f_min` to `f_max`.
    fn mode_sum(&self, f_min: f64, f_max: f64, phi_ref: f64, inclination: f64) -> ModeSum {
        let wf = self
            .waveform_struct(f_min, f_max, phi_ref)
            .unwrap_or_else(|e| panic!("IMRPhenomXHM: {e}"));

        // The (2,2) mode from XLALSimIMRPhenomXASFrequencySequence, which
        // sets up its own waveform variables from the LAL ones.
        let pos22 = self.is_active(2, 2);
        let neg22 = self.is_active(2, -2);
        let xas = if pos22 || neg22 {
            let xas_wf = WaveformStruct::from_si(
                wf.m1_si,
                wf.m2_si,
                wf.chi1l,
                wf.chi2l,
                0.0,
                wf.f_ref,
                wf.phi_ref_in,
                f_min,
                f_max,
                wf.distance,
                0.0,
            )
            .unwrap_or_else(|e| panic!("IMRPhenomXHM: {e}"));
            Some(Xas22::new(xas_wf))
        } else {
            None
        };

        // Higher modes, in the order of the mode loop of LAL.
        let phase22 = PhaseCoefficients::new(&wf);
        let mut modes = Vec::new();
        for ell in 2..=4u32 {
            for emm in 1..=ell {
                if (ell, emm) == (2, 2) {
                    continue;
                }
                let pos = self.is_active(ell, emm as i32);
                let neg = self.is_active(ell, -(emm as i32));
                if !pos && !neg {
                    continue;
                }
                let mode = HMMode::new(ell, emm, &wf, &phase22);
                if mode.wfhm.amp_zero {
                    continue;
                }
                let amp0 = minus_one_to(ell) * wf.amp0;
                let (factor_plus, factor_cross) = add_mode_factors(inclination, ell, emm, pos, neg);
                modes.push((mode, amp0, factor_plus, factor_cross));
            }
        }
        ModeSum {
            wf,
            xas,
            modes,
            factors22: add_mode_factors(inclination, 2, 2, pos22, neg22),
        }
    }

    /// The mode sum of `IMRPhenomXHM_MultiMode` on `frequencies`, with the
    /// powers of `Mf` at frequency `i` given by `powers_at(i, Mf)`.
    fn fill_polarizations(
        &self,
        frequencies: &[f64],
        phi_ref: f64,
        inclination: f64,
        powers_at: impl Fn(usize, f64) -> UsefulPowers,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) {
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return;
        };
        let sum = self.mode_sum(f_min, f_max, phi_ref, inclination);
        let zero = Complex::new(0.0, 0.0);
        let (factor_plus_22, factor_cross_22) = sum.factors22;
        for (i, &f) in frequencies.iter().enumerate() {
            let mf = sum.wf.m_sec * f;
            let powers = powers_at(i, mf);
            let h22 = sum.xas.as_ref().map(|xas| xas.evaluate(f, mf, &powers));
            let mut hp = zero;
            let mut hc = zero;
            if let Some(h22) = h22 {
                hp += factor_plus_22 * h22;
                hc += factor_cross_22 * h22;
            }
            for (mode, amp0, factor_plus, factor_cross) in &sum.modes {
                let (amp, phi) = mode.amplitude_phase(&powers, h22);
                let h = Complex::from_polar(amp0 * amp, phi);
                hp += factor_plus * h;
                hc += factor_cross * h;
            }
            hplus[i] = hp;
            hcross[i] = hc;
        }
    }

    /// [`Self::fill_polarizations`] on a grid, one mode at a time so that
    /// the conversion from amplitude and phase vectorises.
    fn fill_polarizations_grid(
        &self,
        grid: &FrequencyGrid,
        phi_ref: f64,
        inclination: f64,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) {
        let frequencies = grid.frequencies();
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return;
        };
        let sum = self.mode_sum(f_min, f_max, phi_ref, inclination);
        let zero = Complex::new(0.0, 0.0);
        hplus.fill(zero);
        hcross.fill(zero);
        let scale = ScalePowers::new(sum.wf.m_sec);
        let accumulate = |h: &[Complex<f64>],
                          factor_plus: Complex<f64>,
                          factor_cross: Complex<f64>,
                          hplus: &mut [Complex<f64>],
                          hcross: &mut [Complex<f64>]| {
            for ((hp, hc), h) in hplus.iter_mut().zip(hcross.iter_mut()).zip(h) {
                *hp += factor_plus * h;
                *hc += factor_cross * h;
            }
        };

        // h22 = A exp(i phi), stored as (A, -phi) for the vectorised conversion
        let h22 = sum.xas.as_ref().map(|xas| {
            let scale22 = ScalePowers::new(xas.m_sec);
            let mut h22: Vec<Complex<f64>> = (0..grid.len())
                .map(|index| {
                    let powers = grid.powers(&scale22, index);
                    let mf = powers.itself;
                    if mf > xas.mf_cut {
                        return zero;
                    }
                    let (amp, phi) = xas.amplitude_phase(mf, &UsefulPowers::from_grid(&powers));
                    Complex::new(amp, -phi)
                })
                .collect();
            fastmath::polar_to_complex_conj(&mut h22);
            h22
        });
        if let Some(h22) = &h22 {
            let (factor_plus, factor_cross) = sum.factors22;
            accumulate(h22, factor_plus, factor_cross, hplus, hcross);
        }

        let mut h = vec![zero; grid.len()];
        for (mode, amp0, factor_plus, factor_cross) in &sum.modes {
            for (index, value) in h.iter_mut().enumerate() {
                let powers = UsefulPowers::from_grid(&grid.powers(&scale, index));
                let (amp, phi) = mode.amplitude_phase(&powers, h22.as_ref().map(|h22| h22[index]));
                *value = Complex::new(amp0 * amp, -phi);
            }
            fastmath::polar_to_complex_conj(&mut h);
            accumulate(&h, *factor_plus, *factor_cross, hplus, hcross);
        }
    }

    /// The mode `h_lm` at `frequencies` (Hz, increasing) for the orbital
    /// reference phase `phi_ref`, reproducing
    /// `XLALSimIMRPhenomXHMFrequencySequenceOneMode`.
    ///
    /// For `m < 0` this is `h_{l,m}(f)`, which lives at positive frequencies.
    /// For `m > 0` LAL returns `(-1)^l conj(h_{l,-m}(f)) = h_{l,m}(-f)`, the
    /// mode at the negative frequency `-f`. The (3,2) mode is computed
    /// without recycling the (2,2) mode, as in LAL's one-mode functions.
    /// Each mode vanishes above `min(last frequency, f_cut)`.
    ///
    /// # Panics
    ///
    /// If `(l, |m|)` is not a mode of IMRPhenomXHM.
    pub fn mode(&self, ell: u32, emm: i32, frequencies: &[f64], phi_ref: f64) -> Vec<Complex<f64>> {
        assert!(
            IMRPHENOMXHM_MODES.contains(&(ell, emm.unsigned_abs())),
            "IMRPhenomXHM: mode ({ell},{emm}) is not available"
        );
        let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
            return Vec::new();
        };
        let wf = self
            .waveform_struct(f_min, f_max, phi_ref)
            .unwrap_or_else(|e| panic!("IMRPhenomXHM: {e}"));
        let mut h = generate_fd_one_mode(&wf, ell, emm.unsigned_abs(), frequencies);
        if emm > 0 {
            // Transform to positive m: (-1)^l conj(h_{l,-m}).
            let sign = if ell == 2 && emm == 2 {
                1.0
            } else {
                minus_one_to(ell)
            };
            for value in &mut h {
                *value = sign * value.conj();
            }
        }
        h
    }
}

/// `(-1)^l`
fn minus_one_to(ell: u32) -> f64 {
    if ell % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

/// Everything the mode sum needs for one waveform.
struct ModeSum {
    wf: WaveformStruct,
    xas: Option<Xas22>,
    modes: Vec<(HMMode, f64, Complex<f64>, Complex<f64>)>,
    factors22: (Complex<f64>, Complex<f64>),
}

/// The (2,2) mode of `XLALSimIMRPhenomXASFrequencySequence`.
struct Xas22 {
    xas: IMRPhenomXAS,
    phif_ref: f64,
    m_sec: f64,
    mf_cut: f64,
}

impl Xas22 {
    fn new(wf: WaveformStruct) -> Self {
        let xas = IMRPhenomXAS::from_waveform_struct(wf, true);
        Self {
            phif_ref: xas.phif_ref(2.0 * wf.phi0),
            m_sec: wf.m_sec,
            mf_cut: wf.f_max_prime * wf.m_sec,
            xas,
        }
    }

    /// Amplitude and phase of the (2,2) mode, `h22 = A exp(i phi)`, below
    /// the cutoff.
    #[inline]
    fn amplitude_phase(&self, mf: f64, powers: &UsefulPowers) -> (f64, f64) {
        (
            self.xas.amplitude_22_full(mf, powers),
            self.xas.phase_22_full(mf, powers, self.phif_ref),
        )
    }

    /// The (2,2) mode at `f`, reusing the powers of `mf` when the
    /// dimensionless frequencies agree.
    #[inline]
    fn evaluate(&self, f: f64, mf: f64, powers: &UsefulPowers) -> Complex<f64> {
        let mf22 = self.m_sec * f;
        if mf22 > self.mf_cut {
            Complex::new(0.0, 0.0)
        } else if mf22 == mf {
            self.xas.h22_with_powers(mf22, powers, self.phif_ref)
        } else {
            self.xas
                .h22_with_powers(mf22, &UsefulPowers::new(mf22), self.phif_ref)
        }
    }
}

/// One mode of IMRPhenomXHM for given waveform variables: the setup part of
/// `IMRPhenomXASGenerateFD` for the (2,2) and of
/// `IMRPhenomXHMGenerateFDOneMode` for the higher modes.
///
/// [`ModeGenerator::h`] returns the mode with negative `m`, `h_{l,-m}(f)`
/// for `f > 0`, in strain units, exactly as these LAL functions do. The
/// waveform variables may be modified first (e.g. the final spin and the
/// reference phase), as the precessing models do.
#[derive(Clone, Debug)]
pub(crate) enum ModeGenerator {
    /// The (2,2) mode of IMRPhenomXAS with its `phifRef`.
    Xas22 {
        xas: Box<IMRPhenomXAS>,
        phif_ref: f64,
    },
    /// A higher mode with `(-1)^l amp0`.
    HigherMode { mode: Box<HMMode>, amp0: f64 },
}

impl ModeGenerator {
    /// Set up mode (`ell`, `emm`), `emm > 0`, from the waveform variables
    /// `wf`. `phase22` must be `PhaseCoefficients::new(wf)` (the 22 phase
    /// coefficients without connection coefficients, as returned by
    /// `IMRPhenomXGetPhaseCoefficients`).
    pub(crate) fn new(
        wf: &WaveformStruct,
        phase22: &PhaseCoefficients,
        ell: u32,
        emm: u32,
    ) -> Self {
        if (ell, emm) == (2, 2) {
            let xas = Box::new(IMRPhenomXAS::from_waveform_struct(*wf, true));
            let phif_ref = xas.phif_ref(2.0 * wf.phi0);
            Self::Xas22 { xas, phif_ref }
        } else {
            Self::HigherMode {
                mode: Box::new(HMMode::new(ell, emm, wf, phase22)),
                amp0: minus_one_to(ell) * wf.amp0,
            }
        }
    }

    /// `h_{l,-m}` at the dimensionless frequency `mf` (no cutoff), with
    /// `powers` the powers of `mf`.
    #[inline]
    pub(crate) fn h(&self, mf: f64, powers: &UsefulPowers) -> Complex<f64> {
        match self {
            Self::Xas22 { xas, phif_ref } => xas.h22_with_powers(mf, powers, *phif_ref),
            Self::HigherMode { mode, amp0 } => {
                if mode.wfhm.amp_zero {
                    return Complex::new(0.0, 0.0);
                }
                let (amp, phi) = mode.amplitude_phase(powers, None);
                Complex::from_polar(amp0 * amp, phi)
            }
        }
    }

    /// Whether the mode vanishes identically (odd modes of equal-mass,
    /// equal-spin binaries).
    pub(crate) fn vanishes(&self) -> bool {
        matches!(self, Self::HigherMode { mode, .. } if mode.wfhm.amp_zero)
    }

    /// Amplitude and phase of [`Self::h`], `h = amplitude exp(i phase)`, for
    /// a mode that does not vanish.
    #[inline]
    pub(crate) fn amplitude_phase(&self, mf: f64, powers: &UsefulPowers) -> (f64, f64) {
        match self {
            Self::Xas22 { xas, phif_ref } => (
                xas.amplitude_22_full(mf, powers),
                xas.phase_22_full(mf, powers, *phif_ref),
            ),
            Self::HigherMode { mode, amp0 } => {
                let (amp, phi) = mode.amplitude_phase(powers, None);
                (amp0 * amp, phi)
            }
        }
    }
}

/// `IMRPhenomXHMGenerateFDOneMode` (and `IMRPhenomXASGenerateFD` for the
/// (2,2)) on a non-uniform grid: `h_{l,-m}(f)` at `frequencies` (Hz) for the
/// waveform variables `wf`, zero above `wf.f_max_prime`.
pub(crate) fn generate_fd_one_mode(
    wf: &WaveformStruct,
    ell: u32,
    emm: u32,
    frequencies: &[f64],
) -> Vec<Complex<f64>> {
    let phase22 = PhaseCoefficients::new(wf);
    let generator = ModeGenerator::new(wf, &phase22, ell, emm);
    let mf_cut = wf.f_max_prime * wf.m_sec;
    frequencies
        .iter()
        .map(|&f| {
            let mf = wf.m_sec * f;
            if mf <= mf_cut {
                generator.h(mf, &UsefulPowers::new(mf))
            } else {
                Complex::new(0.0, 0.0)
            }
        })
        .collect()
}

/// The factors multiplying `h_{l,-m}` in `hplus` and `hcross`,
/// `IMRPhenomXHMFDAddMode` with azimuth `pi / 2`.
fn add_mode_factors(
    inclination: f64,
    ell: u32,
    emm: u32,
    positive: bool,
    negative: bool,
) -> (Complex<f64>, Complex<f64>) {
    let l = ell as i32;
    let m = emm as i32;
    let minus1l = minus_one_to(ell);
    let i = Complex::new(0.0, 1.0);
    if positive && negative {
        let ym = spin_weighted_spherical_harmonic(inclination, FRAC_PI_2, l, -m);
        let ystar = spin_weighted_spherical_harmonic(inclination, FRAC_PI_2, l, m).conj();
        let factor_plus = 0.5 * (ym + minus1l * ystar);
        let factor_cross = i * 0.5 * (ym - minus1l * ystar);
        (factor_plus, factor_cross)
    } else if positive {
        let ylm = spin_weighted_spherical_harmonic(inclination, FRAC_PI_2, l, m).conj();
        (0.5 * minus1l * ylm, -i * 0.5 * minus1l * ylm)
    } else if negative {
        let ylm = spin_weighted_spherical_harmonic(inclination, FRAC_PI_2, l, -m);
        (0.5 * ylm, i * 0.5 * ylm)
    } else {
        (Complex::new(0.0, 0.0), Complex::new(0.0, 0.0))
    }
}

/// `XLALSpinWeightedSphericalHarmonic` for spin weight -2 and `l` = 2, 3, 4:
/// the real harmonic at azimuth zero times `cpolar(1, m phi)`.
pub(crate) fn spin_weighted_spherical_harmonic(
    theta: f64,
    phi: f64,
    l: i32,
    m: i32,
) -> Complex<f64> {
    let harmonics = SphericalHarmonics::new(theta);
    let index = (m + l) as usize;
    let fac = match l {
        2 => harmonics.y2[index],
        3 => harmonics.y3[index],
        4 => harmonics.y4[index],
        _ => panic!("spin-weighted spherical harmonic ({l},{m}) not implemented"),
    };
    if m != 0 {
        let (s, c) = (m as f64 * phi).sin_cos();
        Complex::new(c * fac, s * fac)
    } else {
        Complex::new(fac, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn spherical_harmonics_are_normalised() {
        // sum_m |Y_{-2}^{lm}|^2 = (2l + 1) / (4 pi)
        for l in 2..=4 {
            for &theta in &[0.0, 0.4, 1.3, 2.9] {
                let sum: f64 = (-l..=l)
                    .map(|m| spin_weighted_spherical_harmonic(theta, 0.7, l, m).norm_sqr())
                    .sum();
                let expected = (2 * l + 1) as f64 / (4.0 * PI);
                assert!((sum - expected).abs() < 1e-14, "l={l}: {sum} != {expected}");
            }
        }
    }

    #[test]
    fn odd_modes_vanish_for_equal_masses_and_spins() {
        let xhm = IMRPhenomXHM::new(60.0, 1.0, 0.3, 0.3, 100.0);
        let frequencies: Vec<f64> = (0..50).map(|i| 20.0 + 4.0 * i as f64).collect();
        for (l, m) in [(2, -1), (3, -3)] {
            assert!(xhm
                .mode(l, m, &frequencies, 0.0)
                .iter()
                .all(|h| h.norm() == 0.0));
        }
        assert!(xhm
            .mode(4, -4, &frequencies, 0.0)
            .iter()
            .all(|h| h.norm() > 0.0));
    }

    #[test]
    fn positive_m_is_conjugate_of_negative_m() {
        let xhm = IMRPhenomXHM::new(50.0, 0.4, 0.2, -0.3, 100.0).with_reference_frequency(20.0);
        let frequencies: Vec<f64> = (0..50).map(|i| 20.0 + 5.0 * i as f64).collect();
        for (l, m) in IMRPHENOMXHM_MODES {
            let m = m as i32;
            let neg = xhm.mode(l, -m, &frequencies, 0.3);
            let pos = xhm.mode(l, m, &frequencies, 0.3);
            let sign = if (l, m) == (2, 2) {
                1.0
            } else {
                minus_one_to(l)
            };
            for (p, n) in pos.iter().zip(&neg) {
                assert_eq!(*p, sign * n.conj());
            }
        }
    }
}
