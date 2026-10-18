//! [`FrequencyDomainModel`] for the frequency-domain models of the crate.

use num_complex::Complex;

use crate::error::{Error, Result};
use crate::grid::FrequencyGrid;
use crate::imrphenomd::IMRPhenomD;
use crate::imrphenomx::{
    IMRPhenomXAS, IMRPhenomXError, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, Polarizations,
    IMRPHENOMXHM_MODES,
};
use crate::model::{Approximant, Binary, Extrinsic, FrequencyDomainModel, WaveformOptions};
use crate::taylor::TaylorF2;
use crate::waveform::Waveform;

/// Panic unless both outputs have length `n`.
fn check_lengths(n: usize, hplus: &[Complex<f64>], hcross: &[Complex<f64>]) {
    assert_eq!(hplus.len(), n, "output length must match the frequencies");
    assert_eq!(hcross.len(), n, "output length must match the frequencies");
}

/// `checkTransverseSpinsZero` of `SimInspiralChooseFDWaveform`.
fn check_aligned(approximant: Approximant, binary: &Binary) -> Result<()> {
    if binary.has_aligned_spins() {
        Ok(())
    } else {
        Err(Error::InPlaneSpins(approximant))
    }
}

/// The masses (solar masses) and aligned spins with the heavier body first.
fn heavier_first(binary: &Binary) -> (f64, f64, f64, f64) {
    let (m1, m2) = (binary.mass_1, binary.mass_2);
    let (chi1, chi2) = (binary.spin_1[2], binary.spin_2[2]);
    if m1 >= m2 {
        (m1, m2, chi1, chi2)
    } else {
        (m2, m1, chi2, chi1)
    }
}

/// `(h_+, h_x)` from the face-on `h_+` in `hplus`, in place: `h_+ = (1 +
/// cos^2 i) / 2 h`, `h_x = -i cos i h`, exactly as
/// [`Waveform::waveform_modes`] (and `SimInspiralChooseFDWaveform`).
fn orient(inclination: f64, hplus: &mut [Complex<f64>], hcross: &mut [Complex<f64>]) {
    let cos_theta_jn = inclination.cos();
    let plus_factor = (1.0 + cos_theta_jn * cos_theta_jn) / 2.0;
    for (p, c) in hplus.iter_mut().zip(hcross.iter_mut()) {
        *c = Complex::new(p.im * cos_theta_jn, -p.re * cos_theta_jn);
        *p *= plus_factor;
    }
}

/// `x > 0` (false for NaN).
fn is_positive(x: f64) -> bool {
    x > 0.0
}

/// `x >= 0` (false for NaN).
fn is_non_negative(x: f64) -> bool {
    x >= 0.0
}

fn invalid(approximant: Approximant, message: &'static str) -> Error {
    Error::InvalidParameter {
        approximant,
        message,
    }
}

/// TaylorF2 as `SimInspiralChooseFDWaveformSequence`: the default (3.5PN)
/// phasing, the Newtonian amplitude and `f_ref = 0` for the phase at
/// coalescence. The masses are not checked against the Kerr bound or for
/// the ISCO, as in LAL.
impl FrequencyDomainModel for TaylorF2 {
    fn from_binary(binary: &Binary) -> Result<Self> {
        let approximant = Approximant::TaylorF2;
        check_aligned(approximant, binary)?;
        if !(binary.mass_1 > 0.0 && binary.mass_2 > 0.0) {
            return Err(invalid(approximant, "masses must be positive"));
        }
        if !is_positive(binary.luminosity_distance) {
            return Err(invalid(approximant, "distance must be positive"));
        }
        if !is_non_negative(binary.f_ref) {
            return Err(invalid(approximant, "f_ref must be non-negative"));
        }
        let (m1, m2, chi1, chi2) = heavier_first(binary);
        let mut model = TaylorF2::new(m1 + m2, m2 / m1, chi1, chi2, binary.luminosity_distance)
            .with_reference_frequency(binary.f_ref);
        model.extrinsic = Extrinsic::of(binary);
        Ok(model)
    }

    fn approximant(&self) -> Approximant {
        Approximant::TaylorF2
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(frequencies.len(), hplus, hcross);
        self.waveform_into(frequencies, self.extrinsic.phi_ref, hplus);
        orient(self.extrinsic.inclination, hplus, hcross);
        Ok(())
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(grid.len(), hplus, hcross);
        let Extrinsic {
            inclination,
            phi_ref,
        } = self.extrinsic;
        self.waveform_modes_on_grid(grid, phi_ref, inclination, hplus, hcross);
        Ok(())
    }
}

impl IMRPhenomD {
    /// The model with the phase referenced as LAL does for a grid starting
    /// at `f_min` (only differs from `self` for `f_ref = 0`).
    fn referenced_to(&self, f_min: Option<f64>) -> Self {
        match f_min {
            Some(f_min) if self.reference_at_first_frequency => {
                let mut model = self.with_reference_frequency(f_min);
                model.reference_at_first_frequency = true;
                model
            }
            _ => *self,
        }
    }
}

/// IMRPhenomD as `SimInspiralChooseFDWaveformSequence`
/// (`XLALSimIMRPhenomDFrequencySequence`), with `f_ref = 0` for the first
/// frequency of the grid.
impl FrequencyDomainModel for IMRPhenomD {
    fn from_binary(binary: &Binary) -> Result<Self> {
        let approximant = Approximant::IMRPhenomD;
        check_aligned(approximant, binary)?;
        if !is_non_negative(binary.f_ref) {
            return Err(invalid(
                approximant,
                "fRef_in must be positive (or 0 for 'ignore')",
            ));
        }
        if !(binary.mass_1 > 0.0 && binary.mass_2 > 0.0) {
            return Err(invalid(approximant, "masses must be positive"));
        }
        if !is_positive(binary.luminosity_distance) {
            return Err(invalid(approximant, "distance must be positive"));
        }
        let (m1, m2, chi1, chi2) = heavier_first(binary);
        if !((-1.0..=1.0).contains(&chi1) && (-1.0..=1.0).contains(&chi2)) {
            return Err(invalid(
                approximant,
                "spins outside the range [-1,1] are not supported",
            ));
        }
        let mut model = IMRPhenomD::new(m1 + m2, m2 / m1, chi1, chi2, binary.luminosity_distance);
        if binary.f_ref > 0.0 {
            model = model.with_reference_frequency(binary.f_ref);
        } else {
            model.reference_at_first_frequency = true;
        }
        model.extrinsic = Extrinsic::of(binary);
        Ok(model)
    }

    fn approximant(&self) -> Approximant {
        Approximant::IMRPhenomD
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(frequencies.len(), hplus, hcross);
        let model = self.referenced_to(frequencies.first().copied());
        model.waveform_into(frequencies, self.extrinsic.phi_ref, hplus);
        orient(self.extrinsic.inclination, hplus, hcross);
        Ok(())
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(grid.len(), hplus, hcross);
        let Extrinsic {
            inclination,
            phi_ref,
        } = self.extrinsic;
        self.referenced_to(grid.frequencies().first().copied())
            .waveform_modes_on_grid(grid, phi_ref, inclination, hplus, hcross);
        Ok(())
    }
}

/// IMRPhenomXAS as `SimInspiralChooseFDWaveformSequence`
/// (`XLALSimIMRPhenomXASFrequencySequence`).
impl FrequencyDomainModel for IMRPhenomXAS {
    fn from_binary(binary: &Binary) -> Result<Self> {
        check_aligned(Approximant::IMRPhenomXAS, binary)?;
        if !is_non_negative(binary.f_ref) {
            return Err(IMRPhenomXError::NegativeReferenceFrequency.into());
        }
        let mut model = IMRPhenomXAS::from_components(
            binary.mass_1,
            binary.mass_2,
            binary.spin_1[2],
            binary.spin_2[2],
            binary.luminosity_distance,
        )?
        .with_reference_frequency(binary.f_ref);
        model.extrinsic = Extrinsic::of(binary);
        Ok(model)
    }

    fn approximant(&self) -> Approximant {
        Approximant::IMRPhenomXAS
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(frequencies.len(), hplus, hcross);
        self.waveform_into(frequencies, self.extrinsic.phi_ref, hplus);
        orient(self.extrinsic.inclination, hplus, hcross);
        Ok(())
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(grid.len(), hplus, hcross);
        let Extrinsic {
            inclination,
            phi_ref,
        } = self.extrinsic;
        self.waveform_modes_on_grid(grid, phi_ref, inclination, hplus, hcross);
        Ok(())
    }
}

/// `IMRPhenomXHM::with_mode_array` returning an error for modes the model
/// does not have.
pub(super) fn xhm_with_mode_array(
    model: IMRPhenomXHM,
    modes: &[(u32, i32)],
) -> Result<IMRPhenomXHM> {
    let valid = |&(ell, emm): &(u32, i32)| IMRPHENOMXHM_MODES.contains(&(ell, emm.unsigned_abs()));
    if modes.is_empty() || !modes.iter().all(valid) {
        return Err(IMRPhenomXError::InvalidModeArray.into());
    }
    Ok(model.with_mode_array(modes))
}

/// IMRPhenomXHM as `SimInspiralChooseFDWaveformSequence`
/// (`XLALSimIMRPhenomXHMFrequencySequence`) with
/// `PhenomXHMThresholdMband = 0`.
impl FrequencyDomainModel for IMRPhenomXHM {
    fn from_binary(binary: &Binary) -> Result<Self> {
        check_aligned(Approximant::IMRPhenomXHM, binary)?;
        if !is_non_negative(binary.f_ref) {
            return Err(IMRPhenomXError::NegativeReferenceFrequency.into());
        }
        let mut model = IMRPhenomXHM::from_components(
            binary.mass_1,
            binary.mass_2,
            binary.spin_1[2],
            binary.spin_2[2],
            binary.luminosity_distance,
        )?
        .with_reference_frequency(binary.f_ref);
        model.extrinsic = Extrinsic::of(binary);
        Ok(model)
    }

    fn approximant(&self) -> Approximant {
        Approximant::IMRPhenomXHM
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(frequencies.len(), hplus, hcross);
        let Extrinsic {
            inclination,
            phi_ref,
        } = self.extrinsic;
        IMRPhenomXHM::polarizations_into(self, frequencies, phi_ref, inclination, hplus, hcross);
        Ok(())
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(grid.len(), hplus, hcross);
        let Extrinsic {
            inclination,
            phi_ref,
        } = self.extrinsic;
        IMRPhenomXHM::polarizations_on_grid(self, grid, phi_ref, inclination, hplus, hcross);
        Ok(())
    }
}

/// IMRPhenomXP for [`crate::model::choose_fd_model`].
pub(super) fn xp(binary: &Binary, options: &WaveformOptions) -> Result<IMRPhenomXP> {
    Ok(IMRPhenomXP::with_version(
        binary,
        options.precession_version,
    )?)
}

/// IMRPhenomXP as `SimInspiralChooseFDWaveformSequence`
/// (`XLALSimIMRPhenomXPFrequencySequence`) with the crate's default Euler
/// angles ([`crate::imrphenomx::PrecessionVersion::SpinTaylor`]).
impl FrequencyDomainModel for IMRPhenomXP {
    fn from_binary(binary: &Binary) -> Result<Self> {
        Ok(IMRPhenomXP::new(binary)?)
    }

    fn approximant(&self) -> Approximant {
        Approximant::IMRPhenomXP
    }

    fn polarizations(&self, frequencies: &[f64]) -> Result<Polarizations> {
        Ok(IMRPhenomXP::polarizations(self, frequencies)?)
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(frequencies.len(), hplus, hcross);
        let (hp, hc) = IMRPhenomXP::polarizations(self, frequencies)?;
        hplus.copy_from_slice(&hp);
        hcross.copy_from_slice(&hc);
        Ok(())
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(grid.len(), hplus, hcross);
        Ok(IMRPhenomXP::polarizations_on_grid(
            self, grid, hplus, hcross,
        )?)
    }
}

/// IMRPhenomXPHM for [`crate::model::choose_fd_model`].
pub(super) fn xphm(binary: &Binary, options: &WaveformOptions) -> Result<IMRPhenomXPHM> {
    let mut model = IMRPhenomXPHM::with_version(binary, options.precession_version)?;
    if let Some(modes) = &options.mode_array {
        // LAL twists up both signs of m' whichever sign is listed
        let mut positive: Vec<(u32, u32)> = modes
            .iter()
            .map(|&(ell, emm)| (ell, emm.unsigned_abs()))
            .collect();
        positive.sort_unstable();
        positive.dedup();
        model = model.with_mode_array(&positive)?;
    }
    Ok(model)
}

/// IMRPhenomXPHM as `SimInspiralChooseFDWaveformSequence`
/// (`XLALSimIMRPhenomXPHMFrequencySequence`) with
/// `PhenomXHMThresholdMband = PhenomXPHMThresholdMband = 0` and the crate's
/// default Euler angles ([`crate::imrphenomx::PrecessionVersion::SpinTaylor`]).
impl FrequencyDomainModel for IMRPhenomXPHM {
    fn from_binary(binary: &Binary) -> Result<Self> {
        Ok(IMRPhenomXPHM::new(binary)?)
    }

    fn approximant(&self) -> Approximant {
        Approximant::IMRPhenomXPHM
    }

    fn polarizations(&self, frequencies: &[f64]) -> Result<Polarizations> {
        Ok(IMRPhenomXPHM::polarizations(self, frequencies)?)
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(frequencies.len(), hplus, hcross);
        let (hp, hc) = IMRPhenomXPHM::polarizations(self, frequencies)?;
        hplus.copy_from_slice(&hp);
        hcross.copy_from_slice(&hc);
        Ok(())
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<()> {
        check_lengths(grid.len(), hplus, hcross);
        Ok(IMRPhenomXPHM::polarizations_on_grid(
            self, grid, hplus, hcross,
        )?)
    }
}
