//! Waveform variables and useful powers, the parts of
//! `LALSimIMRPhenomX_internals.c` that do not depend on the mode.

use std::f64::consts::PI;
use std::fmt;

use once_cell::sync::Lazy;

use crate::constants::{MPC_SI, MRSUN_SI, MSUN_SI, MTSUN_SI};
use crate::imrphenomx::qnm::{fdamp22, fring22};
use crate::imrphenomx::utilities::{
    chi_eff, chi_pn_hat, dchi, f_isco, f_meco, final_mass_2017, final_spin_2017, s_tot_r,
};

/// Errors raised for parameters at which LALSimulation refuses to generate
/// IMRPhenomX waveforms.
#[derive(Clone, Debug, PartialEq)]
pub enum IMRPhenomXError {
    /// A component mass is not positive.
    NonPositiveMass,
    /// The mass ratio exceeds 1000.
    MassRatioTooLarge(f64),
    /// An aligned spin is outside [-1, 1].
    UnphysicalSpin(f64),
    /// The luminosity distance is negative.
    NegativeDistance,
    /// The reference frequency is negative.
    NegativeReferenceFrequency,
    /// The fitted final mass exceeds the total mass.
    UnphysicalFinalMass(f64),
    /// The fitted final spin exceeds one in magnitude.
    UnphysicalFinalSpin(f64),
    /// The MECO frequency exceeds the ISCO frequency of the remnant.
    MecoAboveIsco { f_meco: f64, f_isco: f64 },
    /// A spin vector exceeds the Kerr bound, |chi| > 1 (precessing models).
    UnphysicalSpinMagnitude(f64),
    /// The frequency grid is invalid (empty, non-positive or not increasing
    /// where LALSimulation requires it).
    InvalidFrequencies(&'static str),
    /// The numerical SpinTaylor Euler angles could not be computed in a
    /// case where LALSimulation silently falls back to the MSA angles
    /// (`PhenomXPrecVersion` 223), which are not implemented here.
    SpinTaylorFallbackToMsa(&'static str),
    /// The numerical SpinTaylor Euler angles could not be computed in a
    /// case where LALSimulation raises an error.
    SpinTaylorAngles(&'static str),
}

impl fmt::Display for IMRPhenomXError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositiveMass => write!(f, "masses must be positive"),
            Self::MassRatioTooLarge(q) => {
                write!(f, "model not valid at mass ratios beyond 1000 (q = {q})")
            }
            Self::UnphysicalSpin(chi) => write!(
                f,
                "unphysical spin {chi} requested: must obey the Kerr bound [-1, 1]"
            ),
            Self::NegativeDistance => write!(f, "distance must be non-negative"),
            Self::NegativeReferenceFrequency => {
                write!(f, "reference frequency must be positive or zero to ignore")
            }
            Self::UnphysicalFinalMass(m) => write!(f, "final mass {m} > 1 is not physical"),
            Self::UnphysicalFinalSpin(a) => write!(f, "final spin {a} > 1 is not physical"),
            Self::MecoAboveIsco { f_meco, f_isco } => write!(
                f,
                "f_MECO ({f_meco}) cannot be greater than f_ISCO ({f_isco})"
            ),
            Self::UnphysicalSpinMagnitude(chi) => {
                write!(f, "spin magnitude {chi} violates the Kerr bound |chi| <= 1")
            }
            Self::InvalidFrequencies(why) => write!(f, "invalid frequencies: {why}"),
            Self::SpinTaylorFallbackToMsa(why) => write!(
                f,
                "SpinTaylor Euler angles failed ({why}); lalsimulation falls back to the \
                 MSA angles (PrecVersion 223), which are not implemented"
            ),
            Self::SpinTaylorAngles(why) => write!(f, "SpinTaylor Euler angles: {why}"),
        }
    }
}

impl std::error::Error for IMRPhenomXError {}

/// Cached powers of a number, `IMRPhenomX_UsefulPowers`.
///
/// [`UsefulPowers::new`] reproduces `IMRPhenomX_Initialize_Powers`
/// operation by operation (all fractional powers derived from
/// `pow(x, 1/6)`) and is used wherever LAL evaluates something once per
/// waveform. [`UsefulPowers::fast`] computes the same quantities from a
/// single `cbrt`, `sqrt`, `ln` and one division, for use in the
/// per-frequency loop; the two agree to a few units in the last place.
#[derive(Clone, Copy, Debug, Default)]
#[allow(dead_code)] // the full LAL struct is kept for the higher-mode models
pub(crate) struct UsefulPowers {
    pub(crate) one_sixth: f64,
    pub(crate) m_one_sixth: f64,
    pub(crate) itself: f64,
    pub(crate) m_one: f64,
    pub(crate) one_third: f64,
    pub(crate) m_one_third: f64,
    pub(crate) two_thirds: f64,
    pub(crate) m_two_thirds: f64,
    pub(crate) four_thirds: f64,
    pub(crate) m_four_thirds: f64,
    pub(crate) five_thirds: f64,
    pub(crate) m_five_thirds: f64,
    pub(crate) seven_thirds: f64,
    pub(crate) m_seven_thirds: f64,
    pub(crate) eight_thirds: f64,
    pub(crate) m_eight_thirds: f64,
    pub(crate) two: f64,
    pub(crate) three: f64,
    pub(crate) four: f64,
    pub(crate) five: f64,
    pub(crate) m_two: f64,
    pub(crate) m_three: f64,
    pub(crate) m_four: f64,
    pub(crate) m_five: f64,
    pub(crate) seven_sixths: f64,
    pub(crate) m_seven_sixths: f64,
    pub(crate) log: f64,
    pub(crate) sqrt: f64,
}

impl UsefulPowers {
    /// `IMRPhenomX_Initialize_Powers`, reproduced operation by operation.
    pub(crate) fn new(number: f64) -> Self {
        let sixth = number.powf(1.0 / 6.0);
        let m_sixth = 1.0 / sixth;
        let one_third = sixth * sixth;
        let two_thirds = one_third * one_third;
        let four_thirds = two_thirds * two_thirds;
        let five_thirds = four_thirds * one_third;
        let seven_thirds = four_thirds * number;
        let eight_thirds = seven_thirds * one_third;
        let two = number * number;
        let three = two * number;
        let four = three * number;
        let five = four * number;
        let m_one = 1.0 / number;
        Self {
            one_sixth: sixth,
            m_one_sixth: m_sixth,
            itself: number,
            m_one,
            one_third,
            m_one_third: 1.0 / one_third,
            two_thirds,
            m_two_thirds: 1.0 / two_thirds,
            four_thirds,
            m_four_thirds: 1.0 / four_thirds,
            five_thirds,
            m_five_thirds: 1.0 / five_thirds,
            seven_thirds,
            m_seven_thirds: 1.0 / seven_thirds,
            eight_thirds,
            m_eight_thirds: 1.0 / eight_thirds,
            two,
            three,
            four,
            five,
            m_two: 1.0 / two,
            m_three: 1.0 / three,
            m_four: 1.0 / four,
            m_five: 1.0 / five,
            seven_sixths: sixth * number,
            m_seven_sixths: m_sixth * m_one,
            log: number.ln(),
            sqrt: sixth * sixth * sixth,
        }
    }

    /// The same powers as [`UsefulPowers::new`] using one `cbrt`, one
    /// `sqrt`, one `ln` and a single division.
    #[inline]
    pub(crate) fn fast(number: f64) -> Self {
        let one_third = number.cbrt();
        let sixth = one_third.sqrt();
        let m_one = 1.0 / number;
        let two_thirds = one_third * one_third;
        let m_one_third = two_thirds * m_one;
        let m_one_sixth = sixth * m_one_third;
        let four_thirds = two_thirds * two_thirds;
        let five_thirds = four_thirds * one_third;
        let seven_thirds = four_thirds * number;
        let eight_thirds = seven_thirds * one_third;
        let m_two_thirds = m_one_third * m_one_third;
        let m_four_thirds = m_two_thirds * m_two_thirds;
        let m_five_thirds = m_four_thirds * m_one_third;
        let m_seven_thirds = m_four_thirds * m_one;
        let m_eight_thirds = m_seven_thirds * m_one_third;
        let two = number * number;
        let three = two * number;
        let m_two = m_one * m_one;
        let m_three = m_two * m_one;
        let m_four = m_two * m_two;
        Self {
            one_sixth: sixth,
            m_one_sixth,
            itself: number,
            m_one,
            one_third,
            m_one_third,
            two_thirds,
            m_two_thirds,
            four_thirds,
            m_four_thirds,
            five_thirds,
            m_five_thirds,
            seven_thirds,
            m_seven_thirds,
            eight_thirds,
            m_eight_thirds,
            two,
            three,
            four: two * two,
            five: three * two,
            m_two,
            m_three,
            m_four,
            m_five: m_three * m_two,
            seven_sixths: sixth * number,
            m_seven_sixths: m_one_sixth * m_one,
            log: number.ln(),
            sqrt: sixth * one_third,
        }
    }
}

/// `powers_of_lalpi`: useful powers of pi, initialised as in LAL.
pub(crate) static POWERS_OF_LALPI: Lazy<UsefulPowers> = Lazy::new(|| UsefulPowers::new(PI));

/// Default dimensionless cutoff frequency, `fCutDef` (0.3, or 0.33 for
/// very large effective spin).
const F_CUT_DEF: f64 = 0.3;
const F_CUT_DEF_HIGH_SPIN: f64 = 0.33;

/// Equivalent of `IMRPhenomXWaveformStruct`: the frequency independent
/// quantities shared by all IMRPhenomX modes, built as in
/// `IMRPhenomXSetWaveformVariables` for the default model versions
/// (BBH, no tidal effects, no PNR deviations, GR).
///
/// Field names follow the LAL fields in snake case. Masses `m1`, `m2` are in
/// units of the total mass with `m1 >= m2` (the components are swapped if
/// necessary, exactly as LAL does); frequencies prefixed with `mf` or
/// without a unit (`f_ring`, `f_meco`, ...) are dimensionless (in units of
/// the total mass) while `f_min`, `f_max`, `f_ref`, `f_cut`, ... are in Hz.
#[derive(Clone, Copy, Debug)]
#[allow(dead_code)] // the full LAL struct is kept for the higher-mode models
pub(crate) struct WaveformStruct {
    // Mass parameters
    /// Larger mass in kg.
    pub(crate) m1_si: f64,
    /// Smaller mass in kg.
    pub(crate) m2_si: f64,
    /// Total mass in kg.
    pub(crate) m_tot_si: f64,
    /// Larger mass in units of the total mass.
    pub(crate) m1: f64,
    /// Smaller mass in units of the total mass.
    pub(crate) m2: f64,
    /// Total mass in solar masses.
    pub(crate) m_tot: f64,
    /// Mass ratio >= 1.
    pub(crate) q: f64,
    /// Symmetric mass ratio.
    pub(crate) eta: f64,
    /// PN symmetry parameter |m1 - m2| / M.
    pub(crate) delta: f64,
    pub(crate) eta2: f64,
    pub(crate) eta3: f64,

    // Spin parameters
    pub(crate) chi1l: f64,
    pub(crate) chi2l: f64,
    pub(crate) chi1l2l: f64,
    pub(crate) chi1l2: f64,
    pub(crate) chi1l3: f64,
    pub(crate) chi2l2: f64,
    pub(crate) chi2l3: f64,
    pub(crate) chi_eff: f64,
    pub(crate) chi_pn_hat: f64,
    pub(crate) s_tot_r: f64,
    pub(crate) dchi: f64,
    pub(crate) dchi_half: f64,
    /// (M/m2)(S2.L) - (M/m1)(S1.L) in units of M^2.
    pub(crate) sigma_l: f64,
    /// S1.L + S2.L in units of M^2.
    pub(crate) sl: f64,

    // Matter parameters (always the BBH values)
    pub(crate) lambda1: f64,
    pub(crate) lambda2: f64,
    pub(crate) quadparam1: f64,
    pub(crate) quadparam2: f64,

    // Frequencies and reference values
    /// Reference frequency in Hz (`f_min` if zero was requested).
    pub(crate) f_ref: f64,
    /// Orbital phase at the reference frequency as passed in.
    pub(crate) phi_ref_in: f64,
    /// Orbital phase at the reference frequency.
    pub(crate) phi0: f64,
    /// Azimuthal angle of the binary at the reference frequency.
    pub(crate) beta: f64,
    /// Phase at the reference frequency, set when generating the 22 mode.
    pub(crate) phif_ref: f64,
    /// Dimensionless reference frequency.
    pub(crate) mf_ref: f64,
    /// pi times the total mass in seconds.
    pub(crate) pi_m: f64,
    /// Orbital speed at the reference frequency.
    pub(crate) v_ref: f64,
    /// Total mass in seconds, converts Hz to dimensionless frequency.
    pub(crate) m_sec: f64,
    /// Frequency spacing in Hz, zero for non-uniform grids.
    pub(crate) delta_f: f64,
    pub(crate) delta_mf: f64,
    /// Dimensionless end frequency of the model.
    pub(crate) f_cut_def: f64,
    /// Minimum frequency in Hz.
    pub(crate) f_min: f64,
    /// Maximum frequency in Hz as requested (zero for the default).
    pub(crate) f_max: f64,
    pub(crate) mf_max: f64,
    /// Cutoff frequency in Hz, `f_cut_def / m_sec`.
    pub(crate) f_cut: f64,
    /// Effective maximum frequency, `min(f_max, f_cut)` (or `f_cut`).
    pub(crate) f_max_prime: f64,

    // Final state
    pub(crate) m_final: f64,
    pub(crate) a_final: f64,
    pub(crate) a_final_nonprec: f64,
    pub(crate) a_final_prec: f64,
    /// Ringdown frequency of the (2,2) mode.
    pub(crate) f_ring: f64,
    /// Damping frequency of the (2,2) mode.
    pub(crate) f_damp: f64,
    pub(crate) f_meco: f64,
    pub(crate) f_isco: f64,

    // Normalisation
    /// Luminosity distance in m.
    pub(crate) distance: f64,
    pub(crate) inclination: f64,
    /// Amplitude scale M^2 / D (geometric units, seconds).
    pub(crate) amp0: f64,
    /// Leading order amplitude prefactor sqrt(2 eta / 3) pi^(-1/6).
    pub(crate) amp_norm: f64,
    /// Phase normalisation 5 / (128 pi^(5/3)).
    pub(crate) dphase0: f64,
}

impl WaveformStruct {
    /// `IMRPhenomXSetWaveformVariables` for the default model versions.
    ///
    /// Masses are in solar masses (converted to kg and back as LAL receives
    /// them from `ChooseFDWaveform`), the distance in Mpc, frequencies in Hz
    /// and angles in radians. `delta_f = 0` denotes a non-uniform frequency
    /// grid, `f_ref = 0` selects `f_min` and `f_max = 0` the default cutoff.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        mass_1: f64,
        mass_2: f64,
        chi1l_in: f64,
        chi2l_in: f64,
        delta_f: f64,
        f_ref: f64,
        phi0: f64,
        f_min: f64,
        f_max: f64,
        distance: f64,
        inclination: f64,
    ) -> Result<Self, IMRPhenomXError> {
        if !(mass_1 > 0.0 && mass_2 > 0.0) {
            return Err(IMRPhenomXError::NonPositiveMass);
        }
        if distance < 0.0 {
            return Err(IMRPhenomXError::NegativeDistance);
        }
        if f_ref < 0.0 {
            return Err(IMRPhenomXError::NegativeReferenceFrequency);
        }

        // LAL receives the masses in kg and rescales to solar masses.
        let m1_in = (mass_1 * MSUN_SI) / MSUN_SI;
        let m2_in = (mass_2 * MSUN_SI) / MSUN_SI;

        let mass_ratio = if m1_in > m2_in {
            m1_in / m2_in
        } else {
            m2_in / m1_in
        };
        if mass_ratio > 1000. && (mass_ratio - 1000.).abs() > 1e-12 {
            return Err(IMRPhenomXError::MassRatioTooLarge(mass_ratio));
        }

        // Check if m1 >= m2, if not then swap masses and spins.
        let (m1, m2, chi1l, chi2l) = if m1_in >= m2_in {
            (m1_in, m2_in, chi1l_in, chi2l_in)
        } else {
            (m2_in, m1_in, chi2l_in, chi1l_in)
        };

        // LAL attempts to nudge spins within 1e-6 of +-1 onto the Kerr bound
        // with IMRPhenomX_InternalNudge, but that function takes its argument
        // by value and so has no effect: any spin outside [-1, 1] is an error.
        for chi in [chi1l, chi2l] {
            if !(-1.0..=1.0).contains(&chi) {
                return Err(IMRPhenomXError::UnphysicalSpin(chi));
            }
        }

        // Symmetric mass ratio
        let delta = ((m1 - m2) / (m1 + m2)).abs();
        // use fabs to prevent negative sign due to roundoff; 1 - delta^2 is
        // contracted to a fused multiply-add like LAL's arm64 build.
        let mut eta = (0.25 * crate::spintaylor::fp::fnms(delta, delta, 1.0)).abs();
        let mut q = if m1 > m2 { m1 / m2 } else { m2 / m1 };
        if eta > 0.25 {
            eta = 0.25;
        }
        if eta == 0.25 {
            q = 1.;
        }

        let m1_si = m1 * MSUN_SI;
        let m2_si = m2 * MSUN_SI;
        let m_tot = m1 + m2;
        let m_sec = m_tot * MTSUN_SI;

        let chi_eff = chi_eff(eta, chi1l, chi2l);
        let dchi = dchi(chi1l, chi2l);
        let m1 = m1 / m_tot;
        let m2 = m2 / m_tot;

        // If no reference frequency is passed, set it to the starting frequency.
        let f_ref = if f_ref == 0.0 { f_min } else { f_ref };
        let pi_m = PI * m_sec;

        // Default end of the waveform: 0.3 Mf, raised to 0.33 if chi_eff is
        // very high so that the ringdown of the 44 mode is not cut out.
        let f_cut_def = if chi_eff > 0.99 {
            F_CUT_DEF_HIGH_SPIN
        } else {
            F_CUT_DEF
        };
        let f_cut = f_cut_def / m_sec;
        // By default f_max_prime is f_max. If f_cut < f_max, use f_cut.
        let mut f_max_prime = if f_max != 0.0 { f_max } else { f_cut };
        if f_max_prime > f_cut {
            f_max_prime = f_cut;
        }

        let distance = distance * MPC_SI;

        let mut wf = Self {
            m1_si,
            m2_si,
            m_tot_si: m1_si + m2_si,
            m1,
            m2,
            m_tot,
            q,
            eta,
            delta,
            eta2: eta * eta,
            eta3: eta * (eta * eta),
            chi1l,
            chi2l,
            chi1l2l: chi1l * chi2l,
            chi1l2: chi1l * chi1l,
            chi1l3: chi1l * chi1l * chi1l,
            chi2l2: chi2l * chi2l,
            chi2l3: chi2l * chi2l * chi2l,
            chi_eff,
            chi_pn_hat: chi_pn_hat(eta, chi1l, chi2l),
            s_tot_r: s_tot_r(eta, chi1l, chi2l),
            dchi,
            dchi_half: dchi * 0.5,
            sigma_l: (chi2l * m2) - (chi1l * m1),
            sl: chi1l * (m1 * m1) + chi2l * (m2 * m2),
            lambda1: 0.0,
            lambda2: 0.0,
            quadparam1: 1.0,
            quadparam2: 1.0,
            f_ref,
            phi_ref_in: phi0,
            phi0,
            beta: PI * 0.5 - phi0,
            phif_ref: 0.0,
            mf_ref: hz_to_mf(f_ref, m_tot),
            pi_m,
            v_ref: (pi_m * f_ref).cbrt(),
            m_sec,
            delta_f,
            delta_mf: hz_to_mf(delta_f, m_tot),
            f_cut_def,
            f_min,
            f_max,
            mf_max: hz_to_mf(f_max, m_tot),
            f_cut,
            f_max_prime,
            m_final: 0.0,
            a_final: 0.0,
            a_final_nonprec: 0.0,
            a_final_prec: 0.0,
            f_ring: 0.0,
            f_damp: 0.0,
            f_meco: 0.0,
            f_isco: 0.0,
            distance,
            inclination,
            amp0: m_tot * MRSUN_SI * m_tot * MTSUN_SI / distance,
            amp_norm: (2.0f64 / 3.0).sqrt() * eta.sqrt() * POWERS_OF_LALPI.m_one_sixth,
            dphase0: 5.0 / (128.0 * PI.powf(5.0 / 3.0)),
        };
        wf.set_final_state()?;
        Ok(wf)
    }

    /// Default final mass and spin (`XLALSimIMRPhenomXFinalMass2017`,
    /// `XLALSimIMRPhenomXFinalSpin2017`), the corresponding (2,2) QNM
    /// frequencies and the MECO/ISCO frequencies.
    fn set_final_state(&mut self) -> Result<(), IMRPhenomXError> {
        self.m_final = final_mass_2017(self.eta, self.chi1l, self.chi2l);
        self.a_final = final_spin_2017(self.eta, self.chi1l, self.chi2l);
        self.a_final_nonprec = self.a_final;
        self.a_final_prec = self.a_final;
        self.set_ringdown_frequencies();

        if self.m_final > 1.0 {
            return Err(IMRPhenomXError::UnphysicalFinalMass(self.m_final));
        }
        if self.a_final.abs() > 1.0 {
            return Err(IMRPhenomXError::UnphysicalFinalSpin(self.a_final));
        }

        self.f_meco = f_meco(self.eta, self.chi1l, self.chi2l);
        self.f_isco = f_isco(self.a_final);
        if self.f_meco > self.f_isco {
            return Err(IMRPhenomXError::MecoAboveIsco {
                f_meco: self.f_meco,
                f_isco: self.f_isco,
            });
        }
        Ok(())
    }

    /// Recompute `f_ring` and `f_damp` from the current `a_final` and
    /// `m_final` (`evaluate_QNMfit_fring22 / Mfinal`).
    pub(crate) fn set_ringdown_frequencies(&mut self) {
        self.f_ring = fring22(self.a_final) / self.m_final;
        self.f_damp = fdamp22(self.a_final) / self.m_final;
    }

    /// Override the final spin and update the (2,2) QNM frequencies, as the
    /// precessing models do. Like LAL, `f_isco` keeps its aligned-spin value.
    #[allow(dead_code)]
    pub(crate) fn set_final_spin(&mut self, a_final: f64) {
        self.a_final = a_final;
        self.set_ringdown_frequencies();
    }

    /// Set the reference frequency in Hz (zero selects `f_min`).
    pub(crate) fn set_reference_frequency(&mut self, f_ref: f64) {
        self.f_ref = if f_ref == 0.0 { self.f_min } else { f_ref };
        self.mf_ref = hz_to_mf(self.f_ref, self.m_tot);
        self.v_ref = (self.pi_m * self.f_ref).cbrt();
    }

    /// Effective maximum frequency for a requested maximum frequency
    /// `f_max` in Hz: `f_max_prime` of `IMRPhenomXSetWaveformVariables`.
    pub(crate) fn f_max_prime(&self, f_max: f64) -> f64 {
        let f_max_prime = if f_max != 0.0 { f_max } else { self.f_cut };
        if f_max_prime > self.f_cut {
            self.f_cut
        } else {
            f_max_prime
        }
    }
}

/// `XLALSimIMRPhenomXUtilsHztoMf`: frequency in Hz to dimensionless.
pub(crate) fn hz_to_mf(f_hz: f64, m_tot_msun: f64) -> f64 {
    f_hz * (MTSUN_SI * m_tot_msun)
}

#[cfg(test)]
mod tests {
    use super::UsefulPowers;

    #[test]
    fn fast_powers_agree_with_lal_powers() {
        for &x in &[1e-4, 3.7e-3, 0.018, 0.0712, 0.3] {
            let a = UsefulPowers::new(x);
            let b = UsefulPowers::fast(x);
            let pairs = [
                (a.one_sixth, b.one_sixth),
                (a.m_one_sixth, b.m_one_sixth),
                (a.one_third, b.one_third),
                (a.m_one_third, b.m_one_third),
                (a.two_thirds, b.two_thirds),
                (a.m_two_thirds, b.m_two_thirds),
                (a.four_thirds, b.four_thirds),
                (a.m_four_thirds, b.m_four_thirds),
                (a.five_thirds, b.five_thirds),
                (a.m_five_thirds, b.m_five_thirds),
                (a.seven_thirds, b.seven_thirds),
                (a.m_seven_thirds, b.m_seven_thirds),
                (a.eight_thirds, b.eight_thirds),
                (a.m_eight_thirds, b.m_eight_thirds),
                (a.two, b.two),
                (a.three, b.three),
                (a.four, b.four),
                (a.five, b.five),
                (a.m_one, b.m_one),
                (a.m_two, b.m_two),
                (a.m_three, b.m_three),
                (a.m_four, b.m_four),
                (a.m_five, b.m_five),
                (a.seven_sixths, b.seven_sixths),
                (a.m_seven_sixths, b.m_seven_sixths),
                (a.log, b.log),
                (a.sqrt, b.sqrt),
            ];
            for (i, (a, b)) in pairs.iter().enumerate() {
                assert!(((a - b) / a).abs() < 2e-15, "{i}: {a} != {b} at {x}");
            }
        }
    }
}
