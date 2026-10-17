//! Mode-dependent variables and coefficients of IMRPhenomXHM,
//! `LALSimIMRPhenomXHM_internals.c`, for the default model versions
//! (`PhenomXHMReleaseVersion` 122022 with the 122019 phase).
//!
//! Every coefficient is computed operation by operation as in LALSimulation,
//! including its fused multiply-adds (see [`fma`]), so that they agree with
//! LAL to the last bit. This matters for the (3,2) mode, whose intermediate
//! phase is fixed by a second derivative computed by finite differences with
//! a step of 1e-7 (amplifying rounding differences by ~1e14), and for the
//! numerical derivatives of the amplitude at the region boundaries.

use std::f64::consts::{FRAC_PI_4, PI};

use num_complex::Complex;

use crate::imrphenomx::amplitude::AmpCoefficients;
use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct};
use crate::imrphenomx::linalg::lu_solve;
use crate::imrphenomx::phase::{
    dphase_22, phase_22, phase_22_connection_coefficients, time_shift_22, PhaseCoefficients,
};
use crate::imrphenomx::utilities::psi4_to_strain;
use crate::imrphenomx::{hm_inspiral, hm_intermediate, hm_qnm, hm_ringdown, libm, ringdown};

/// `FALSE_ZERO`: value replacing negative amplitudes.
const FALSE_ZERO: f64 = 1.0e-15;

/// Finite-difference step of the (3,2) ringdown phase derivatives in
/// `IMRPhenomXHM_GetPhaseCoefficients`.
const FSTEP_MIXING: f64 = 0.0000001;

/// A parameter-space fit of IMRPhenomXHM, `ParameterSpaceFit`.
pub(crate) type Fit = fn(&WaveformStruct) -> f64;

/// A quasi-normal mode fit as a function of the final spin.
type QnmFit = fn(f64) -> f64;

/// The higher modes implemented in IMRPhenomXHM, in the order of LAL's
/// `modeInt`: (2,1), (3,3), (3,2), (4,4).
pub(crate) const HIGHER_MODES: [(u32, u32); 4] = [(2, 1), (3, 3), (3, 2), (4, 4)];

/// Equivalent of `IMRPhenomXHMWaveformStruct`: the variables of one higher
/// mode, set as in `IMRPhenomXHM_SetHMWaveformVariables` for the default
/// release (122022). Frequencies are in units of the total mass.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HMWaveformStruct {
    pub(crate) ell: u32,
    pub(crate) emm: u32,
    /// `10 ell + emm`.
    pub(crate) mode_tag: u32,
    /// Index of the mode in [`HIGHER_MODES`].
    pub(crate) mode_int: usize,
    /// Whether the mode is built with spheroidal-spherical mode mixing (3,2).
    pub(crate) mixing_on: bool,
    /// Odd modes vanish for equal masses and equal spins.
    pub(crate) amp_zero: bool,
    /// Collocation points used for the intermediate amplitude, one digit per
    /// point (`IMRPhenomXHMIntermediateAmpVersion`): 0 skip, 1 value,
    /// 2 value and derivative.
    pub(crate) intermediate_amp_version: [u8; 6],
    /// `IMRPhenomXHMRingdownAmpVersion`: 2 (collocation points) or 1 (fits of
    /// the coefficients, used by the (3,2) mode).
    pub(crate) ringdown_amp_version: u32,
    /// Number of intermediate phase collocation points (6 for the (3,2)).
    pub(crate) n_colloc_pts_inter_phase: usize,
    /// Limit between comparable and extreme mass ratios for the phase.
    pub(crate) eta_emr: f64,
    pub(crate) chi_s: f64,
    pub(crate) chi_a: f64,
    /// MECO frequency of the mode, `m fMECO / 2`.
    pub(crate) f_meco_lm: f64,
    /// Ringdown frequency of the mode.
    pub(crate) f_ring: f64,
    /// Damping frequency of the mode.
    pub(crate) f_damp: f64,
    /// Linear-in-frequency term of the 22 phase (`IMRPhenomX_TimeShift_22`),
    /// set with the phase coefficients.
    pub(crate) timeshift: f64,
    /// Constant term of the 22 phase (zero).
    pub(crate) phaseshift: f64,
    /// Phase offset `phifRef` of the 22, set with the phase coefficients.
    pub(crate) phiref22: f64,
    /// Mixing coefficients (rounded to single precision as in LAL's
    /// `COMPLEX8`): a222, a223, -a322, -a323.
    pub(crate) mixing_coeffs: [Complex<f64>; 4],
    /// Leading-order normalisation of the 22, `ampNorm`.
    pub(crate) amp_norm: f64,
    /// Amplitude scale M^2 / D, `amp0` (without `ampNorm`).
    pub(crate) amp0: f64,
    /// Time shift of the hybrids, `-2 pi (500 + psi4ToStrain)`.
    pub(crate) delta_t: f64,
    /// Start of the exponential fall-off of the ringdown amplitude.
    pub(crate) f_amp_rd_falloff: f64,
}

impl HMWaveformStruct {
    /// `IMRPhenomXHM_SetHMWaveformVariables` for mode (`ell`, `emm`), `emm > 0`.
    pub(crate) fn new(ell: u32, emm: u32, wf: &WaveformStruct) -> Self {
        let mode_tag = 10 * ell + emm;
        let mode_int = HIGHER_MODES
            .iter()
            .position(|&lm| lm == (ell, emm))
            .unwrap_or_else(|| panic!("IMRPhenomXHM: mode ({ell},{emm}) is not available"));
        let odd = mode_tag == 21 || mode_tag == 33;
        let amp_zero = odd && wf.q == 1. && wf.chi1l == wf.chi2l;
        let mixing_on = mode_tag == 32;
        let intermediate_amp_version = if mode_tag == 21 {
            [1, 1, 0, 1, 0, 2]
        } else {
            [2, 1, 1, 1, 1, 2]
        };
        let ringdown_amp_version = if mode_tag == 32 { 1 } else { 2 };
        let (fring, fdamp): (QnmFit, QnmFit) = match mode_tag {
            21 => (hm_qnm::fring21, hm_qnm::fdamp21),
            33 => (hm_qnm::fring33, hm_qnm::fdamp33),
            32 => (hm_qnm::fring32, hm_qnm::fdamp32),
            _ => (hm_qnm::fring44, hm_qnm::fdamp44),
        };
        let mixing_coeffs = if mixing_on {
            mixing_coefficients(wf.a_final)
        } else {
            [Complex::new(0.0, 0.0); 4]
        };
        let psi4tostrain = psi4_to_strain(wf.eta, wf.s_tot_r, wf.dchi);
        Self {
            ell,
            emm,
            mode_tag,
            mode_int,
            mixing_on,
            amp_zero,
            intermediate_amp_version,
            ringdown_amp_version,
            n_colloc_pts_inter_phase: if mixing_on { 6 } else { 5 },
            eta_emr: 0.05,
            chi_s: (wf.chi1l + wf.chi2l) * 0.5,
            chi_a: (wf.chi1l - wf.chi2l) * 0.5,
            f_meco_lm: wf.f_meco * emm as f64 * 0.5,
            f_ring: fring(wf.a_final) / wf.m_final,
            f_damp: fdamp(wf.a_final) / wf.m_final,
            timeshift: 0.0,
            phaseshift: 0.0,
            phiref22: 0.0,
            mixing_coeffs,
            amp_norm: wf.amp_norm,
            amp0: wf.amp0,
            delta_t: -2. * PI * (500. + psi4tostrain),
            f_amp_rd_falloff: 0.0,
        }
    }
}

/// `IMRPhenomXHM_Initialize_MixingCoeffs`: the mixing coefficients stored in
/// single precision, with the signs of the (3,2) ones flipped to match the
/// conventions of the hybrids.
fn mixing_coefficients(af: f64) -> [Complex<f64>; 4] {
    let single = |re: f64, im: f64| Complex::new(re as f32 as f64, im as f32 as f64);
    [
        single(hm_qnm::re_l2m2lp2(af), hm_qnm::im_l2m2lp2(af)),
        single(hm_qnm::re_l2m2lp3(af), hm_qnm::im_l2m2lp3(af)),
        -single(hm_qnm::re_l3m2lp2(af), hm_qnm::im_l3m2lp2(af)),
        -single(hm_qnm::re_l3m2lp3(af), hm_qnm::im_l3m2lp3(af)),
    ]
}

/// Equivalent of `IMRPhenomXHMAmpCoefficients` for the 122022 release.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct HMAmpCoefficients {
    /// Inspiral-intermediate transition, `fAmpMatchIN`.
    pub(crate) f_amp_match_in: f64,
    /// Intermediate-ringdown transition, `fAmpMatchIM`.
    pub(crate) f_amp_match_im: f64,
    /// End of the polynomial ringdown region of the (3,2) mode, `fRDAux`.
    pub(crate) f_rd_aux: f64,
    /// PN amplitude coefficients `pnInitial` to `pnSixThirds`.
    pub(crate) pn: [Complex<f64>; 7],
    /// `ampNorm (2/m)^(-7/6)`.
    pub(crate) pn_dominant: f64,
    /// Leading PN amplitude of the mode relative to the 22.
    pub(crate) pn_global_factor: f64,
    pub(crate) amp_norm: f64,
    pub(crate) colloc_freqs_insp: [f64; 3],
    pub(crate) colloc_values_insp: [f64; 3],
    pub(crate) pn_amplitude_insp: [f64; 3],
    /// Pseudo-PN coefficients of f^(7/3), f^(8/3) and f^3.
    pub(crate) inspiral_coefficients: [f64; 3],
    pub(crate) fcut_insp_seven_thirds: f64,
    pub(crate) fcut_insp_eight_thirds: f64,
    pub(crate) fcut_insp_three: f64,
    pub(crate) n_coefficients_inter: usize,
    pub(crate) colloc_freqs_inter: [f64; 8],
    pub(crate) colloc_values_inter: [f64; 8],
    /// Coefficients of the intermediate polynomial (times f^(-7/6)).
    pub(crate) inter_coefficients: [f64; 8],
    pub(crate) colloc_freqs_rd: [f64; 3],
    pub(crate) colloc_values_rd: [f64; 3],
    /// Ringdown coefficients: a_lambda, lambda, sigma, and the value and
    /// decay rate of the exponential fall-off.
    pub(crate) rd_coefficients: [f64; 5],
    pub(crate) n_coefficients_rd_aux: usize,
    pub(crate) colloc_freqs_rd_aux: [f64; 4],
    pub(crate) colloc_values_rd_aux: [f64; 4],
    /// Cubic polynomial of the auxiliary ringdown region of the (3,2).
    pub(crate) rd_aux_coefficients: [f64; 4],
}

/// Equivalent of `IMRPhenomXHMPhaseCoefficients` for the 122019 phase.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct HMPhaseCoefficients {
    pub(crate) f_phase_match_in: f64,
    pub(crate) f_phase_match_im: f64,
    /// Constant phase aligning the mode with the 22 in the inspiral.
    pub(crate) deltaphi_lm: f64,
    /// Ringdown coefficients without mode mixing.
    pub(crate) alpha0: f64,
    pub(crate) alpha2: f64,
    pub(crate) alpha_l: f64,
    /// Ringdown phase and derivative at the intermediate-ringdown boundary.
    pub(crate) phi0_rd: f64,
    pub(crate) dphi0_rd: f64,
    /// Spheroidal ringdown coefficients of the (3,2) mode.
    pub(crate) phi0_s: f64,
    pub(crate) alpha0_s: f64,
    pub(crate) alpha2_s: f64,
    pub(crate) alpha4_s: f64,
    pub(crate) alpha_l_s: f64,
    /// Intermediate coefficients.
    pub(crate) c0: f64,
    pub(crate) c1: f64,
    pub(crate) c2: f64,
    pub(crate) c3: f64,
    pub(crate) c4: f64,
    pub(crate) c_l: f64,
    /// Connection coefficients.
    pub(crate) c_rd: f64,
    pub(crate) c1_rd: f64,
    pub(crate) c_insp: f64,
    pub(crate) c1_insp: f64,
    /// Rescaled inspiral coefficients of f^((i-5)/3) and their logarithms.
    pub(crate) phi: [f64; 13],
    pub(crate) phi_l: [f64; 13],
    pub(crate) lambda_pn: f64,
    pub(crate) colloc_freqs_rd: [f64; 4],
    pub(crate) colloc_values_rd: [f64; 4],
    pub(crate) colloc_freqs_inter: [f64; 6],
    pub(crate) colloc_values_inter: [f64; 6],
}

/// The 22-mode data needed by the spheroidal-to-spherical rotation of the
/// (3,2) ringdown when the 22 is not recycled.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Mixing22 {
    /// Amplitude coefficients of the 22 (`pAmp22`).
    pub(crate) amp22: AmpCoefficients,
    /// Phase coefficients of the 22 including the connection coefficients.
    pub(crate) phase22: PhaseCoefficients,
    pub(crate) f_ring22: f64,
    pub(crate) f_damp22: f64,
    /// 1 / eta
    pub(crate) inveta: f64,
}

/// All frequency-independent data of one higher mode.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HMMode {
    pub(crate) wfhm: HMWaveformStruct,
    pub(crate) amp: HMAmpCoefficients,
    pub(crate) phase: HMPhaseCoefficients,
    /// Present for the (3,2) mode.
    pub(crate) mixing: Option<Mixing22>,
}

impl HMMode {
    /// Build mode (`ell`, `emm`) as in the mode loop of
    /// `IMRPhenomXHM_MultiMode2` / `IMRPhenomXHMGenerateFDOneMode`, from the
    /// waveform variables and the 22 phase coefficients returned by
    /// `IMRPhenomXGetPhaseCoefficients` (without connection coefficients).
    pub(crate) fn new(
        ell: u32,
        emm: u32,
        wf: &WaveformStruct,
        phase22: &PhaseCoefficients,
    ) -> Self {
        let mut wfhm = HMWaveformStruct::new(ell, emm, wf);
        let mut amp = HMAmpCoefficients::default();
        let mut phase = HMPhaseCoefficients::default();
        if wfhm.amp_zero {
            return Self {
                wfhm,
                amp,
                phase,
                mixing: None,
            };
        }
        let mut phase22 = *phase22;
        let mixing = if wfhm.mixing_on {
            get_spheroidal_coefficients(&mut phase, &mut phase22, &mut wfhm, wf);
            Some(Mixing22 {
                amp22: AmpCoefficients::new(wf),
                phase22,
                f_ring22: wf.f_ring,
                f_damp22: wf.f_damp,
                inveta: 1.0 / wf.eta,
            })
        } else {
            None
        };
        get_amplitude_coefficients(&mut amp, &phase, mixing.as_ref(), &mut wfhm, wf);
        get_phase_coefficients(
            &amp,
            &mut phase,
            mixing.as_ref(),
            &mut phase22,
            &mut wfhm,
            wf,
        );
        Self {
            wfhm,
            amp,
            phase,
            mixing,
        }
    }

    /// `IMRPhenomXHM_Amplitude_noModeMixing` / `_ModeMixing` /
    /// `_ModeMixingRecycle`: the amplitude (strain, without `amp0`) at the
    /// dimensionless frequency of `powers`, with `rotation` the spherical
    /// (3,2) ringdown.
    #[inline]
    fn amplitude_with<F: Fn(&UsefulPowers) -> Complex<f64>>(
        &self,
        powers: &UsefulPowers,
        rotation: F,
    ) -> f64 {
        if self.wfhm.amp_zero {
            return 0.0;
        }
        let mf = powers.itself;
        let amp = if mf < self.amp.f_amp_match_in {
            hm_inspiral::amp_ansatz(powers, &self.wfhm, &self.amp)
        } else if mf >= self.amp.f_amp_match_im {
            if self.wfhm.mixing_on {
                let h = rotation(powers);
                h.re.hypot(h.im)
            } else {
                hm_ringdown::amp_ansatz(powers, &self.wfhm, &self.amp)
            }
        } else {
            hm_intermediate::amp_ansatz(powers, &self.amp)
        };
        if amp < 0.0 {
            FALSE_ZERO
        } else {
            amp
        }
    }

    /// `IMRPhenomXHM_Phase_noModeMixing` / `_ModeMixing` /
    /// `_ModeMixingRecycle`.
    #[inline]
    fn phase_with<F: Fn(&UsefulPowers) -> Complex<f64>>(
        &self,
        powers: &UsefulPowers,
        rotation: F,
    ) -> f64 {
        let mf = powers.itself;
        let p = &self.phase;
        if mf < p.f_phase_match_in {
            let phi_ins = hm_inspiral::phase_ansatz_int(mf, powers, p);
            fma(p.c1_insp, mf, phi_ins) + p.c_insp + p.deltaphi_lm
        } else if mf >= p.f_phase_match_im {
            let phi_mrd = if self.wfhm.mixing_on {
                let h = rotation(powers);
                h.im.atan2(h.re)
            } else {
                hm_ringdown::phase_ansatz_int(mf, powers, &self.wfhm, p)
            };
            fma(p.c1_rd, mf, phi_mrd) + p.c_rd + p.deltaphi_lm
        } else {
            hm_intermediate::phase_ansatz_int(mf, powers, &self.wfhm, p) + p.deltaphi_lm
        }
    }

    /// Amplitude and phase of the mode at the dimensionless frequency of
    /// `powers`, evaluating the spheroidal-to-spherical rotation of the (3,2)
    /// ringdown at most once. `wf22` is the (2,2) mode (strain, including
    /// `amp0`) to recycle in the (3,2) ringdown as `IMRPhenomXHM_MultiMode2`
    /// does; without it the ringdown 22 is computed internally as in
    /// `IMRPhenomXHMGenerateFDOneMode`.
    #[inline]
    pub(crate) fn amplitude_phase(
        &self,
        powers: &UsefulPowers,
        wf22: Option<Complex<f64>>,
    ) -> (f64, f64) {
        let rotation = std::cell::OnceCell::new();
        let rotate = |p: &UsefulPowers| {
            *rotation.get_or_init(|| match wf22 {
                Some(wf22) => self.spheroidal_to_spherical_recycle(p, wf22),
                None => self.spheroidal_to_spherical(p),
            })
        };
        (
            self.amplitude_with(powers, rotate),
            self.phase_with(powers, rotate),
        )
    }

    /// `SpheroidalToSpherical`: the spherical (3,2) mode in the ringdown from
    /// the spheroidal (3,2) ansatz and the ringdown 22 mode, both rescaled by
    /// the leading-order 22 amplitude.
    pub(crate) fn spheroidal_to_spherical(&self, powers: &UsefulPowers) -> Complex<f64> {
        let mix = self
            .mixing
            .as_ref()
            .expect("IMRPhenomXHM: spheroidal rotation without mixing data");
        let mf = powers.itself;
        let amp22 = ringdown::amp_22_ansatz_exact(mf, mix.f_ring22, mix.f_damp22, &mix.amp22);
        let phi22 = fma(
            mix.inveta,
            phase_22(mf, powers, &mix.phase22),
            self.wfhm.timeshift * mf,
        ) + self.wfhm.phaseshift
            + self.wfhm.phiref22;
        let rotation = libm::cexp_i(phi22);
        let mut wf22r = Complex::new(amp22 * rotation.re, amp22 * rotation.im);
        // RingdownAmpVersion != 0: rescale with the leading order.
        let rescale = self.wfhm.amp_norm * powers.m_seven_sixths;
        wf22r = Complex::new(wf22r.re * rescale, wf22r.im * rescale);
        self.rotate(powers, wf22r)
    }

    /// `SpheroidalToSphericalRecycle`: as
    /// [`HMMode::spheroidal_to_spherical`] with the 22 mode `wf22` (strain,
    /// including `amp0`) computed elsewhere.
    #[inline]
    pub(crate) fn spheroidal_to_spherical_recycle(
        &self,
        powers: &UsefulPowers,
        wf22: Complex<f64>,
    ) -> Complex<f64> {
        let wf22r = Complex::new(wf22.re / self.wfhm.amp0, wf22.im / self.wfhm.amp0);
        self.rotate(powers, wf22r)
    }

    /// `conj(a322) h22 + conj(a323) h32_spheroidal`, with the complex
    /// products expanded as clang does (no contraction).
    #[inline]
    fn rotate(&self, powers: &UsefulPowers, wf22r: Complex<f64>) -> Complex<f64> {
        let mf = powers.itself;
        let amplm = hm_ringdown::amp_ansatz(powers, &self.wfhm, &self.amp);
        let philm = hm_ringdown::phase_ansatz_int(mf, powers, &self.wfhm, &self.phase);
        let [_, _, mix2, mix3] = self.wfhm.mixing_coeffs;
        // conj(mix2) * wf22r
        let a_re = wf22r.re * mix2.re + wf22r.im * mix2.im;
        let a_im = wf22r.im * mix2.re - wf22r.re * mix2.im;
        // conj(mix3) * amplm * exp(i philm)
        let b_re = amplm * mix3.re;
        let b_im = -(mix3.im * amplm);
        let rotation = libm::cexp_i(philm);
        let c_re = b_re * rotation.re - b_im * rotation.im;
        let c_im = b_re * rotation.im + b_im * rotation.re;
        Complex::new(a_re + c_re, a_im + c_im)
    }
}

/* ---------------------------------------------------------------------- */
/*                        Cutting frequencies                              */
/* ---------------------------------------------------------------------- */

/// `IMRPhenomXHM_Amplitude_fcutInsp` (version 122022).
fn amplitude_fcut_insp(wfhm: &HMWaveformStruct, wf: &WaveformStruct) -> f64 {
    let emm = wfhm.emm as f64;
    let eta = wf.eta;
    let chi1 = wf.chi1l;
    let f_meco = wfhm.f_meco_lm;
    let fcut_emr = 1.25
        * emm
        * (fma(
            0.000316707064291237,
            chi1 * chi1,
            fma(-0.0000858396080377194, chi1, 0.011671068725758493),
        ) * fma(6.2873167352395125, eta, 0.8447212540381764))
        / fma(-0.9977728883419751, chi1, 1.2857082764038923);
    if wf.q < 20. {
        f_meco
    } else {
        let transition_eta = 0.0192234; // q = 50
        let sharpness = 0.004;
        let funcs = fma(0.5, ((eta - transition_eta) / sharpness).tanh(), 0.5);
        fma(funcs, f_meco, (1. - funcs) * fcut_emr)
    }
}

/// `IMRPhenomXHM_Amplitude_fcutRD` (version 122022).
fn amplitude_fcut_rd(wfhm: &HMWaveformStruct, wf: &WaveformStruct) -> f64 {
    if wfhm.mixing_on {
        fma(-0.5, wf.f_damp, wf.f_ring)
    } else {
        wfhm.f_ring - wfhm.f_damp
    }
}

/// `GetfcutInsp`: end of the inspiral phase region.
fn get_fcut_insp(wf: &WaveformStruct, wfhm: &HMWaveformStruct) -> f64 {
    let eta_factor = fma(0.001, 0.25 / wf.eta - 1., 1.0);
    eta_factor * wfhm.emm as f64 * 0.5 * wf.f_meco
}

/// `IMRPhenomXHM_Intermediate_CollocPtsFreqs` (version 122019).
fn intermediate_colloc_pts_freqs(
    phase: &mut HMPhaseCoefficients,
    wfhm: &HMWaveformStruct,
    wf: &WaveformStruct,
) {
    let fring = wfhm.f_ring;
    let fdamp = wfhm.f_damp;
    let fcut = get_fcut_insp(wf, wfhm);
    let sqrt3 = 3f64.sqrt();
    let freqs = &mut phase.colloc_freqs_inter;
    freqs[0] = fcut;
    if wfhm.mode_tag == 32 {
        let f_end = fma(-0.5, wf.f_damp, wf.f_ring);
        freqs[1] = fma(sqrt3, fcut - f_end, 2. * (fcut + f_end)) / 4.;
        freqs[2] = fma(3., fcut, f_end) / 4.;
        freqs[3] = (fcut + f_end) / 2.;
        // The first and second derivatives at f_end are used.
        freqs[4] = f_end;
        freqs[5] = f_end;
        phase.f_phase_match_im = f_end;
        // Correct the cutting frequency for EMR with negative spins.
        if wf.eta < 0.01 && wf.chi1l < 0. {
            phase.f_phase_match_im *= fma(-0.25, wf.chi1l, 1.2);
        }
    } else {
        freqs[1] = fma(sqrt3, fcut - fring, 2. * (fcut + fring)) / 4.;
        freqs[2] = fma(3., fcut, fring) / 4.;
        freqs[3] = (fcut + fring) / 2.;
        freqs[4] = fma(3., fring, fcut) / 4.;
        freqs[5] = fma(7., fring, fcut) / 8.;
        phase.f_phase_match_im = fring - fdamp;
    }
    phase.f_phase_match_in = wfhm.f_meco_lm;
}

/* ---------------------------------------------------------------------- */
/*                          Amplitude coefficients                         */
/* ---------------------------------------------------------------------- */

/// `IMRPhenomXHM_GetAmplitudeCoefficients` (release 122022).
fn get_amplitude_coefficients(
    amp: &mut HMAmpCoefficients,
    phase: &HMPhaseCoefficients,
    mixing: Option<&Mixing22>,
    wfhm: &mut HMWaveformStruct,
    wf: &WaveformStruct,
) {
    amp.amp_norm = wf.amp_norm;
    amp.pn_dominant = wf.amp_norm * (2. / wfhm.emm as f64).powf(-7. / 6.);
    wfhm.f_amp_rd_falloff = 0.;
    amp.n_coefficients_inter = wfhm
        .intermediate_amp_version
        .iter()
        .map(|&d| d as usize)
        .sum();
    amp.n_coefficients_rd_aux = 0;
    if wfhm.mixing_on {
        amp.n_coefficients_rd_aux = 4;
        amp.f_rd_aux = wfhm.f_ring - wfhm.f_damp;
    }
    amp.f_amp_match_in = amplitude_fcut_insp(wfhm, wf);
    amp.f_amp_match_im = amplitude_fcut_rd(wfhm, wf);
    wfhm.f_amp_rd_falloff = fma(2., wfhm.f_damp, wfhm.f_ring);

    hm_inspiral::get_pn_amplitude_coefficients(amp, wfhm, wf);
    hm_inspiral::get_inspiral_amp_coefficients(amp, wfhm, wf);
    hm_ringdown::rd_amp_coefficients(amp, wfhm, wf);
    // The value of the (3,2) at the right boundary comes from the rotated
    // spheroidal ansatz with the coefficients computed so far.
    let wfhm: &HMWaveformStruct = wfhm;
    hm_intermediate::intermediate_amp_coefficients(amp, wfhm, wf, |powers, amp| {
        rotated_32(powers, *amp, phase, mixing, wfhm)
    });
}

/// `SpheroidalToSpherical` with the coefficients available while the
/// amplitude coefficients are being computed.
fn rotated_32(
    powers: &UsefulPowers,
    amp: HMAmpCoefficients,
    phase: &HMPhaseCoefficients,
    mixing: Option<&Mixing22>,
    wfhm: &HMWaveformStruct,
) -> Complex<f64> {
    let mode = HMMode {
        wfhm: *wfhm,
        amp,
        phase: *phase,
        mixing: mixing.copied(),
    };
    mode.spheroidal_to_spherical(powers)
}

/* ---------------------------------------------------------------------- */
/*                            Phase coefficients                           */
/* ---------------------------------------------------------------------- */

/// `IMRPhenomXHM_PN21AmpSign`: sign of the real part of the 21 amplitude
/// expanded in frequency.
fn pn21_amp_sign(ff: f64, wf: &WaveformStruct) -> i32 {
    let eta = wf.eta;
    let chi1 = wf.chi1l;
    let chi2 = wf.chi2l;
    let delta = fma(-4., eta, 1.).sqrt();
    let output = (-16. * delta * eta * ff * PI.powf(1.5)) / (3. * 5f64.sqrt())
        + (4.
            * 2f64.powf(0.3333333333333333)
            * (chi1 - chi2 + delta * (chi1 + chi2))
            * eta
            * ff.powf(1.3333333333333333)
            * PI.powf(1.8333333333333333))
            / 5f64.sqrt()
        + (2.
            * 2f64.powf(0.6666666666666666)
            * eta
            * (306. * delta - 360. * delta * eta)
            * ff.powf(1.6666666666666667)
            * PI.powf(2.1666666666666665))
            / (189. * 5f64.sqrt());
    if output >= 0. {
        1
    } else {
        -1
    }
}

/// `IMRPhenomXHM_GetPhaseCoefficients` (version 122019).
fn get_phase_coefficients(
    amp: &HMAmpCoefficients,
    phase: &mut HMPhaseCoefficients,
    mixing: Option<&Mixing22>,
    phase22: &mut PhaseCoefficients,
    wfhm: &mut HMWaveformStruct,
    wf: &WaveformStruct,
) {
    let ell = wfhm.ell;
    let emm = wfhm.emm;
    phase.phi = [0.; 13];
    phase.phi_l = [0.; 13];
    phase.c0 = 0.;
    phase.c1 = 0.;
    phase.c2 = 0.;
    phase.c3 = 0.;
    phase.c4 = 0.;
    phase.c_l = 0.;
    phase.alpha0 = 0.;
    phase.alpha2 = 0.;
    phase.alpha_l = 0.;

    let n_inter = wfhm.n_colloc_pts_inter_phase;
    let eta_m1 = 1. / wf.eta;
    intermediate_colloc_pts_freqs(phase, wfhm, wf);
    for i in 0..6 {
        phase.colloc_values_inter[i] =
            hm_intermediate::phase_fit(wfhm.mode_int, i, wf) + wfhm.delta_t;
    }
    let fcut_rd = phase.f_phase_match_im;
    let fcut_insp = phase.f_phase_match_in;

    /* Inspiral: rescale the 22 coefficients, phi_lm(f) ~ m/2 phi_22(2 f / m) */
    let phenx_non_log = [
        phase22.phi0,
        phase22.phi1,
        phase22.phi2,
        phase22.phi3,
        phase22.phi4,
        phase22.phi5,
        phase22.phi6,
        phase22.phi7,
        phase22.phi8,
        phase22.phi9,
        0.,
        0.,
        0.,
    ];
    let phenx_log = [
        0.,
        0.,
        0.,
        0.,
        0.,
        phase22.phi5l,
        phase22.phi6l,
        0.,
        phase22.phi8l,
        phase22.phi9l,
        0.,
        0.,
        0.,
    ];
    let pseudo_pn = [
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        phase22.sigma1,
        phase22.sigma2,
        phase22.sigma3,
        phase22.sigma4,
        phase22.sigma5,
    ];
    let m_over_2 = emm as f64 * 0.5;
    let two_over_m = 1. / m_over_2;
    let fact = phase22.phi_norm / wf.eta;
    let log_m_over_2 = m_over_2.ln();
    for i in 0..13 {
        let non_log = phenx_non_log[i] * fact;
        let log = phenx_log[i] * fact;
        let pseudo = pseudo_pn[i] * fact;
        let scale = m_over_2.powf((8 - i as i32) as f64 / 3.);
        phase.phi[i] = fma(-log, log_m_over_2, non_log + pseudo) * scale;
        phase.phi_l[i] = log * scale;
    }
    phase.lambda_pn = if wf.eta > 0.01 {
        hm_inspiral::lambda_pn(wf.eta, wfhm.mode_tag)
    } else {
        hm_inspiral::phase_lambda_fit(wfhm.mode_int, wf)
    };
    phase.phi[8] += phase.lambda_pn;

    /* Intermediate-ringdown */
    // For high spins, extend the inspiral of the 21 into the intermediate
    // region. Like LAL, this uses the 22 phase derivative before the 22
    // connection coefficients are computed (they are zero).
    if wfhm.mode_tag == 21 && wf.s_tot_r >= 0.8 {
        let mut insp_vals = [0.; 3];
        for (i, value) in insp_vals.iter_mut().enumerate() {
            let ff = two_over_m * phase.colloc_freqs_inter[i];
            let powers = UsefulPowers::new(ff);
            *value = 1. / wf.eta * dphase_22(ff, &powers, phase22);
        }
        let diff12 = insp_vals[0] - insp_vals[1];
        let diff23 = insp_vals[1] - insp_vals[2];
        phase.colloc_values_inter[1] = phase.colloc_values_inter[2] + diff23;
        phase.colloc_values_inter[0] = phase.colloc_values_inter[1] + diff12;
    }

    // Choose the collocation points according to spin and mass ratio.
    let mut cpoints = [0usize, 1, 0, 0, 5];
    if wf.eta < wfhm.eta_emr
        || (emm == ell && wf.s_tot_r >= 0.8)
        || (wfhm.mode_tag == 33 && wf.s_tot_r < 0.)
    {
        cpoints[2] = 3;
        cpoints[3] = 4;
    } else if wf.s_tot_r >= 0.8 && wfhm.mode_tag == 21 {
        cpoints[2] = 2;
        cpoints[3] = 4;
    } else {
        cpoints[2] = 2;
        cpoints[3] = 3;
    }

    let fdamp = wfhm.f_damp;
    let fring = wfhm.f_ring;
    let row = |ff: f64| -> [f64; 6] {
        let ffm1 = 1. / ff;
        let ffm2 = ffm1 * ffm1;
        // pow(x, 2) calls: not contracted.
        let lorentzian = fdamp / (fdamp * fdamp + (ff - fring) * (ff - fring));
        [1., lorentzian, ffm1, ffm2, ffm2 * ffm2, ffm1 * ffm2]
    };

    let x: [f64; 6] = if !wfhm.mixing_on {
        // Ringdown: rescale the 22 ringdown, alpha0 + fRD^2 alpha2 / f^2 +
        // alphaL fdamp / (fdamp^2 + (f - fRD)^2).
        let wlm = if ell == emm { 2. } else { emm as f64 / 3. };
        phase.alpha2 = 1. / (fring * fring) * wlm * hm_ringdown::phase_22_alpha2(wf);
        phase.alpha_l = eta_m1 * hm_ringdown::phase_22_alpha_l(wf);
        let powers = UsefulPowers::new(fcut_rd);
        phase.phi0_rd = hm_ringdown::phase_ansatz_int(fcut_rd, &powers, wfhm, phase);
        phase.dphi0_rd = hm_ringdown::phase_ansatz(fcut_rd, &powers, wfhm, phase);

        let mut matrix = [[0.; 5]; 5];
        let mut b = [0.; 5];
        for i in 0..5 {
            let ind = cpoints[i];
            b[i] = phase.colloc_values_inter[ind];
            let powers = row(phase.colloc_freqs_inter[ind]);
            matrix[i].copy_from_slice(&powers[..5]);
        }
        let x = lu_solve(matrix, b);
        [x[0], x[1], x[2], x[3], x[4], 0.]
    } else {
        // Derivatives of the rotated spheroidal ringdown phase at fcut_rd
        // by finite differences.
        let mut phase_args = [0.; 3];
        for (i, arg) in phase_args.iter_mut().enumerate() {
            let ff = fma((i as i32 - 1) as f64, FSTEP_MIXING, fcut_rd);
            let powers = UsefulPowers::new(ff);
            let h = rotated_32(&powers, *amp, phase, mixing, wfhm);
            *arg = h.im.atan2(h.re) % (2. * PI);
            // Make sure all the points belong to the same branch.
            if *arg > 0. {
                *arg -= 2. * PI;
            }
        }
        phase.phi0_rd = phase_args[1];
        let fstep_m1 = 1. / FSTEP_MIXING;
        phase.dphi0_rd = 0.5 * fstep_m1 * (phase_args[2] - phase_args[0]);
        let d2phi0_rd =
            fstep_m1 * fstep_m1 * (fma(-2., phase_args[1], phase_args[2]) + phase_args[0]);

        if wf.eta > wfhm.eta_emr {
            phase.colloc_freqs_inter[n_inter - 2] = fcut_rd;
            phase.colloc_values_inter[n_inter - 2] = phase.dphi0_rd;
        }
        phase.colloc_freqs_inter[n_inter - 1] = fcut_rd;
        phase.colloc_values_inter[n_inter - 1] = d2phi0_rd;

        let mut matrix = [[0.; 6]; 6];
        let mut b = [0.; 6];
        for i in 0..6 {
            b[i] = phase.colloc_values_inter[i];
            matrix[i] = row(phase.colloc_freqs_inter[i]);
        }
        // The last point is the second derivative of the rotated ansatz.
        let ff = phase.colloc_freqs_inter[5];
        let ffm1 = 1. / ff;
        let ffm2 = ffm1 * ffm1;
        let ffm3 = ffm2 * ffm1;
        let ffm4 = ffm2 * ffm2;
        let ffm5 = ffm3 * ffm2;
        let den = fdamp * fdamp + (ff - fring) * (ff - fring);
        matrix[5] = [
            0.,
            -2. * fdamp * (ff - fring) / (den * den),
            -ffm2,
            -2. * ffm3,
            -4. * ffm5,
            -3. * ffm4,
        ];
        lu_solve(matrix, b)
    };

    phase.c0 = x[0];
    phase.c_l = x[1];
    phase.c1 = x[2];
    phase.c2 = x[3];
    phase.c4 = x[4];
    if wfhm.mode_tag == 32 {
        phase.c3 = x[5];
        // For extreme mass ratios glue the intermediate and ringdown regions.
        if wf.eta < wfhm.eta_emr {
            let powers = UsefulPowers::new(fcut_rd);
            phase.c0 = phase.c0 + phase.dphi0_rd
                - hm_intermediate::phase_ansatz(fcut_rd, &powers, wfhm, phase);
        }
    }

    // Glue inspiral and intermediate.
    let powers = UsefulPowers::new(fcut_insp);
    phase.c1_insp = hm_intermediate::phase_ansatz(fcut_insp, &powers, wfhm, phase)
        - hm_inspiral::phase_ansatz(fcut_insp, &powers, phase);
    phase.c_insp = fma(
        -phase.c1_insp,
        fcut_insp,
        hm_intermediate::phase_ansatz_int(fcut_insp, &powers, wfhm, phase),
    ) - hm_inspiral::phase_ansatz_int(fcut_insp, &powers, phase);
    // Glue intermediate and ringdown.
    let powers = UsefulPowers::new(fcut_rd);
    phase.c1_rd = hm_intermediate::phase_ansatz(fcut_rd, &powers, wfhm, phase) - phase.dphi0_rd;
    phase.c_rd = fma(
        -phase.c1_rd,
        fcut_rd,
        hm_intermediate::phase_ansatz_int(fcut_rd, &powers, wfhm, phase),
    ) - phase.phi0_rd;

    // Align each mode so that at low frequency its phase relative to the 22
    // agrees with PN, Eq. 4.13 of arXiv:2001.10914.
    let falign = if wf.eta > wfhm.eta_emr {
        0.6 * m_over_2 * wf.f_meco
    } else {
        m_over_2 * wf.f_meco
    };
    let powers_of_falign = UsefulPowers::new(falign);
    let powers_of_f = UsefulPowers::new(two_over_m * falign);

    if !wfhm.mixing_on {
        let powers_of_mf_ref = UsefulPowers::new(wf.mf_ref);
        phase_22_connection_coefficients(phase22);
        wfhm.timeshift = time_shift_22(phase22, wf);
        wfhm.phiref22 = phi_ref_22(wf, phase22, wfhm, &powers_of_mf_ref);
    }

    let phi22 = fma(
        1. / wf.eta,
        phase_22(two_over_m * falign, &powers_of_f, phase22),
        wfhm.phaseshift,
    ) + wfhm.phiref22;
    let phi_lm = fma(
        phase.c1_insp,
        falign,
        hm_inspiral::phase_ansatz_int(falign, &powers_of_falign, phase),
    ) + phase.c_insp;
    let deltaphi_lm = fma(
        -(3. * FRAC_PI_4),
        1. - m_over_2,
        fma(m_over_2, phi22, wfhm.timeshift * falign),
    ) - phi_lm;
    phase.deltaphi_lm = deltaphi_lm % (2. * PI);

    // The 21 amplitude changes sign across the parameter space: the model
    // amplitude is positive by construction, so flip the phase by pi.
    if wfhm.mode_tag == 21 && pn21_amp_sign(0.008, wf) > 0 {
        phase.deltaphi_lm += PI;
    }
}

/// `phiref22`: the phase offset `phifRef` of the 22 mode as computed in
/// `IMRPhenomXHM_GetPhaseCoefficients` and `GetSpheroidalCoefficients`.
fn phi_ref_22(
    wf: &WaveformStruct,
    phase22: &PhaseCoefficients,
    wfhm: &HMWaveformStruct,
    powers_of_mf_ref: &UsefulPowers,
) -> f64 {
    fma(
        -1. / wf.eta,
        phase_22(wf.mf_ref, powers_of_mf_ref, phase22),
        -(wfhm.timeshift * wf.mf_ref),
    ) - wfhm.phaseshift
        + 2.0 * wf.phi0
        + FRAC_PI_4
}

/// `GetSpheroidalCoefficients` (ringdown phase version 122019): the
/// spheroidal ringdown phase of the (3,2) and its alignment with the 22.
fn get_spheroidal_coefficients(
    phase: &mut HMPhaseCoefficients,
    phase22: &mut PhaseCoefficients,
    wfhm: &mut HMWaveformStruct,
    wf: &WaveformStruct,
) {
    // IMRPhenomXHM_Ringdown_CollocPtsFreqs (122019)
    let fringlm = wfhm.f_ring;
    let fdamplm = wfhm.f_damp;
    phase.colloc_freqs_rd = [
        wf.f_ring,
        fma(-1.5, fdamplm, fringlm),
        fma(-0.5, fdamplm, fringlm),
        fma(0.5, fdamplm, fringlm),
    ];
    let mut matrix = [[0.; 4]; 4];
    let mut b = [0.; 4];
    for i in 0..4 {
        phase.colloc_values_rd[i] = hm_ringdown::phase_32_fit(i, wf);
        b[i] = phase.colloc_values_rd[i];
        let ff = phase.colloc_freqs_rd[i];
        let ffm1 = 1. / ff;
        let ffm2 = ffm1 * ffm1;
        let lorentzian = fdamplm / (fdamplm * fdamplm + (ff - fringlm) * (ff - fringlm));
        matrix[i] = [1., lorentzian, ffm2, ffm2 * ffm2];
    }
    let [alpha0_s, alpha_l_s, alpha2_s, alpha4_s] = lu_solve(matrix, b);
    phase.alpha0_s = alpha0_s;
    phase.alpha_l_s = alpha_l_s;
    phase.alpha2_s = alpha2_s;
    phase.alpha4_s = alpha4_s;

    // Time shift of the spheroidal ansatz.
    let fref_rd = wf.f_ring + wf.f_damp;
    let powers_of_fref = UsefulPowers::new(fref_rd);
    let tshift = hm_ringdown::phase_32_spheroidal_time_shift(wf);
    phase_22_connection_coefficients(phase22);
    wfhm.timeshift = time_shift_22(phase22, wf);
    phase.phi0_s = 0.;
    let dphi22_ref = fma(
        1. / wf.eta,
        dphase_22(fref_rd, &powers_of_fref, phase22),
        wfhm.timeshift,
    );
    phase.alpha0_s = phase.alpha0_s + dphi22_ref + tshift
        - hm_ringdown::phase_ansatz(fref_rd, &powers_of_fref, wfhm, phase);

    // Phase shift of the spheroidal ansatz.
    let fref_rd = wf.f_ring;
    let powers_of_fref = UsefulPowers::new(fref_rd);
    let powers_of_mf_ref = UsefulPowers::new(wf.mf_ref);
    wfhm.phiref22 = phi_ref_22(wf, phase22, wfhm, &powers_of_mf_ref);
    let phi22_ref = fma(
        1. / wf.eta,
        phase_22(fref_rd, &powers_of_fref, phase22),
        wfhm.timeshift * fref_rd,
    ) + wfhm.phaseshift
        + wfhm.phiref22;
    let phishift = hm_ringdown::phase_32_spheroidal_phase_shift(wf);
    phase.phi0_s =
        phi22_ref - hm_ringdown::phase_ansatz_int(fref_rd, &powers_of_fref, wfhm, phase) + phishift;
}
