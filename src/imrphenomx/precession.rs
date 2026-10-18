//! Precession variables of IMRPhenomXP/XPHM: the parts of
//! `LALSimIMRPhenomX_precession.c` that set up `IMRPhenomXPrecessionStruct`
//! (`IMRPhenomXGetAndSetPrecessionVariables`,
//! `IMRPhenomX_SetPrecessingRemnantParams`), the frame conventions relating
//! the source (L) frame to the J frame, and the NNLO post-Newtonian Euler
//! angles (`PhenomXPrecVersion` 102). The MSA angles (223) are in
//! [`crate::imrphenomx::msa`].
//!
//! Only the default settings of LALSimulation are supported:
//! `PhenomXPConvention` 1, `PhenomXPFinalSpinMod` 4, no PNR angles or
//! coprecessing deviations, symmetric waveforms and no multibanding.
//!
//! The frame construction and the final spin mirror the fused multiply-adds
//! of the arm64 LALSimulation build (see [`crate::spintaylor::fp`]), which
//! the SpinTaylor angles need to agree with LAL bit for bit. The NNLO angle
//! coefficients are evaluated without fused operations and agree to a few
//! units in the last place.

use std::f64::consts::{FRAC_PI_2, PI};

use crate::constants::{C_SI, G_SI, MSUN_SI, MTSUN_SI};
use crate::grid::{FrequencyGrid, ScalePowers};
use crate::imrphenomx::internals::{IMRPhenomXError, WaveformStruct, POWERS_OF_LALPI};
use crate::imrphenomx::msa::{MsaInputs, MsaSystem};
use crate::imrphenomx::utilities::final_spin_2017;
use crate::imrphenomx::wigner::SphericalHarmonics;
use crate::spintaylor::fp::{fma, fms, fnms};

/// `MAX_TOL_ATAN` of `LALSimIMRPhenomX_precession.h`.
pub(crate) const MAX_TOL_ATAN: f64 = 1.0e-15;

/// The Euler angle prescription, LALSimulation's `PhenomXPrecVersion`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PrecessionVersion {
    /// `PhenomXPrecVersion = 102`: NNLO single-spin post-Newtonian Euler
    /// angles with the 3PN orbital angular momentum
    /// (`IMRPhenomX_PN_Euler_alpha_NNLO`, `IMRPhenomX_PN_Euler_epsilon_NNLO`).
    Nnlo,
    /// `PhenomXPrecVersion = 320` (the default here): Euler angles from the
    /// numerical integration of the SpinTaylorT4 precession equations, with
    /// an analytical continuation through merger and ringdown. Like
    /// LALSimulation, binaries with in-plane spins below `1e-7` use the
    /// NNLO angles instead. Where the PN integration does not provide the
    /// angles (it stops before reaching the first frequency, typically for
    /// total masses above ~150 solar masses), LALSimulation silently uses
    /// the MSA angles (223) instead, and so do we.
    #[default]
    SpinTaylor,
    /// `PhenomXPrecVersion = 223`: double-spin MSA Euler angles
    /// (Chatziioannou et al. 2017, arXiv:1703.03967) as implemented in
    /// `LALSimInspiralFDPrecAngles`, with the precession-averaged final
    /// spin. Like LALSimulation, binaries for which the MSA system cannot
    /// be initialised use the NNLO angles (102) instead.
    Msa,
}

impl PrecessionVersion {
    /// The LALSimulation `PhenomXPrecVersion` flag.
    pub fn lal_flag(&self) -> u32 {
        match self {
            Self::Nnlo => 102,
            Self::SpinTaylor => 320,
            Self::Msa => 223,
        }
    }
}

/// Parameters of a precessing binary in the conventions of
/// `SimInspiralChooseFDWaveform`: the same type as [`crate::model::Binary`].
pub type PrecessingBinary = crate::model::Binary;

impl PrecessingBinary {
    /// `XLALIMRPhenomXPCheckMassesAndSpins`: masses (solar masses) and spins
    /// ordered such that the first body is strictly the heavier one. Note
    /// that, unlike `IMRPhenomXSetWaveformVariables`, equal masses swap the
    /// bodies.
    pub(crate) fn ordered(&self) -> (f64, f64, [f64; 3], [f64; 3]) {
        // Compare the masses in kg as LAL does.
        if self.mass_1 * MSUN_SI > self.mass_2 * MSUN_SI {
            (self.mass_1, self.mass_2, self.spin_1, self.spin_2)
        } else {
            (self.mass_2, self.mass_1, self.spin_2, self.spin_1)
        }
    }
}

/// Rotate `v` about the `z` axis by `angle`, `IMRPhenomX_rotate_z`.
#[inline]
pub(crate) fn rotate_z(angle: f64, v: [f64; 3]) -> [f64; 3] {
    let (sina, cosa) = (angle.sin(), angle.cos());
    rotate_z_cs(cosa, sina, v)
}

/// [`rotate_z`] with precomputed `cos(angle)`, `sin(angle)`.
#[inline]
pub(crate) fn rotate_z_cs(cosa: f64, sina: f64, v: [f64; 3]) -> [f64; 3] {
    [
        fms(v[0], cosa, v[1] * sina),
        fma(v[0], sina, v[1] * cosa),
        v[2],
    ]
}

/// Rotate `v` about the `y` axis by `angle`, `IMRPhenomX_rotate_y`.
#[inline]
pub(crate) fn rotate_y(angle: f64, v: [f64; 3]) -> [f64; 3] {
    let (sina, cosa) = (angle.sin(), angle.cos());
    rotate_y_cs(cosa, sina, v)
}

/// [`rotate_y`] with precomputed `cos(angle)`, `sin(angle)`.
#[inline]
pub(crate) fn rotate_y_cs(cosa: f64, sina: f64, v: [f64; 3]) -> [f64; 3] {
    [
        fma(v[0], cosa, v[2] * sina),
        v[1],
        fma(-v[0], sina, v[2] * cosa),
    ]
}

/// `XLALSimIMRPhenomXatan2tol`: `atan2(a, b)`, zero if both arguments are
/// below `tol` in magnitude.
#[inline]
pub(crate) fn atan2_tol(a: f64, b: f64, tol: f64) -> f64 {
    if a.abs() < tol && b.abs() < tol {
        0.0
    } else {
        a.atan2(b)
    }
}

/// `XLALSimIMRPhenomXUtilsMftoHz`.
#[inline]
pub(crate) fn mf_to_hz(mf: f64, m_tot: f64) -> f64 {
    mf / (MTSUN_SI * m_tot)
}

/// `XLALSimIMRPhenomXUtilsHztoMf`.
#[inline]
pub(crate) fn hz_to_mf(f: f64, m_tot: f64) -> f64 {
    f * (MTSUN_SI * m_tot)
}

/// `XLALSimIMRPhenomXPrecessingFinalSpin2017`: final spin of a precessing
/// binary from the aligned-spin fit and the in-plane spin `chi_inplane`
/// (weighted by `(m1 / M)^2`), carrying the sign of the aligned-spin fit.
pub(crate) fn precessing_final_spin_2017(
    eta: f64,
    chi1l: f64,
    chi2l: f64,
    chi_inplane: f64,
) -> f64 {
    let m1 = 0.5 * (1.0 + fnms(4.0, eta, 1.0).sqrt());
    let m2 = 0.5 * (1.0 - fnms(4.0, eta, 1.0).sqrt());
    let m = m1 + m2;
    // IMRPhenomX_InternalNudge takes its argument by value: no effect.
    let (q_factor, af_parallel) = if m1 >= m2 {
        (m1 / m, final_spin_2017(eta, chi1l, chi2l))
    } else {
        (m2 / m, final_spin_2017(eta, chi2l, chi1l))
    };
    let s_perp = chi_inplane * q_factor * q_factor;
    1.0f64.copysign(af_parallel) * fma(s_perp, s_perp, af_parallel * af_parallel).sqrt()
}

/// `XLALSimInspiralChirpTimeBound`: upper bound on the chirp time (s) from
/// `f_start` (Hz) for masses in kg and aligned spins.
pub(crate) fn chirp_time_bound(f_start: f64, m1: f64, m2: f64, s1: f64, s2: f64) -> f64 {
    let m = m1 + m2;
    let mu = m1 * m2 / m;
    let eta = mu / m;
    let chi = (if s1.abs() > s2.abs() { s1 } else { s2 }).abs();
    // XLALSimInspiralTaylorT2Timing_0PNCoeff
    let total_mass = m * (G_SI / C_SI.powf(3.0));
    let c0 = (-5. * total_mass / (256. * eta)).abs();
    let c2 = 7.43 / 2.52 + 11. / 3. * eta;
    let c3 = (226.0 / 15.0) * chi;
    let c4 = 30.58673 / 5.08032 + 54.29 / 5.04 * eta + 61.7 / 7.2 * eta * eta;
    let v = (PI * G_SI * m * f_start).cbrt() / C_SI;
    c0 * v.powf(-8.0) * (1.0 + (c2 + (c3 + c4 * v) * v) * v * v)
}

/// `get_deltaF_from_wfstruct`: dimensionless frequency spacing of the
/// default uniform grid for the binary, used by the SpinTaylor angles when
/// the waveform is evaluated on a non-uniform grid.
pub(crate) fn delta_mf_from_wfstruct(wf: &WaveformStruct) -> f64 {
    let seglen = chirp_time_bound(wf.f_ref, wf.m1_si, wf.m2_si, wf.chi1l, wf.chi2l);
    let delta_f_v1 = 1. / 4.0f64.max(2.0f64.powf((seglen.ln() / 2.0f64.ln()).ceil()));
    let delta_f = delta_f_v1.min(0.1);
    hz_to_mf(delta_f, wf.m_tot)
}

/// Coefficients of the PN orbital angular momentum, `pPrec->L0 ... L8L`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct OrbitalAngularMomentum {
    pub(crate) l0: f64,
    pub(crate) l1: f64,
    pub(crate) l2: f64,
    pub(crate) l3: f64,
    pub(crate) l4: f64,
    pub(crate) l5: f64,
    pub(crate) l6: f64,
    pub(crate) l7: f64,
    pub(crate) l8: f64,
    pub(crate) l8l: f64,
}

impl OrbitalAngularMomentum {
    /// The 3PN orbital angular momentum with linear-in-spin terms used by
    /// PrecVersions 102 and 3xx.
    fn three_pn(eta: f64, delta: f64, chi1l: f64, chi2l: f64) -> Self {
        Self {
            l0: 1.0,
            l1: 0.0,
            l2: 3.0 / 2. + eta / 6.0,
            l3: (5.0
                * fma(
                    chi1l,
                    fnms(2.0, delta, -2.0) + eta,
                    chi2l * (fma(2.0, delta, -2.0) + eta),
                ))
                / 6.,
            l4: fma(-57.0 + eta, eta, 81.0) / 24.,
            l5: (-7.0
                * fma(
                    chi1l,
                    fma(
                        eta,
                        fma(2.0, eta, -121.0),
                        fma(delta, fnms(31.0, eta, 72.0), 72.0),
                    ),
                    chi2l
                        * fma(
                            delta,
                            fma(31.0, eta, -72.0),
                            fma(eta, fma(2.0, eta, -121.0), 72.0),
                        ),
                ))
                / 144.,
            l6: fma(
                eta,
                fma(
                    2214.0,
                    POWERS_OF_LALPI.two,
                    fma(eta, fma(7.0, eta, 1674.0), -62001.0),
                ),
                10935.0,
            ) / 1296.,
            l7: 0.0,
            l8: 0.0,
            l8l: 0.0,
        }
    }

    /// The 3PN orbital angular momentum of versions 222 and 223, with the
    /// spin-orbit terms of `LALSimInspiralFDPrecAngles`.
    fn three_pn_msa(eta: f64, delta: f64, chi1l: f64, chi2l: f64) -> Self {
        let eta2 = eta * eta;
        // chi1L + chi2L + chi1L delta - chi2L delta
        let a = fnms(chi2l, delta, fma(chi1l, delta, chi1l + chi2l));
        Self {
            l3: fma(-7.0, a, 5.0 * (chi1l + chi2l) * eta) / 6.,
            l5: fma(
                28.0 * (chi1l + chi2l),
                eta2,
                fma(
                    511.0 * (chi1l - chi2l) * delta,
                    eta,
                    fma(-1650.0, a, 1336.0 * (chi1l + chi2l) * eta),
                ),
            ) / 600.,
            l6: fma(
                eta,
                fma(
                    2214.0,
                    POWERS_OF_LALPI.two,
                    fma(7.0, eta2, fma(1674.0, eta, -62001.0)),
                ),
                10935.0,
            ) / 1296.,
            ..Self::three_pn(eta, delta, chi1l, chi2l)
        }
    }

    /// `XLALSimIMRPhenomXLPNAnsatz`: `LNorm sum_a L_a x^(a/2)` with
    /// `x = v^2`.
    pub(crate) fn ansatz(&self, v: f64, l_norm: f64) -> f64 {
        let x = v * v;
        let x2 = x * x;
        let x3 = x * x2;
        let x4 = x * x3;
        let sqx = x.sqrt();
        let mut sum = fma(self.l1, sqx, self.l0);
        sum = fma(self.l2, x, sum);
        sum = fma(self.l3, x * sqx, sum);
        sum = fma(self.l4, x2, sum);
        sum = fma(self.l5, x2 * sqx, sum);
        sum = fma(self.l6, x3, sum);
        sum = fma(self.l7, x3 * sqx, sum);
        sum = fma(self.l8, x4, sum);
        sum = fma(self.l8l * x4, x.ln(), sum);
        l_norm * sum
    }

    /// [`Self::ansatz`] for the 3PN coefficients, skipping the vanishing
    /// terms (which leaves the result unchanged bit for bit).
    #[inline]
    pub(crate) fn ansatz_3pn(&self, v: f64, l_norm: f64) -> f64 {
        let x = v * v;
        let x2 = x * x;
        let x3 = x * x2;
        let sqx = x.sqrt();
        let mut sum = fma(self.l2, x, self.l0);
        sum = fma(self.l3, x * sqx, sum);
        sum = fma(self.l4, x2, sum);
        sum = fma(self.l5, x2 * sqx, sum);
        sum = fma(self.l6, x3, sum);
        l_norm * sum
    }
}

/// Coefficients of the NNLO PN Euler angles `alpha` and `epsilon`
/// (`pPrec->alpha1 ... alpha5`, `epsilon1 ... epsilon5`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct NnloCoefficients {
    pub(crate) alpha1: f64,
    pub(crate) alpha2: f64,
    pub(crate) alpha3: f64,
    pub(crate) alpha4l: f64,
    pub(crate) alpha5: f64,
    pub(crate) epsilon1: f64,
    pub(crate) epsilon2: f64,
    pub(crate) epsilon3: f64,
    pub(crate) epsilon4l: f64,
    pub(crate) epsilon5: f64,
}

/// Equivalent of the parts of `IMRPhenomXPrecessionStruct` used by
/// IMRPhenomXP/XPHM with the default settings, built by
/// [`PrecessionStruct::new`] (`IMRPhenomXGetAndSetPrecessionVariables`).
///
/// Masses are in units of the total mass and the body labels follow
/// `XLALIMRPhenomXPCheckMassesAndSpins` (body 1 strictly heavier, so equal
/// masses are swapped with respect to the aligned-spin waveform struct).
#[derive(Clone, Debug)]
#[allow(dead_code)] // the LAL struct is kept complete for XPHM
pub(crate) struct PrecessionStruct {
    /// Effective `PhenomXPrecVersion` (after the fallbacks to 102 for
    /// negligible in-plane spins with the SpinTaylor angles or a failed
    /// initialisation of the MSA system).
    pub(crate) version: u32,
    /// First digit of the version: 1 NNLO, 2 MSA, 3 SpinTaylor. As in LAL,
    /// this stays 2 when the MSA system falls back to the NNLO angles.
    pub(crate) precessing_tag: u32,
    /// Mass of the larger body in kg as passed to LAL.
    pub(crate) m1_si: f64,
    pub(crate) m2_si: f64,
    pub(crate) m1: f64,
    pub(crate) m2: f64,
    pub(crate) m1_2: f64,
    pub(crate) m2_2: f64,
    /// `pi G M / c^3` in seconds.
    pub(crate) pi_gm: f64,
    /// Dimensionless spins (source frame).
    pub(crate) chi1: [f64; 3],
    pub(crate) chi2: [f64; 3],
    /// Effective precession spin.
    pub(crate) chi_p: f64,
    /// Norm of the total in-plane spin, normalised by `m1^2`.
    pub(crate) chi_tot_perp: f64,
    /// Aligned spin angular momentum `S.L` in units of `M^2`.
    pub(crate) sl: f64,
    /// In-plane spin `chi_p m1^2`.
    pub(crate) sperp: f64,
    /// Orbital angular momentum coefficients.
    pub(crate) l: OrbitalAngularMomentum,
    /// Orbital angular momentum at the reference frequency.
    pub(crate) l_ref: f64,
    /// Total angular momentum at the reference frequency (source frame).
    pub(crate) j0: [f64; 3],
    pub(crate) j0_norm: f64,
    /// Polar and azimuthal angle of J in the source frame.
    pub(crate) theta_j_sf: f64,
    pub(crate) phi_j_sf: f64,
    /// Third rotation angle from the source frame to the J frame.
    pub(crate) kappa: f64,
    /// Angle between J and the line of sight.
    pub(crate) theta_jn: f64,
    /// Initial values of the Euler angles alpha and epsilon.
    pub(crate) alpha0: f64,
    pub(crate) epsilon0: f64,
    /// Angle rotating the polarisations from the J-frame (P, Q) basis to the
    /// LAL wave frame.
    pub(crate) zeta_polarization: f64,
    /// NNLO Euler angle coefficients (zero for the SpinTaylor versions).
    pub(crate) nnlo: NnloCoefficients,
    /// The MSA system (version 223 only).
    pub(crate) msa: Option<MsaSystem>,
    /// NNLO or MSA offsets at the reference frequency, including `alpha0`,
    /// `epsilon0`.
    pub(crate) alpha_offset: f64,
    pub(crate) epsilon_offset: f64,
    /// `-2` spin-weighted spherical harmonics at `(thetaJN, 0)`.
    pub(crate) ylm: SphericalHarmonics,
    /// Smallest and largest `|m|` of the co-precessing modes, which set the
    /// frequency range of the SpinTaylor angles (2 and 2 for XP).
    pub(crate) m_min: u32,
    pub(crate) m_max: u32,
    /// Frequency buffer (Hz) below `f_min` for the SpinTaylor integration.
    pub(crate) integration_buffer: f64,
    /// Aligned-spin (2,2) ringdown and damping frequencies, which set the
    /// end of the PN integration (`IMRPhenomX_InspiralAngles_SpinTaylor` is
    /// called before the precessing final spin is set).
    pub(crate) f_ring_aligned: f64,
    pub(crate) f_damp_aligned: f64,
}

impl PrecessionStruct {
    /// `IMRPhenomXGetAndSetPrecessionVariables` for `PhenomXPConvention` 1
    /// and `PhenomXPFinalSpinMod` 4, except for the numerical integration
    /// of the SpinTaylor equations (version 320), which depends on the
    /// frequency grid and is done separately by
    /// [`crate::imrphenomx::spintaylor_angles`]. If that integration fails,
    /// LAL switches to version 223 and carries on with the same structs;
    /// the SpinTaylor set-up done before (`deltaMF`, the integration
    /// buffer, the PN orders inserted in the LALDict) does not enter the
    /// MSA angles, so the fallback is equivalent to calling this function
    /// again with [`PrecessionVersion::Msa`].
    ///
    /// `wf` is modified as LAL modifies `pWF`: `phi0 = 0`, the precessing
    /// final spin with the corresponding (2,2) QNM frequencies and, for the
    /// SpinTaylor angles on a non-uniform grid, `delta_mf`.
    ///
    /// `m1_si`, `m2_si` (kg) and the spins must be ordered with
    /// [`PrecessingBinary::ordered`].
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        wf: &mut WaveformStruct,
        m1_si: f64,
        m2_si: f64,
        chi1: [f64; 3],
        chi2: [f64; 3],
        requested: PrecessionVersion,
        m_min: u32,
        m_max: u32,
    ) -> Result<Self, IMRPhenomXError> {
        let [chi1x, chi1y, chi1z] = chi1;
        let [chi2x, chi2y, chi2z] = chi2;

        let mut version = requested.lal_flag();
        // Default to NNLO angles if the in-plane spins are negligible: the
        // SpinTaylor solutions would be dominated by numerical noise.
        let chi_in_plane = fma(
            chi2y,
            chi2y,
            fma(chi2x, chi2x, fma(chi1x, chi1x, chi1y * chi1y)),
        )
        .sqrt();
        if chi_in_plane < 1e-7 && version == 320 {
            version = 102;
        }

        let m1 = m1_si / wf.m_tot_si;
        let m2 = m2_si / wf.m_tot_si;
        let m = m1 + m2;
        let m1_2 = m1 * m1;
        let m1_3 = m1 * m1_2;
        let m1_4 = m1 * m1_3;
        let m1_5 = m1 * m1_4;
        let m1_6 = m1 * m1_5;
        let m1_7 = m1 * m1_6;
        let m1_8 = m1 * m1_7;
        let m2_2 = m2 * m2;
        let q = m1 / m2;

        let eta = wf.eta;
        let delta = wf.delta;
        let chi_eff = wf.chi_eff;

        let pi_gm = PI * (m1_si + m2_si) * (G_SI / C_SI) / (C_SI * C_SI);

        let chi1_norm = fma(chi1z, chi1z, fma(chi1x, chi1x, chi1y * chi1y)).sqrt();
        let chi2_norm = fma(chi2z, chi2z, fma(chi2x, chi2x, chi2y * chi2y)).sqrt();
        for norm in [chi1_norm, chi2_norm] {
            if norm.is_nan() || norm.abs() > 1.0 {
                return Err(IMRPhenomXError::UnphysicalSpinMagnitude(norm));
            }
        }

        let s1x = chi1x * m1_2;
        let s1y = chi1y * m1_2;
        let s2x = chi2x * m2_2;
        let s2y = chi2y * m2_2;
        let s1_perp = m1_2 * fma(chi1x, chi1x, chi1y * chi1y).sqrt();
        let s2_perp = m2_2 * fma(chi2x, chi2x, chi2y * chi2y).sqrt();
        let s_tot_perp = fma(s1x + s2x, s1x + s2x, (s1y + s2y) * (s1y + s2y)).sqrt();
        let chi_tot_perp = s_tot_perp * (m * m) / m1_2;

        // Effective precessing spin (Schmidt et al. 2015)
        let a1 = 2.0 + (3.0 * m2) / (2.0 * m1);
        let a2 = 2.0 + (3.0 * m1) / (2.0 * m2);
        let asp1 = a1 * s1_perp;
        let asp2 = a2 * s2_perp;
        let num = if asp2 > asp1 { asp2 } else { asp1 };
        let den = if m2 > m1 { a2 * m2_2 } else { a1 * m1_2 };
        let chip = num / den;
        let chi1l = chi1z;
        let chi2l = chi2z;
        let sl = fma(chi1l, m1_2, chi2l * m2_2);
        let sperp = chip * m1_2;

        // SpinTaylor set-up done before the integration: the buffer below
        // f_min and, on non-uniform grids, the default frequency spacing.
        let precessing_tag = version / 100;
        let mut integration_buffer = 0.0;
        if precessing_tag == 3 {
            if wf.delta_f == 0. {
                wf.delta_mf = delta_mf_from_wfstruct(wf);
            }
            integration_buffer = if wf.delta_f > 0. {
                3. * wf.delta_f
            } else {
                0.5
            };
        }
        let f_ring_aligned = wf.f_ring;
        let f_damp_aligned = wf.f_damp;

        // IMRPhenomX_Initialize_MSA_System; if it fails (MSA_ERROR), LAL
        // continues with the NNLO angles (102) for versions 220, 223, 224.
        let msa = if version == 223 {
            let msa = MsaSystem::new(&MsaInputs {
                eta,
                m1: wf.m1,
                m2: wf.m2,
                m1_2,
                m2_2,
                chi1,
                chi2,
                pi_gm,
                f_ref: wf.f_ref,
            });
            if msa.is_none() {
                version = 102;
            }
            msa
        } else {
            None
        };

        // IMRPhenomX_SetPrecessingRemnantParams: PhenomXPFinalSpinMod 4 uses
        // the total in-plane spin for the SpinTaylor angles; otherwise it
        // becomes 3, the precession-averaged final spin of the MSA system,
        // which for the other versions (and a failed MSA system) defaults
        // to version 0 (chi_p).
        let mut a_final_prec = match (&msa, precessing_tag) {
            (Some(msa), _) => msa.final_spin(
                final_spin_2017(eta, chi1l, chi2l),
                m,
                m1_2,
                m2_2,
                chi1l,
                chi2l,
            ),
            (None, 3) => precessing_final_spin_2017(eta, chi1l, chi2l, chi_tot_perp),
            (None, _) => precessing_final_spin_2017(eta, chi1l, chi2l, chip),
        };
        let mut a_final = a_final_prec;
        if a_final.abs() > 1.0 {
            a_final = 1.0f64.copysign(a_final);
        }
        if a_final_prec.abs() > 1.0 {
            a_final_prec = 1.0f64.copysign(a_final_prec);
        }
        wf.a_final_prec = a_final_prec;
        wf.set_final_spin(a_final);

        let l = if msa.is_some() {
            OrbitalAngularMomentum::three_pn_msa(eta, delta, chi1l, chi2l)
        } else {
            OrbitalAngularMomentum::three_pn(eta, delta, chi1l, chi2l)
        };
        let l_ref = m * m * l.ansatz(wf.v_ref, wf.eta / wf.v_ref);

        // Source frame J = L + S1 + S2 (L along z).
        let j0 = [
            fma(m1_2, chi1x, m2_2 * chi2x),
            fma(m1_2, chi1y, m2_2 * chi2y),
            fma(m1_2, chi1z, m2_2 * chi2z) + l_ref,
        ];
        let j0_norm = fma(j0[2], j0[2], fma(j0[0], j0[0], j0[1] * j0[1])).sqrt();
        let theta_j_sf = if j0_norm < 1e-10 {
            0.0
        } else {
            (j0[2] / j0_norm).acos()
        };
        let phi_ref = wf.phi_ref_in;
        // Convention 1
        let phi_j_sf = if j0[0].abs() < MAX_TOL_ATAN && j0[1].abs() < MAX_TOL_ATAN {
            0.0
        } else {
            j0[1].atan2(j0[0])
        };
        wf.phi0 = 0.0;

        // kappa: rotate N with RY[-thetaJ].RZ[-phiJ], which brings J to z,
        // and take minus the azimuth of the rotated N.
        let inclination = wf.inclination;
        let n_sf = [
            inclination.sin() * (FRAC_PI_2 - phi_ref).cos(),
            inclination.sin() * (FRAC_PI_2 - phi_ref).sin(),
            inclination.cos(),
        ];
        let tmp = rotate_y(-theta_j_sf, rotate_z(-phi_j_sf, n_sf));
        let kappa = atan2_tol(tmp[1], tmp[0], MAX_TOL_ATAN);
        let alpha0 = PI - kappa;

        let j0_dot_n = fma(j0[2], n_sf[2], fma(j0[0], n_sf[0], j0[1] * n_sf[1]));
        let theta_jn = (j0_dot_n / j0_norm).acos();
        let nz_jf = theta_jn.cos();
        let nx_jf = theta_jn.sin();

        // Polarisation basis: X of the LAL wave frame in the J frame, and its
        // angle in the P, Q basis of Arun et al.
        let x_sf = [
            -inclination.cos() * phi_ref.sin(),
            -inclination.cos() * phi_ref.cos(),
            inclination.sin(),
        ];
        let tmp = rotate_z(-kappa, rotate_y(-theta_j_sf, rotate_z(-phi_j_sf, x_sf)));
        let p_arun = [nz_jf, 0.0, -nx_jf];
        let q_arun = [0.0, 1.0, 0.0];
        let x_dot_p = fma(
            tmp[2],
            p_arun[2],
            fma(tmp[0], p_arun[0], tmp[1] * p_arun[1]),
        );
        let x_dot_q = fma(
            tmp[2],
            q_arun[2],
            fma(tmp[0], q_arun[0], tmp[1] * q_arun[1]),
        );
        let zeta_polarization = x_dot_q.atan2(x_dot_p);

        let nnlo = if version / 100 == 1 {
            nnlo_coefficients(
                eta, delta, q, chi_eff, chip, m1, m1_2, m1_3, m1_4, m1_5, m1_6, m1_8,
            )
        } else {
            NnloCoefficients::default()
        };

        let epsilon0 = phi_j_sf - PI;

        // Get_alphaepsilon_atfref with mprime = 2.
        let omega_ref = wf.pi_m * wf.f_ref * 2. / 2.;
        let (alpha_offset, epsilon_offset) = match &msa {
            Some(msa) => {
                let angles = msa.phi_zeta_costhetal(omega_ref.cbrt());
                (angles.phiz - alpha0, angles.zeta - epsilon0)
            }
            None => {
                let logomega_ref = omega_ref.ln();
                let omega_ref_cbrt = omega_ref.cbrt();
                let omega_ref_cbrt2 = omega_ref_cbrt * omega_ref_cbrt;
                let (alpha, epsilon) =
                    nnlo.evaluate(omega_ref, logomega_ref, omega_ref_cbrt, omega_ref_cbrt2);
                (alpha - alpha0, epsilon - epsilon0)
            }
        };

        Ok(Self {
            version,
            precessing_tag,
            m1_si,
            m2_si,
            m1,
            m2,
            m1_2,
            m2_2,
            pi_gm,
            chi1,
            chi2,
            chi_p: chip,
            chi_tot_perp,
            sl,
            sperp,
            l,
            l_ref,
            j0,
            j0_norm,
            theta_j_sf,
            phi_j_sf,
            kappa,
            theta_jn,
            alpha0,
            epsilon0,
            zeta_polarization,
            nnlo,
            msa,
            alpha_offset,
            epsilon_offset,
            ylm: SphericalHarmonics::new(theta_jn),
            m_min,
            m_max,
            integration_buffer,
            f_ring_aligned,
            f_damp_aligned,
        })
    }

    /// NNLO Euler angles `(alpha, epsilon, cos(beta))` at the orbital
    /// frequency `omega = pi Mf` (`IMRPhenomX_PN_Euler_alpha_NNLO`,
    /// `IMRPhenomX_PN_Euler_epsilon_NNLO` and the opening angle of
    /// `IMRPhenomXPTwistUp22` / `IMRPhenomXWignerdCoefficients`), including
    /// the reference offsets.
    #[inline]
    pub(crate) fn nnlo_angles(&self, omega: f64, eta: f64) -> (f64, f64, f64) {
        let logomega = omega.ln();
        let omega_cbrt = omega.cbrt();
        let omega_cbrt2 = omega_cbrt * omega_cbrt;
        let (alpha, epsilon) = self.nnlo.evaluate(omega, logomega, omega_cbrt, omega_cbrt2);
        let alpha = alpha - self.alpha_offset;
        let epsilon = epsilon - self.epsilon_offset;
        let v = omega_cbrt;
        let l = self.l.ansatz_3pn(v, eta / v);
        // The sign of L + SL is ignored in s := Sp / (L + SL).
        let s = self.sperp / (l + self.sl);
        let s2 = s * s;
        let cos_beta = 1.0f64.copysign(l + self.sl) / (1.0 + s2).sqrt();
        (alpha, epsilon, cos_beta)
    }

    /// [`Self::nnlo_angles`] at `omega = scale f` for the first `alpha.len()`
    /// frequencies `f` of `grid`, written to `alpha`, `epsilon` and
    /// `cos_beta` (for the grid paths): the powers of `omega` are products
    /// of those stored in the grid, which leaves a loop without
    /// transcendental functions (agreeing with [`Self::nnlo_angles`] to
    /// rounding).
    pub(crate) fn nnlo_angles_on_grid(
        &self,
        eta: f64,
        grid: &FrequencyGrid,
        scale: f64,
        alpha: &mut [f64],
        epsilon: &mut [f64],
        cos_beta: &mut [f64],
    ) {
        let scale = ScalePowers::new(scale);
        let n = alpha.len();
        let (epsilon, cos_beta) = (&mut epsilon[..n], &mut cos_beta[..n]);
        let c = &self.nnlo;
        let l = &self.l;
        let sperp2 = self.sperp * self.sperp;
        for i in 0..n {
            let powers = grid.powers(&scale, i);
            let (inv, cbrt, log) = (powers.inv, powers.cbrt, powers.ln);
            // omega^(-1/3) and omega^(-2/3)
            let inv_cbrt = cbrt * cbrt * inv;
            let inv_cbrt2 = cbrt * inv;
            alpha[i] = fma(
                c.alpha5,
                cbrt,
                fma(
                    c.alpha4l,
                    log,
                    c.alpha1 * inv + c.alpha2 * inv_cbrt2 + c.alpha3 * inv_cbrt,
                ),
            ) - self.alpha_offset;
            epsilon[i] = fma(
                c.epsilon5,
                cbrt,
                fma(
                    c.epsilon4l,
                    log,
                    c.epsilon1 * inv + c.epsilon2 * inv_cbrt2 + c.epsilon3 * inv_cbrt,
                ),
            ) - self.epsilon_offset;
            // OrbitalAngularMomentum::ansatz_3pn at v = omega^(1/3)
            let v = cbrt;
            let x = v * v;
            let x2 = x * x;
            let mut sum = fma(l.l2, x, l.l0);
            sum = fma(l.l3, x * v, sum);
            sum = fma(l.l4, x2, sum);
            sum = fma(l.l5, x2 * v, sum);
            sum = fma(l.l6, x2 * x, sum);
            let l_total = fma(eta * inv_cbrt, sum, self.sl);
            // sign(L) / sqrt(1 + (Sp / L)^2)
            cos_beta[i] = l_total / fma(l_total, l_total, sperp2).sqrt();
        }
    }

    /// MSA Euler angles `(alpha, epsilon, cos(beta))` at the orbital
    /// velocity `v`, including the reference offsets: the MSA branch of
    /// `IMRPhenomXPTwistUp22` with `v = (pi Mf)^(1/3)`. IMRPhenomXPHM uses
    /// the same offsets for every co-precessing mode `m'`, evaluated at
    /// `v = (2 pi Mf / m')^(1/3)` (`Get_alpha_epsilon_offset`).
    ///
    /// Panics if the struct was not built for the MSA angles.
    #[inline]
    pub(crate) fn msa_angles(&self, v: f64) -> (f64, f64, f64) {
        let msa = self.msa.as_ref().expect("MSA system");
        let angles = msa.phi_zeta_costhetal(v);
        (
            angles.phiz - self.alpha_offset,
            angles.zeta - self.epsilon_offset,
            angles.cos_theta_l,
        )
    }

    /// [`Self::msa_angles`] at `v = x^(1/3)` for each of `x`, passed in
    /// order to `out(alpha, epsilon, cos(beta))` (for the grid paths).
    pub(crate) fn msa_angles_many(&self, x: &[f64], mut out: impl FnMut(f64, f64, f64)) {
        let msa = self.msa.as_ref().expect("MSA system");
        msa.phi_zeta_costhetal_many(x, |angles| {
            out(
                angles.phiz - self.alpha_offset,
                angles.zeta - self.epsilon_offset,
                angles.cos_theta_l,
            )
        });
    }
}

impl NnloCoefficients {
    /// The NNLO sums `alpha1 / omega + ... + alpha5 omega^(1/3)` (and the
    /// same for epsilon) before the offsets, with the products contracted
    /// into fused multiply-adds as in LAL's arm64 build.
    #[inline]
    fn evaluate(&self, omega: f64, logomega: f64, omega_cbrt: f64, omega_cbrt2: f64) -> (f64, f64) {
        let alpha = fma(
            self.alpha5,
            omega_cbrt,
            fma(
                self.alpha4l,
                logomega,
                self.alpha1 / omega + self.alpha2 / omega_cbrt2 + self.alpha3 / omega_cbrt,
            ),
        );
        let epsilon = fma(
            self.epsilon5,
            omega_cbrt,
            fma(
                self.epsilon4l,
                logomega,
                self.epsilon1 / omega + self.epsilon2 / omega_cbrt2 + self.epsilon3 / omega_cbrt,
            ),
        );
        (alpha, epsilon)
    }
}

/// NNLO Euler angle coefficients of `IMRPhenomXGetAndSetPrecessionVariables`
/// (versions 101-104), written exactly as in LAL.
#[allow(clippy::too_many_arguments)]
fn nnlo_coefficients(
    eta: f64,
    delta: f64,
    q: f64,
    chi_eff: f64,
    chip: f64,
    m1: f64,
    m1_2: f64,
    m1_3: f64,
    m1_4: f64,
    m1_5: f64,
    m1_6: f64,
    m1_8: f64,
) -> NnloCoefficients {
    let eta2 = eta * eta;
    let eta3 = eta * eta2;
    let eta4 = eta * eta3;
    let eta5 = eta * eta4;
    let eta6 = eta * eta5;
    let delta2 = delta * delta;
    let delta3 = delta * delta2;
    let chip2 = chip * chip;
    let chil = (1.0 + q) * (chi_eff / q);
    let chil2 = chil * chil;

    let alpha1 = -35. / 192. + (5. * delta) / (64. * m1);
    let alpha2 = ((15. * chil * delta * m1) / 128. - (35. * chil * m1_2) / 128.) / eta;
    // -5515/3072 + eta * (...) is contracted to a fused multiply-add
    let alpha3 = fma(
        eta,
        -515. / 384. - (15. * delta2) / (256. * m1_2) + (175. * delta) / (256. * m1),
        -5515. / 3072.,
    ) + (4555. * delta) / (7168. * m1)
        + ((15. * chip2 * delta * m1_3) / 128. - (35. * chip2 * m1_4) / 128.) / eta2;
    let alpha4l = (5. * chil * delta2) / 16. - (5. * chil * delta * m1) / 3.
        + (2545. * chil * m1_2) / 1152.
        + ((-2035. * chil * delta * m1) / 21504. + (2995. * chil * m1_2) / 9216.) / eta
        + ((5. * chil * chip2 * delta * m1_5) / 128. - (35. * chil * chip2 * m1_6) / 384.) / eta3
        - (35. * PI) / 48.
        + (5. * delta * PI) / (16. * m1);
    // alpha5 and epsilon5 with the fused multiply-adds of LAL's build.
    let b_in = fma(
        eta2,
        fma(784., eta, 323.),
        336. * fma(25., chil2, chip2) * m1_4,
    );
    let mut c_in = fma(8024297., eta4, 857412. * eta5);
    c_in = fma(3080448., eta6, c_in);
    c_in = fma(143640. * chip2 * eta2, m1_4, c_in);
    c_in = fnms(127008. * chip2 * fma(-4., chil2, chip2), m1_8, c_in);
    c_in = fma(
        6048. * eta3,
        fms(
            fma(2632., chil2, 115. * chip2),
            m1_4,
            672. * chil * m1_2 * PI,
        ),
        c_in,
    );
    let mut d_in = fma(-5579177., eta4, 80136. * eta5);
    d_in = fnms(3845520., eta6, d_in);
    d_in = fma(146664. * chip2 * eta2, m1_4, d_in);
    d_in = fma(127008. * chip2 * fma(-4., chil2, chip2), m1_8, d_in);
    d_in = fnms(
        42336. * eta3,
        fms(fma(726., chil2, 29. * chip2), m1_4, 96. * chil * m1_2 * PI),
        d_in,
    );
    let mut sum = fma(-190512. * delta3, eta6, 2268. * delta2 * eta3 * m1 * b_in);
    sum = fma(7. * m1_3, c_in, sum);
    sum = fma(3. * delta * m1_2, d_in, sum);
    let alpha5 = (5. * sum) / (6.5028096e7 * eta4 * m1_3);

    let epsilon1 = -35. / 192. + (5. * delta) / (64. * m1);
    let epsilon2 = ((15. * chil * delta * m1) / 128. - (35. * chil * m1_2) / 128.) / eta;
    let epsilon3 = fma(
        eta,
        -515. / 384. - (15. * delta2) / (256. * m1_2) + (175. * delta) / (256. * m1),
        -5515. / 3072.,
    ) + (4555. * delta) / (7168. * m1);
    let epsilon4l = (5. * chil * delta2) / 16. - (5. * chil * delta * m1) / 3.
        + (2545. * chil * m1_2) / 1152.
        + ((-2035. * chil * delta * m1) / 21504. + (2995. * chil * m1_2) / 9216.) / eta
        - (35. * PI) / 48.
        + (5. * delta * PI) / (16. * m1);
    let b_in = fma(eta2, fma(784., eta, 323.), 8400. * chil2 * m1_4);
    let c_in = fma(
        eta,
        fma(504. * eta, fma(7630., eta, -159.), 5579177.),
        254016. * chil * m1_2 * fma(121. * chil, m1_2, -16. * PI),
    );
    let d_in = fma(
        eta,
        fma(36. * eta, fma(85568., eta, 23817.), 8024297.),
        338688. * chil * m1_2 * fma(47. * chil, m1_2, -12. * PI),
    );
    let mut sum = fma(-190512. * delta3, eta3, 2268. * delta2 * m1 * b_in);
    sum = fnms(3. * delta * m1_2, c_in, sum);
    sum = fma(7. * m1_3, d_in, sum);
    let epsilon5 = (5. * sum) / (6.5028096e7 * eta * m1_3);

    NnloCoefficients {
        alpha1,
        alpha2,
        alpha3,
        alpha4l,
        alpha5,
        epsilon1,
        epsilon2,
        epsilon3,
        epsilon4l,
        epsilon5,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `SimIMRPhenomXPCalculateModelParametersFromSourceFrame` with
    /// `PhenomXPrecVersion = 102` from the lalsimulation 6.2.0 python bindings:
    /// (m1, m2, f_ref, phi_ref, inclination, chi1, chi2) and
    /// (chi1L, chi2L, chi_p, thetaJN, alpha0, phi_aligned, zeta_polarization).
    // The values are printed by lalsimulation: some (e.g. chi_p = 1/sqrt(2)
    // for the equal-mass case) happen to be close to mathematical constants.
    #[rustfmt::skip]
    #[allow(clippy::type_complexity, clippy::approx_constant)]
    const LAL_VALUES: [(f64, f64, f64, f64, f64, [f64; 3], [f64; 3], [f64; 7]); 4] = [
        (36.0, 29.0, 20.0, 1.2, 0.7, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1],
         [0.3, -0.1, 0.4472135954999579, 0.5936609957517057, 3.8490739691900155, -0.9704448051639335, -3.0067282752337894]),
        (8.0, 30.0, 35.0, 4.0, 2.3, [0.6, -0.3, 0.5], [0.1, 0.1, -0.7],
         [-0.7, 0.5, 0.15155571847498572, 2.621980790956642, 5.967133701779648, -0.5039188648694057, -3.005040273408551]),
        (30.0, 30.0, 25.0, 0.3, 1.6, [0.5, 0.5, 0.5], [0.5, -0.5, -0.5],
         [-0.5, 0.5, 0.7071067811865476, 1.527665180056583, 1.8691624719776798, -0.0, 2.9081862957742257]),
        (50.0, 5.0, 15.0, 5.5, 0.1, [0.0, 0.0, 0.6], [0.0, 0.0, -0.2],
         [0.6, -0.2, 0.0, 0.09999999999999945, 0.7876110196153103, -0.0, -3.1415926535897927]),
    ];

    #[test]
    fn frame_matches_lalsimulation() {
        for (m1, m2, f_ref, phi_ref, inclination, spin_1, spin_2, expected) in LAL_VALUES {
            let binary = PrecessingBinary {
                mass_1: m1,
                mass_2: m2,
                spin_1,
                spin_2,
                luminosity_distance: 1.0,
                inclination,
                phi_ref,
                f_ref,
            };
            let (m1, m2, chi1, chi2) = binary.ordered();
            // As set up by XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame
            let mut wf = WaveformStruct::new(
                m1,
                m2,
                chi1[2],
                chi2[2],
                0.125,
                f_ref,
                phi_ref,
                30.,
                1024.,
                1.0,
                inclination,
            )
            .unwrap();
            let prec = PrecessionStruct::new(
                &mut wf,
                m1 * MSUN_SI,
                m2 * MSUN_SI,
                chi1,
                chi2,
                PrecessionVersion::Nnlo,
                2,
                2,
            )
            .unwrap();
            let values = [
                prec.chi1[2],
                prec.chi2[2],
                prec.chi_p,
                prec.theta_jn,
                prec.alpha0,
                -prec.phi_j_sf,
                prec.zeta_polarization,
            ];
            for (k, (value, expected)) in values.into_iter().zip(expected).enumerate() {
                let mut diff = (value - expected).abs();
                if k >= 4 && !crate::reference::BIT_EXACT {
                    // Angles: for aligned spins zeta is +pi or -pi depending
                    // on the sign of a rounding error, which differs with
                    // glibc's libm.
                    diff = diff.min((diff - 2.0 * PI).abs());
                }
                assert!(
                    diff <= 4.0 * f64::EPSILON * expected.abs().max(1.0),
                    "{value} != {expected}"
                );
            }
        }
    }

    #[test]
    fn negligible_in_plane_spins_use_nnlo() {
        let mut wf =
            WaveformStruct::new(25., 20., 0.5, -0.3, 0., 20., 0.7, 20., 1000., 100., 0.9).unwrap();
        let prec = PrecessionStruct::new(
            &mut wf,
            25. * MSUN_SI,
            20. * MSUN_SI,
            [1e-8, 0.0, 0.5],
            [0.0, 5e-8, -0.3],
            PrecessionVersion::SpinTaylor,
            2,
            2,
        )
        .unwrap();
        assert_eq!(prec.version, 102);
        assert_eq!(prec.precessing_tag, 1);
    }

    #[test]
    fn rejects_spins_beyond_the_kerr_bound() {
        let mut wf =
            WaveformStruct::new(25., 20., 0.5, -0.3, 0., 20., 0.7, 20., 1000., 100., 0.9).unwrap();
        let error = PrecessionStruct::new(
            &mut wf,
            25. * MSUN_SI,
            20. * MSUN_SI,
            [0.9, 0.0, 0.5],
            [0.0, 0.0, -0.3],
            PrecessionVersion::SpinTaylor,
            2,
            2,
        )
        .unwrap_err();
        assert!(matches!(error, IMRPhenomXError::UnphysicalSpinMagnitude(_)));
    }
}
