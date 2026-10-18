//! Multiple-scale-analysis (MSA) Euler angles of IMRPhenomXP/XPHM
//! (`PhenomXPrecVersion` 223; Chatziioannou et al. 2017, arXiv:1703.03967,
//! as implemented in `LALSimInspiralFDPrecAngles`), following the MSA
//! routines of `LALSimIMRPhenomX_precession.c`:
//!
//! - [`MsaSystem::new`] ports `IMRPhenomX_Initialize_MSA_System`: the
//!   precession-averaged spin couplings, the PN coefficients of
//!   `d omega / dt`, the phase of the spin magnitude oscillations and the
//!   coefficients of `phi_z` and `zeta`, including their values at the
//!   reference frequency.
//! - [`MsaSystem::phi_zeta_costhetal`] ports
//!   `IMRPhenomX_Return_phi_zeta_costhetaL_MSA`: `phi_z` (= `alpha`), `zeta`
//!   (= `-gamma`, before the reference offsets) and `cos(theta_L)` (=
//!   `cos(beta)`) at the orbital velocity `v`.
//! - [`MsaSystem::final_spin`] is the precession-averaged final spin of
//!   `IMRPhenomX_SetPrecessingRemnantParams` (`PhenomXPFinalSpinMod` 3, the
//!   version used when the default 4 is requested with MSA angles).
//!
//! The Jacobi elliptic function `sn` and the incomplete elliptic integral of
//! the first kind are ports of GSL's `gsl_sf_elljac_e` and
//! `gsl_sf_ellint_F` (GSL 2.8).
//!
//! # Floating point
//!
//! Like the rest of the IMRPhenomXP port, the expressions mirror the fused
//! multiply-adds that clang emits for the arm64 builds of LALSimulation and
//! GSL (see [`crate::spintaylor::fp`]): for each `+`/`-` whose left operand is
//! a (possibly negated) product within the same C expression, the product is
//! fused with the other operand; otherwise a product on the right is fused.
//!
//! # State
//!
//! LAL's `IMRPhenomX_Return_phi_zeta_costhetaL_MSA` stores the roots
//! `S_3^2`, `S_-^2`, `S_+^2` and the spin magnitude at each frequency in the
//! shared `IMRPhenomXPrecessionStruct`. As it is called from inside the
//! `#pragma omp parallel for` of `IMRPhenomXPGenerateFD` (the IMRPhenomXPHM
//! loops are not parallelised), multi-threaded LAL builds race on these members and give
//! run-to-run differences of up to ~1e-2 in the MSA waveforms (this is also
//! what made the SpinTaylor -> MSA fallback appear to differ from direct
//! `PhenomXPrecVersion` 223 calls). Here the per-frequency roots are local
//! values and [`MsaSystem`] is immutable after construction; the results
//! agree with single-threaded LAL (`OMP_NUM_THREADS=1`).

use std::f64::consts::{PI, TAU};

use crate::cmath::LibmCbrt;
use crate::fastmath;
use crate::imrphenomx::precession::OrbitalAngularMomentum;
use crate::spintaylor::fp::{fma, fms, fnms};

/// `GSL_DBL_EPSILON`.
const GSL_DBL_EPSILON: f64 = f64::EPSILON;
/// `GSL_SQRT_DBL_EPSILON`.
const GSL_SQRT_DBL_EPSILON: f64 = 1.0 / 67108864.0; // 2^-26
/// `GSL_DBL_MIN` and `GSL_DBL_MAX`.
const GSL_DBL_MIN: f64 = f64::MIN_POSITIVE;
const GSL_DBL_MAX: f64 = f64::MAX;

/// `IMRPhenomX_vector_dot_product` as compiled (with fused multiply-adds).
#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    fma(a[2], b[2], fma(a[0], b[0], a[1] * b[1]))
}

/// `IMRPhenomX_vector_cross_product`.
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        fms(a[1], b[2], a[2] * b[1]),
        fms(a[2], b[0], a[0] * b[2]),
        fms(a[0], b[1], a[1] * b[0]),
    ]
}

/// `IMRPhenomX_vector_L2_norm`.
#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// The roots `S_3^2 <= S_-^2 <= S_+^2` of the cubic for the total spin
/// magnitude, Eq. 21 of arXiv:1703.03967 (`IMRPhenomX_Return_Roots_MSA`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Roots {
    pub(crate) s32: f64,
    pub(crate) smi2: f64,
    pub(crate) spl2: f64,
}

impl Roots {
    #[inline]
    fn spl2_m_smi2(&self) -> f64 {
        self.spl2 - self.smi2
    }

    /// Whether `|S_-^2 - S_+^2| < 1e-5`, the tolerance of
    /// `LALSimInspiralFDPrecAngles` below which the spin magnitude is taken
    /// as constant (false for NaN, as the comparisons in LAL).
    #[inline]
    fn constant_magnitude(&self) -> bool {
        (self.smi2 - self.spl2).abs() < 1.0e-5
    }

    /// Whether the MSA corrections are evaluated, `|S_-^2 - S_+^2| > 1e-5`.
    #[inline]
    fn has_corrections(&self) -> bool {
        (self.smi2 - self.spl2).abs() > 1.0e-5
    }
}

/// Powers of the orbital velocity `v` shared by the MSA functions, each
/// computed by the same operations as in LAL.
#[derive(Clone, Copy, Debug, Default)]
struct VPowers {
    v: f64,
    /// `v * v`
    v2: f64,
    /// `v * v2`
    v3: f64,
    /// `v2 * v2`
    v4: f64,
    /// `v3 * v3`
    v6: f64,
    /// `1 / v`
    invv: f64,
    /// `invv * invv`
    invv2: f64,
    /// `invv * invv2`
    invv3: f64,
    /// `log(v)`
    logv: f64,
}

impl VPowers {
    #[inline]
    fn new(v: f64) -> Self {
        let v2 = v * v;
        let v3 = v * v2;
        let v4 = v2 * v2;
        let invv = 1.0 / v;
        let invv2 = invv * invv;
        Self {
            v,
            v2,
            v3,
            v4,
            v6: v3 * v3,
            invv,
            invv2,
            invv3: invv * invv2,
            logv: v.ln(),
        }
    }
}

/// `phi_z`, `zeta` and `cos(theta_L)` at one frequency
/// (`IMRPhenomX_Return_phi_zeta_costhetaL_MSA`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MsaAngles {
    /// Precession angle of L about J, the Euler angle `alpha` before the
    /// reference offset.
    pub(crate) phiz: f64,
    /// `zeta = -gamma`, the Euler angle `epsilon` before the reference
    /// offset.
    pub(crate) zeta: f64,
    /// Cosine of the angle between L and J, `cos(beta)`.
    pub(crate) cos_theta_l: f64,
}

/// The MSA system of one binary, the MSA members of
/// `IMRPhenomXPrecessionStruct` after `IMRPhenomX_Initialize_MSA_System`
/// for `PhenomXPrecVersion` 223 and `ExpansionOrder` 5.
#[derive(Clone, Debug)]
pub(crate) struct MsaSystem {
    eta: f64,
    eta2: f64,
    eta3: f64,
    eta4: f64,
    inveta: f64,
    inveta2: f64,
    inveta3: f64,
    inveta4: f64,
    sqrt_inveta: f64,
    /// Mass ratio `m2 / m1 <= 1` and `(1 - q) / (1 + q)`.
    qq: f64,
    delta_qq: f64,
    /// Squared spin norms (in units of `M^2`), individual and total at the
    /// reference frequency.
    s1_norm_2: f64,
    s2_norm_2: f64,
    s_0_norm_2: f64,
    /// Whether the roots of the spin magnitude cubic take their degenerate
    /// values for all frequencies: a spin (anti-)aligned with L or
    /// vanishing (`dotS1Ln == +-1`, ..., `S1_norm_2 == 0`, ... in
    /// `IMRPhenomX_Return_Roots_MSA`).
    degenerate_roots: bool,
    /// Effective spin `xi`.
    seff: f64,
    /// Coefficients of the 3PN orbital angular momentum.
    constants_l: [f64; 5],
    c1: f64,
    c1_over_eta: f64,
    sav2: f64,
    sav: f64,
    inv_sav: f64,
    inv_sav2: f64,
    /// Precession-averaged `S1.L` and `S2.L` (for the final spin).
    s1l_pav: f64,
    s2l_pav: f64,
    /// Phase of the spin-magnitude oscillations, Eq. 51.
    psi0: f64,
    psi1: f64,
    psi2: f64,
    g0: f64,
    /// `<Omega_z>^(n)` and `<Omega_zeta>^(n)` coefficients (Eqs. D16-D21,
    /// F6-F11), `n = 5` set to zero by the expansion order.
    omegaz_coeff: [f64; 6],
    omegazeta_coeff: [f64; 6],
    /// Integration constants of `phi_z` and `zeta`.
    phiz_0: f64,
    zeta_0: f64,
    /// `MSA_ERROR`: `|Omega_z^(5)| > 1000` (the system is complete anyway).
    msa_error: bool,
}

/// Inputs of `IMRPhenomX_Initialize_MSA_System` taken from the waveform
/// and precession structs.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MsaInputs {
    /// Symmetric mass ratio.
    pub(crate) eta: f64,
    /// `pWF->m1`, `pWF->m2`: component masses of the aligned-spin struct in
    /// units of the total mass.
    pub(crate) m1: f64,
    pub(crate) m2: f64,
    /// `m1_SI / Mtot_SI` squared and `m2_SI / Mtot_SI` squared (the
    /// normalisation of the spin norms in the precession struct).
    pub(crate) m1_2: f64,
    pub(crate) m2_2: f64,
    /// Dimensionless spins in the L frame at `f_ref`.
    pub(crate) chi1: [f64; 3],
    pub(crate) chi2: [f64; 3],
    /// `pi G M / c^3` (s) and the reference frequency (Hz).
    pub(crate) pi_gm: f64,
    pub(crate) f_ref: f64,
}

impl MsaSystem {
    /// `IMRPhenomX_Initialize_MSA_System` for `PhenomXPrecVersion` 223.
    ///
    /// Returns `None` when LAL sets `MSA_ERROR` (`|Omega_z^(5)| > 1000`),
    /// in which case LAL falls back to the NNLO angles (version 102).
    pub(crate) fn new(input: &MsaInputs) -> Option<Self> {
        let system = Self::new_unchecked(input);
        (!system.msa_error).then_some(system)
    }

    /// `IMRPhenomX_Initialize_MSA_System` without the check of `MSA_ERROR`:
    /// like LAL, the system is completed even when it sets `MSA_ERROR`
    /// (which IMRPhenomXPNR ignores in some of its uses).
    pub(crate) fn new_unchecked(input: &MsaInputs) -> Self {
        let eta = input.eta;
        let eta2 = eta * eta;
        let eta3 = eta * eta2;
        let eta4 = eta * eta3;
        let inveta = 1.0 / eta;
        let inveta2 = 1.0 / eta2;
        let inveta3 = 1.0 / eta3;
        let inveta4 = 1.0 / eta4;
        let sqrt_inveta = 1.0 / eta.sqrt();
        let (m1, m2) = (input.m1, input.m2);

        // Spin norms of the precession struct (pPrec->S1_norm_2, ...).
        let chi1_norm = dot(input.chi1, input.chi1).sqrt();
        let chi2_norm = dot(input.chi2, input.chi2).sqrt();
        let s1_norm = chi1_norm.abs() * input.m1_2;
        let s2_norm = chi2_norm.abs() * input.m2_2;
        let s1_norm_2 = s1_norm * s1_norm;
        let s2_norm_2 = s2_norm * s2_norm;

        // PN coefficients of d omega / dt (LALSimInspiralFDPrecAngles).
        let domegadt_ns: [f64; 9] = [
            96. / 5.,
            -1486. / 35.,
            -264. / 5.,
            384. * PI / 5.,
            34103. / 945.,
            13661. / 105.,
            944. / 15.,
            PI * (-4159. / 35.),
            PI * (-2268. / 5.),
        ];
        let domegadt_so: [f64; 6] = [
            -904. / 5.,
            -120.,
            -62638. / 105.,
            4636. / 5.,
            -6472. / 35.,
            3372. / 5.,
        ];
        let domegadt_ss: [f64; 4] = [-494. / 5., -1442. / 5., -233. / 5., -719. / 5.];
        let l_nonspin: [f64; 9] = [
            3. / 2.,
            1. / 6.,
            27. / 8.,
            -19. / 8.,
            1. / 24.,
            135. / 16.,
            -6889. / 144. + 41. / 24. * PI * PI,
            31. / 24.,
            7. / 1296.,
        ];
        let l_spinorbit: [f64; 6] = [
            -14. / 6.,
            -3. / 2.,
            -11. / 2.,
            133. / 72.,
            -33. / 8.,
            7. / 4.,
        ];

        // Chatziioannou et al use q = m2 / m1 < 1.
        let q = m2 / m1;
        let qq = q;
        let delta_qq = (1.0 - qq) / (1.0 + qq);

        // Dimensionful spins: eta / q = m1^2, eta q = m2^2.
        let s1v = input.chi1.map(|c| c * eta / q);
        let s2v = input.chi2.map(|c| c * eta * q);
        let s1_0_norm = norm(s1v);
        let s2_0_norm = norm(s2v);

        let v_0 = (input.pi_gm * input.f_ref).libm_cbrt();
        let v_0_2 = v_0 * v_0;
        let lhat = [0.0, 0.0, 1.0];
        let l_0 = lhat.map(|c| c * (eta / v_0));

        let dot_s1l = dot(s1v, lhat);
        let dot_s2l = dot(s2v, lhat);
        let dot_s1s2 = dot(s1v, s2v);
        let dot_s1_ln = dot_s1l / s1_0_norm;
        let dot_s2_ln = dot_s2l / s2_0_norm;

        // IMRPhenomX_Get_PN_beta, _sigma, _tau
        let pn_beta = |a: f64, b: f64| fma(dot_s1l, fma(b, qq, a), dot_s2l * (a + b / qq));
        let pn_sigma = |a: f64, b: f64| inveta * fms(a, dot_s1s2, b * dot_s1l * dot_s2l);
        let pn_tau = |a: f64, b: f64| {
            fma(
                qq,
                fms(s1_norm_2, a, b * dot_s1l * dot_s1l),
                fms(a, s2_norm_2, b * dot_s2l * dot_s2l) / qq,
            ) / eta
        };

        let constants_l = [
            fma(eta, l_nonspin[1], l_nonspin[0]),
            pn_beta(l_spinorbit[0], l_spinorbit[1]),
            fma(
                eta * eta,
                l_nonspin[4],
                fma(eta, l_nonspin[3], l_nonspin[2]),
            ),
            pn_beta(
                fma(l_spinorbit[3], eta, l_spinorbit[2]),
                fma(l_spinorbit[5], eta, l_spinorbit[4]),
            ),
            fma(
                l_nonspin[8] * eta * eta,
                eta,
                fma(
                    l_nonspin[7] * eta,
                    eta,
                    fma(l_nonspin[6], eta, l_nonspin[5]),
                ),
            ),
        ];

        // Effective spin
        let seff = fma(1.0 + q, dot_s1l, (1. + (1.0 / q)) * dot_s2l);

        // Total spin and angular momentum at the reference frequency.
        let s_0: [f64; 3] = std::array::from_fn(|k| s1v[k] + s2v[k]);
        let j_0: [f64; 3] = std::array::from_fn(|k| l_0[k] + s_0[k]);
        let s_0_norm = norm(s_0);
        let s_0_norm_2 = s_0_norm * s_0_norm;
        let l_0_norm = norm(l_0);
        let j_0_norm = norm(j_0);

        let mut system = Self {
            eta,
            eta2,
            eta3,
            eta4,
            inveta,
            inveta2,
            inveta3,
            inveta4,
            sqrt_inveta,
            qq,
            delta_qq,
            s1_norm_2,
            s2_norm_2,
            s_0_norm_2,
            degenerate_roots: dot_s1_ln == 1.
                || dot_s2_ln == 1.
                || dot_s1_ln == -1.
                || dot_s2_ln == -1.
                || s1_norm_2 == 0.
                || s2_norm_2 == 0.,
            seff,
            constants_l,
            c1: 0.0,
            c1_over_eta: 0.0,
            sav2: 0.0,
            sav: 0.0,
            inv_sav: 0.0,
            inv_sav2: 0.0,
            s1l_pav: 0.0,
            s2l_pav: 0.0,
            psi0: 0.0,
            psi1: 0.0,
            psi2: 0.0,
            g0: 0.0,
            omegaz_coeff: [0.0; 6],
            omegazeta_coeff: [0.0; 6],
            phiz_0: 0.0,
            zeta_0: 0.0,
            msa_error: false,
        };

        let roots = system.roots(l_0_norm, j_0_norm);

        // Eq. 45, from the initial conditions
        let sav2 = 0.5 * (roots.spl2 + roots.smi2);
        let sav = sav2.sqrt();
        system.sav2 = sav2;
        system.sav = sav;
        system.inv_sav2 = 1.0 / sav2;
        system.inv_sav = 1.0 / sav;

        // c_1, Eq. 41
        let c_1 = 0.5 * (fms(j_0_norm, j_0_norm, l_0_norm * l_0_norm) - sav2) / l_0_norm * eta;
        let c1_2 = c_1 * c_1;
        system.c1 = c_1;
        system.c1_over_eta = c_1 / eta;

        // Precession averaged spin couplings, Eqs. A9 - A14
        let omq2 = fnms(q, q, 1.0) + 1e-16;
        let s1l_pav = fms(c_1, 1.0 + q, q * eta * seff) / (eta * omq2);
        let s2l_pav = -q * fms(c_1, 1.0 + q, eta * seff) / (eta * omq2);
        system.s1l_pav = s1l_pav;
        system.s2l_pav = s2l_pav;
        // The other averaged couplings (S1.S2, (S1.L)^2, ...) only enter
        // the PN coefficients a3, a4, a5 of versions 220/221/224 (through
        // beta3, beta5, sigma4), which version 223 replaces below.

        // Version 223: PN coefficients from the initial spin couplings, as
        // in LALSimInspiralFDPrecAngles.
        let a0 = eta * domegadt_ns[0];
        let a2 = eta * fma(eta, domegadt_ns[2], domegadt_ns[1]);
        let a3 = eta * (domegadt_ns[3] + pn_beta(domegadt_so[0], domegadt_so[1]));
        let a4 = eta
            * (fma(
                eta,
                fma(eta, domegadt_ns[6], domegadt_ns[5]),
                domegadt_ns[4],
            ) + pn_sigma(domegadt_ss[0], domegadt_ss[1])
                + pn_tau(domegadt_ss[2], domegadt_ss[3]));

        let a0_2 = a0 * a0;
        let a0_3 = a0_2 * a0;
        let a2_2 = a2 * a2;

        // Appendix A of arXiv:1703.03967
        let g0 = 1. / a0;
        let g2 = -(a2 / a0_2);
        let g3 = -(a3 / a0_2);
        let g4 = -fms(a4, a0, a2_2) / a0_3;
        // a5 and g5 only enter the n = 5 coefficients, which ExpansionOrder
        // 5 drops.
        system.g0 = g0;

        let delta = delta_qq;
        let delta2 = delta * delta;

        // Eq. C1
        system.psi1 = 3.0 * fms(2.0 * eta2, seff, c_1) / (eta * delta2);

        let c_1_over_nu = system.c1_over_eta;
        let c_1_over_nu_2 = c_1_over_nu * c_1_over_nu;
        let one_p_q_sq = (1. + q) * (1. + q);
        let seff_2 = seff * seff;
        let q_2 = q * q;
        let one_m_q_sq = (1. - q) * (1. - q);
        let one_m_q2_2 = (1. - q_2) * (1. - q_2);
        let one_m_q_4 = one_m_q_sq * one_m_q_sq;

        // Delta as in LALSimInspiralFDPrecAngles
        let del1 = 4. * c_1_over_nu_2 * one_p_q_sq;
        let del2 = 8. * c_1_over_nu * q * (1. + q) * seff;
        let del3 = 4. * fms(one_m_q2_2, s1_norm_2, q_2 * seff_2);
        let del4 = 4. * c_1_over_nu_2 * q_2 * one_p_q_sq;
        let del5 = 8. * c_1_over_nu * q_2 * (1. + q) * seff;
        let del6 = 4. * fms(one_m_q2_2, s2_norm_2, q_2 * seff_2);
        let big_delta = ((del1 - del2 - del3) * (del4 - del5 - del6)).abs().sqrt();

        // Eq. C2 as in LALSimInspiralFDPrecAngles
        let u1 = 3. * g2 / g0;
        let u2 = 0.75 * one_p_q_sq / one_m_q_4;
        let u3 = -20. * c_1_over_nu_2 * q_2 * one_p_q_sq;
        let u4 = 2.
            * one_m_q2_2
            * fnms(
                2. * q,
                sav2,
                fma(q * (2. + q), s1_norm_2, fma(2., q, 1.) * s2_norm_2),
            );
        let u5 = 2. * q_2 * fma(7., q_2, fma(6., q, 7.)) * 2. * c_1_over_nu * seff;
        let u6 = 2. * q_2 * fma(3., q_2, fma(4., q, 3.)) * seff_2;
        let u7 = q * big_delta;
        system.psi2 = fma(u2, u3 + u4 + u5 - u6 + u7, u1);

        // Appendix D
        let rm = roots.spl2 - roots.smi2;
        let rm_2 = rm * rm;
        let cp = fms(roots.spl2, eta2, c1_2);
        let cm = fms(roots.smi2, eta2, c1_2);
        let cpcm = (cp * cm).abs();
        let sqrt_cpcm = cpcm.sqrt();
        let a1dd = 0.5 + 0.75 / eta;
        let a2dd = -0.75 * seff / eta;
        let d2rmsq = (cp - sqrt_cpcm) / eta2;
        let d4rmsq = fnms(cp / eta4, sqrt_cpcm - cp, -0.5 * rm * sqrt_cpcm / eta2);
        let s0m = s1_norm_2 - s2_norm_2;
        let aw = -3. * (1. + q) / q
            * fma(
                (1. - q) * eta2,
                s0m,
                fms(2. * (1. + q) * eta2 * seff, c_1, (1. + q) * c1_2),
            );
        let cw = 3. / 32. / eta * rm_2;
        let dw = fms(4.0, cp, 4.0 * d2rmsq * eta2);
        let hw = -2.0 * fms(2.0, d2rmsq, rm) * c_1;
        let fw = fnms(0.25, rm_2, fms(rm, d2rmsq, d4rmsq));
        let add = aw / dw;
        let hdd = hw / dw;
        let cdd = cw / dw;
        let fdd = fw / dw;
        let gw = 3. / 16. / eta2 / eta * rm_2 * fnms(eta2, seff, c_1);
        let gdd = gw / dw;
        let hdd_2 = hdd * hdd;
        let addfdd = add * fdd;
        let addfddhdd = addfdd * hdd;
        let addhdd_2 = add * hdd_2;

        // Eqs. D10 - D15
        let omegaz0 = a1dd + add;
        let omegaz1 = fnms(add, hdd, fnms(add, seff, a2dd));
        let omegaz2 = fma(add, hdd_2, fnms(add, fdd, fma(add * hdd, seff, cdd)));
        let omegaz3 = fma(addfdd - cdd - addhdd_2, seff + hdd, addfddhdd);
        let omegaz4 = fms(
            fnms(2.0, addfdd, cdd + addhdd_2),
            fma(hdd, seff, hdd_2) - fdd,
            add * fdd * fdd,
        );
        let omegaz5 = fnms(
            addfdd * fdd,
            hdd,
            fms(
                (cdd - addfdd + addhdd_2) * fdd,
                fma(2.0, hdd, seff),
                fnms(2.0, addfdd, cdd + addhdd_2) * hdd_2 * (seff + hdd),
            ),
        );
        // MSA_ERROR: "|Omegaz5| is larger than expected and may be
        // pathological. Triggering MSA failure."
        system.msa_error = omegaz5.abs() > 1000.0;

        // Eqs. D16 - D21; ExpansionOrder 5 drops the n = 5 terms.
        system.omegaz_coeff = [
            3.0 * g0 * omegaz0,
            3.0 * g0 * omegaz1,
            3.0 * fma(g0, omegaz2, g2 * omegaz0),
            3.0 * fma(g3, omegaz0, fma(g0, omegaz3, g2 * omegaz1)),
            3.0 * fma(
                g4,
                omegaz0,
                fma(g3, omegaz1, fma(g0, omegaz4, g2 * omegaz2)),
            ),
            0.0,
        ];

        // Appendix E/F
        let c1oveta2 = c_1 / eta2;
        let omegazeta0 = omegaz0;
        let omegazeta1 = fma(omegaz0, c1oveta2, omegaz1);
        let omegazeta2 = fma(omegaz1, c1oveta2, omegaz2);
        let omegazeta3 = fma(omegaz2, c1oveta2, omegaz3) + gdd;
        let omegazeta4 = fnms(gdd, hdd, fnms(gdd, seff, fma(omegaz3, c1oveta2, omegaz4)));
        system.omegazeta_coeff = [
            -g0 * omegazeta0,
            -1.5 * g0 * omegazeta1,
            -3.0 * fma(g0, omegazeta2, g2 * omegazeta0),
            3.0 * fma(g3, omegazeta0, fma(g0, omegazeta3, g2 * omegazeta1)),
            3.0 * fma(
                g4,
                omegazeta0,
                fma(g3, omegazeta1, fma(g0, omegazeta4, g2 * omegazeta2)),
            ),
            0.0,
        ];

        // psi0
        if !roots.constant_magnitude() {
            let mm = ((roots.smi2 - roots.spl2) / (roots.s32 - roots.spl2)).sqrt();
            let tmp_b = fms(s_0_norm, s_0_norm, roots.spl2) / (roots.smi2 - roots.spl2);
            let volume_element = dot(cross(l_0, s1v), s2v);
            let vol_sign = if volume_element > 0.0 {
                1.0
            } else if volume_element < 0.0 {
                -1.0
            } else {
                0.0
            };
            let psi_of_v0 = system.psi_of_v(v_0, v_0_2, 0.0);
            if !(0.0..=1.0).contains(&tmp_b) && !tmp_b.is_nan() {
                if tmp_b > 1.0 && (tmp_b - 1.) < 0.00001 {
                    system.psi0 = ellint_f((vol_sign * 1.0f64.sqrt()).asin(), mm) - psi_of_v0;
                }
                if tmp_b < 0.0 && tmp_b > -0.00001 {
                    system.psi0 = ellint_f((vol_sign * 0.0f64.sqrt()).asin(), mm) - psi_of_v0;
                }
            } else {
                system.psi0 = ellint_f((vol_sign * tmp_b.sqrt()).asin(), mm) - psi_of_v0;
            }
        }

        // Integration constants phi_z,0 and zeta_0
        let v0 = VPowers::new(v_0);
        let corrections = if roots.has_corrections() {
            system.msa_corrections(&v0, l_0_norm, j_0_norm, &roots)
        } else {
            [0.0, 0.0]
        };
        let phiz_0 = system.phiz(&v0, j_0_norm);
        let zeta_0 = system.zeta(&v0);
        system.phiz_0 = -phiz_0 - corrections[0];
        system.zeta_0 = -zeta_0 - corrections[1];
        system
    }

    /// `cos(theta_L)` of [`Self::phi_zeta_costhetal`] as LAL evaluates it
    /// when the precession struct does not carry version 222 or 223 (as in
    /// the PNR angles): the 3PN orbital angular momentum is
    /// `XLALSimIMRPhenomXLPNAnsatz` with the coefficients `l` instead of the
    /// MSA expansion.
    pub(crate) fn cos_theta_l_with(&self, v: f64, l: &OrbitalAngularMomentum) -> f64 {
        let powers = VPowers::new(v);
        let l_norm = self.eta / v;
        let j_norm = self.j_norm(l_norm);
        let l_norm_3pn = l.ansatz(v, l_norm);
        let j_norm_3pn = self.j_norm(l_norm_3pn);
        let roots = self.roots(l_norm, j_norm);
        let s_norm = self.s_norm(&powers, &roots);
        cos_theta_lj(l_norm_3pn, j_norm_3pn, s_norm)
    }

    /// `MSA_ERROR` of `IMRPhenomX_Initialize_MSA_System`.
    pub(crate) fn msa_error(&self) -> bool {
        self.msa_error
    }

    /// `g0`, `delta_qq`, `psi1` and `psi2` of the precession struct (used by
    /// the PNR angles to estimate the period of the two-spin oscillations).
    pub(crate) fn psi_terms(&self) -> (f64, f64, f64, f64) {
        (self.g0, self.delta_qq, self.psi1, self.psi2)
    }

    /// `IMRPhenomX_psiofv`
    #[inline]
    fn psi_of_v(&self, v: f64, v2: f64, psi0: f64) -> f64 {
        psi0 - 0.75 * self.g0 * self.delta_qq * fma(self.psi2, v2, fma(self.psi1, v, 1.0))
            / (v2 * v)
    }

    /// `IMRPhenomX_JNorm_MSA`, Eq. 41
    #[inline]
    fn j_norm(&self, l_norm: f64) -> f64 {
        (fma(l_norm, l_norm, 2.0 * l_norm * self.c1_over_eta) + self.sav2).sqrt()
    }

    /// `IMRPhenomX_L_norm_3PN_of_v`
    #[inline]
    fn l_norm_3pn(&self, v: f64, v2: f64, l_norm: f64) -> f64 {
        let c = &self.constants_l;
        let inner = fma(v2, c[4], fma(v, c[3], c[2]));
        l_norm * fma(v2, fma(v2, inner, fma(v, c[1], c[0])), 1.)
    }

    /// `IMRPhenomX_Return_Roots_MSA` with
    /// `IMRPhenomX_Return_Spin_Evolution_Coefficients_MSA`.
    fn roots(&self, l_norm: f64, j_norm: f64) -> Roots {
        let degenerate = || {
            let smi2 = self.s_0_norm_2;
            // numerical perturbation to keep alpha from diverging
            Roots {
                s32: 0.0,
                smi2,
                spl2: smi2 + 1e-9,
            }
        };
        if self.degenerate_roots {
            return degenerate();
        }
        let j_norm2 = j_norm * j_norm;
        let l_norm2 = l_norm * l_norm;
        let s1 = self.s1_norm_2;
        let s2 = self.s2_norm_2;
        let q = self.qq;
        let eta = self.eta;
        let j2ml2 = j_norm2 - l_norm2;
        let j2ml2sq = j2ml2 * j2ml2;
        let delta = self.delta_qq;
        let delta_sq = delta * delta;
        let seff = self.seff;

        // Eqs. B2 - B4
        let b = fnms(2.0, j_norm2, fma(l_norm2 + s1, q, 2.0 * l_norm * seff)) - s1 - s2
            + (l_norm2 + s2) / q;
        let one_m_q_over_q = (1.0 - q) / q;
        let mut c = fnms(2.0 * l_norm * seff, j2ml2, j2ml2sq);
        c = fnms(2.0 * one_m_q_over_q * l_norm2, fnms(q, s2, s1), c);
        c = fma(4.0 * eta * l_norm2 * seff, seff, c);
        c = fnms(2.0 * delta * (s1 - s2) * seff, l_norm, c);
        c = fma(2.0 * one_m_q_over_q * fms(q, s1, s2), j_norm2, c);
        let d = fma(
            2.0 * delta * l_norm * seff * (s1 - s2),
            j2ml2,
            fma(
                one_m_q_over_q * fnms(q, s1, s2),
                j2ml2sq,
                delta_sq * (s1 - s2) * (s1 - s2) * l_norm2 / eta,
            ),
        );

        let b2 = b * b;
        let b3 = b2 * b;
        let bc = b * c;
        let p = c - b2 / 3.;
        let qc = fms(2.0 / 27.0, b3, bc / 3.0) + d;
        let sqrtarg = (-p / 3.0).sqrt();
        let acosarg = (1.5 * qc / p / sqrtarg).clamp(-1.0, 1.0);
        // clamp keeps NaN, as the comparisons in LAL do
        let theta = acosarg.acos() / 3.0;

        if theta.is_nan() || sqrtarg.is_nan() {
            degenerate()
        } else {
            let cos_theta = theta.cos();
            let tmp1 = fms(2.0 * sqrtarg, (theta - 2.0 * TAU / 3.0).cos(), b / 3.0);
            let tmp2 = fms(2.0 * sqrtarg, (theta - TAU / 3.0).cos(), b / 3.0);
            let tmp3 = fms(2.0 * sqrtarg, cos_theta, b / 3.0);
            let tmp4 = tmp1.max(tmp2).max(tmp3);
            let tmp5 = tmp1.min(tmp2).min(tmp3);
            let tmp6 = if (tmp4 - tmp3) > 0.0 && (tmp5 - tmp3) < 0.0 {
                tmp3
            } else if (tmp4 - tmp1) > 0.0 && (tmp5 - tmp1) < 0.0 {
                tmp1
            } else {
                tmp2
            };
            Roots {
                s32: tmp5,
                smi2: tmp6.abs(),
                spl2: tmp4.abs(),
            }
        }
    }

    /// `IMRPhenomX_Return_SNorm_MSA`, Eq. 23
    #[inline]
    fn s_norm(&self, v: &VPowers, roots: &Roots) -> f64 {
        self.s_norm_with(v, roots, elljac_sn)
    }

    /// [`Self::s_norm`] with the Jacobi elliptic function `sn(u | m)`.
    #[inline(always)]
    fn s_norm_with(&self, v: &VPowers, roots: &Roots, elljac_sn: impl Fn(f64, f64) -> f64) -> f64 {
        let sn = if roots.constant_magnitude() {
            0.0
        } else {
            // Eq. 25
            let m = (roots.smi2 - roots.spl2) / (roots.s32 - roots.spl2);
            let psi = self.psi_of_v(v.v, v.v2, self.psi0);
            elljac_sn(psi, m)
        };
        fma((roots.smi2 - roots.spl2) * sn, sn, roots.spl2).sqrt()
    }

    /// `IMRPhenomX_Return_phiz_MSA`, Eq. 66 (without the integration
    /// constant `phi_z,0`, which is added by the caller; the NaN check is
    /// done after adding it, as in LAL).
    #[inline]
    fn phiz_sum(&self, v: &VPowers, j_norm: f64) -> f64 {
        let eta = self.eta;
        let invv = v.invv;
        let invv2 = invv * invv;
        let l_newt = eta / v.v;
        let c1 = self.c1;
        let c12 = c1 * c1;
        let sav2 = self.sav2;
        let sav = self.sav;
        let inv_sav = self.inv_sav;
        let inv_sav2 = self.inv_sav2;
        let vv = v.v;

        // Eqs. D27, D28
        let log1 = fma(eta, l_newt, fma(j_norm, eta, c1)).abs().ln();
        let log2 = fma(sav2, vv, fma(j_norm * sav, vv, c1)).abs().ln();

        // Eqs. D22 - D27
        let phiz_0_coeff = fms(
            j_norm * self.inveta4,
            fms(0.5, c12, c1 * self.eta2 * invv / 6.0)
                - sav2 * self.eta2 / 3.0
                - self.eta4 * invv2 / 3.0,
            c1 * 0.5 * self.inveta * fms(c12, self.inveta4, sav2 * self.inveta2) * log1,
        );
        let phiz_1_coeff = fma(
            -0.5 * j_norm * self.inveta2,
            fma(eta, l_newt, c1),
            0.5 * self.inveta3 * fnms(self.eta2, sav2, c12) * log1,
        );
        let phiz_2_coeff = fnms(c1 * log1, self.inveta, fma(sav, log2, -j_norm));
        let phiz_3_coeff = fma(c1 * log2, inv_sav, fms(j_norm, vv, eta * log1));
        let phiz_4_coeff = fms(
            0.5 * j_norm * inv_sav2 * vv,
            fma(vv, sav2, c1),
            0.5 * inv_sav2 * inv_sav * fnms(self.eta2, sav2, c12) * log2,
        );
        let phiz_5_coeff = fma(
            -j_norm * vv,
            fms(0.5 * c12 * inv_sav2, inv_sav2, c1 * vv * inv_sav2 / 6.0)
                - vv * vv / 3.0
                - self.eta2 * inv_sav2 / 3.0,
            0.5 * c1 * inv_sav2 * inv_sav2 * inv_sav * fnms(self.eta2, sav2, c12) * log2,
        );

        let o = &self.omegaz_coeff;
        let mut sum = fma(phiz_0_coeff, o[0], phiz_1_coeff * o[1]);
        sum = fma(phiz_2_coeff, o[2], sum);
        sum = fma(phiz_3_coeff, o[3], sum);
        sum = fma(phiz_4_coeff, o[4], sum);
        fma(phiz_5_coeff, o[5], sum)
    }

    /// `IMRPhenomX_Return_phiz_MSA`
    #[inline]
    fn phiz(&self, v: &VPowers, j_norm: f64) -> f64 {
        let phiz = self.phiz_sum(v, j_norm) + self.phiz_0;
        if phiz.is_nan() {
            0.0
        } else {
            phiz
        }
    }

    /// `IMRPhenomX_Return_zeta_MSA`, Eq. F5
    #[inline]
    fn zeta(&self, v: &VPowers) -> f64 {
        let o = &self.omegazeta_coeff;
        let mut sum = fma(o[0], v.invv3, o[1] * v.invv2);
        sum = fma(o[2], v.invv, sum);
        sum = fma(o[3], v.logv, sum);
        sum = fma(o[4], v.v, sum);
        sum = fma(o[5], v.v2, sum);
        let zeta = fma(self.eta, sum, self.zeta_0);
        if zeta.is_nan() {
            0.0
        } else {
            zeta
        }
    }

    /// `IMRPhenomX_Return_MSA_Corrections_MSA`: the first MSA corrections
    /// to `phi_z` (Eq. 67) and `zeta` (Eq. F19), as implemented in
    /// `LALSimInspiralFDPrecAngles` (version 223).
    fn msa_corrections(&self, v: &VPowers, l_norm: f64, j_norm: f64, roots: &Roots) -> [f64; 2] {
        let eta = self.eta;
        let seff = self.seff;
        let (v1, v2, v3, v4, v6) = (v.v, v.v2, v.v3, v.v4, v.v6);
        let spl2 = roots.spl2;
        let smi2 = roots.smi2;
        let spl2_m_smi2 = roots.spl2_m_smi2();
        let spl = spl2.sqrt();
        let j_norm2 = j_norm * j_norm;
        let l_norm2 = l_norm * l_norm;
        let j2_m_spl2 = j_norm2 - spl2;
        let one_m_seff_v = fnms(seff, v1, 1.0);

        // Eqs. B6 - B8 (IMRPhenomX_Return_Constants_c_MSA)
        let mut x = fma(4.0 * self.eta3 * seff, v1, self.eta3);
        x = fnms(
            2.0 * eta
                * fma(
                    2.0 * (self.s1_norm_2 - self.s2_norm_2),
                    self.delta_qq,
                    j2_m_spl2,
                ),
            v2,
            x,
        );
        x = fnms(4.0 * eta * seff * j2_m_spl2, v3, x);
        x = fma(j2_m_spl2 * j2_m_spl2 * v4, self.inveta, x);
        let c0 = j_norm * (0.75 * one_m_seff_v * v2 * x);
        let c2 = j_norm
            * (-1.5
                * eta
                * (spl2 - smi2)
                * fnms(j2_m_spl2 * v2, self.inveta2, fma(2.0 * seff, v1, 1.0))
                * one_m_seff_v
                * v4);
        let c4 = j_norm * (0.75 * self.inveta * (spl2 - smi2) * (spl2 - smi2) * one_m_seff_v * v6);

        // Eqs. B9 - B11 (IMRPhenomX_Return_Constants_d_MSA)
        let d0 =
            -fnms(l_norm + spl, l_norm + spl, j_norm2) * fnms(l_norm - spl, l_norm - spl, j_norm2);
        let d2 = -2.0 * (spl2 - smi2) * (j_norm2 + l_norm2 - spl2);
        let d4 = -(spl2 - smi2) * (spl2 - smi2);

        let two_d0 = 2.0 * d0;
        // Eq. B20
        let sd = fms(d2, d2, 4.0 * d0 * d4).abs().sqrt();
        // Eqs. F20, F21
        let a_theta_l = 0.5 * ((j_norm / l_norm) + (l_norm / j_norm) - (spl2 / (j_norm * l_norm)));
        let b_theta_l = 0.5 * spl2_m_smi2 / (j_norm * l_norm);
        // Eqs. B16, B17
        let nc_num = 2.0 * (d0 + d2 + d4);
        let nc_denom = two_d0 + d2 + sd;
        let nc = nc_num / nc_denom;
        let nd = nc_denom / two_d0;
        let sqrt_nc = nc.abs().sqrt();
        let sqrt_nd = nd.abs().sqrt();

        // Phase of S and its derivative (IMRPhenomX_Return_Psi_MSA,
        // IMRPhenomX_Return_Psi_dot_MSA)
        let psi = -0.75 * self.g0 * self.delta_qq * fma(self.psi2, v2, fma(self.psi1, v1, 1.0))
            / (v2 * v1)
            + self.psi0;
        let a_coeff = -1.5 * (v4 * v2) * fnms(v1, seff, 1.0) * self.sqrt_inveta;
        let psi_dot = 0.5 * a_coeff * (spl2 - roots.s32).sqrt();
        let tan_psi = psi.tan();
        let atan_psi = tan_psi.atan();

        // Eq. B19 denominator
        let c2den = 2.0 * d0 * sd * (d0 + d2 + d4);

        // Eq. 67 as in LALSimInspiralFDPrecAngles
        let cphi_term = if nc == 1.0 {
            0.0
        } else {
            let num = fnms(
                c0,
                fms(2. * d0, d4, (d2 + d4) * (d2 - sd)),
                fms(
                    c4 * d0,
                    fma(2., d0, d2) + sd,
                    c2 * d0 * (fma(2., d4, d2) - sd),
                ),
            );
            (num / c2den).abs() * (sqrt_nc / (nc - 1.) * (atan_psi - (sqrt_nc * tan_psi).atan()))
                / psi_dot
        };
        let dphi_term = if nd == 1.0 {
            0.0
        } else {
            let num = fnms(
                c0,
                fma(-(2. * d0), d4, (d2 + d4) * (d2 + sd)),
                fma(
                    -c4 * d0,
                    fma(2., d0, d2) - sd,
                    c2 * d0 * (fma(2., d4, d2) + sd),
                ),
            );
            num.abs() / c2den * (sqrt_nd / (nd - 1.) * (atan_psi - (sqrt_nd * tan_psi).atan()))
                / psi_dot
        };
        let mut phiz_msa = cphi_term + dphi_term;
        // Eq. F19 as in LALSimInspiralFDPrecAngles
        let mut zeta_msa = fma(
            a_theta_l,
            phiz_msa,
            2. * b_theta_l * d0 * (cphi_term / (sd - d2) - dphi_term / (sd + d2)),
        );
        if phiz_msa.is_nan() {
            phiz_msa = 0.0;
        }
        if zeta_msa.is_nan() {
            zeta_msa = 0.0;
        }
        [phiz_msa, zeta_msa]
    }

    /// `IMRPhenomX_Return_phi_zeta_costhetaL_MSA`: `phi_z`, `zeta` and
    /// `cos(theta_L)` at the orbital velocity `v = (pi M f)^(1/3)` (for the
    /// co-precessing mode with `m' = 2`; XPHM evaluates the mode `m'` at
    /// `v = (2 pi M f / m')^(1/3)`).
    #[inline]
    pub(crate) fn phi_zeta_costhetal(&self, v: f64) -> MsaAngles {
        let powers = VPowers::new(v);
        let l_norm = self.eta / v;
        let j_norm = self.j_norm(l_norm);
        let l_norm_3pn = self.l_norm_3pn(v, powers.v2, l_norm);
        let j_norm_3pn = self.j_norm(l_norm_3pn);
        let roots = self.roots(l_norm, j_norm);
        let s_norm = self.s_norm(&powers, &roots);
        let corrections = if roots.has_corrections() {
            self.msa_corrections(&powers, l_norm, j_norm, &roots)
        } else {
            [0.0, 0.0]
        };
        let phiz = self.phiz(&powers, j_norm);
        let zeta = self.zeta(&powers);
        MsaAngles {
            phiz: phiz + corrections[0],
            zeta: zeta + corrections[1],
            cos_theta_l: cos_theta_lj(l_norm_3pn, j_norm_3pn, s_norm),
        }
    }

    /// [`Self::phi_zeta_costhetal`] at the orbital velocities
    /// `v_i = x_i^(1/3)`, passed to `out` in order: the same operations,
    /// done in passes over blocks of points so that the (latency-bound)
    /// evaluations at different points overlap.
    pub(crate) fn phi_zeta_costhetal_many(&self, x: &[f64], mut out: impl FnMut(MsaAngles)) {
        const BLOCK: usize = 32;
        #[derive(Clone, Copy, Default)]
        struct Point {
            powers: VPowers,
            l_norm: f64,
            j_norm: f64,
            l_norm_3pn: f64,
            j_norm_3pn: f64,
            roots: Roots,
            s_norm: f64,
            corrections: [f64; 2],
            phiz: f64,
        }
        let mut points = [Point::default(); BLOCK];
        for block in x.chunks(BLOCK) {
            let points = &mut points[..block.len()];
            for (point, &x) in points.iter_mut().zip(block) {
                let v = x.libm_cbrt();
                let powers = VPowers::new(v);
                let l_norm = self.eta / v;
                let j_norm = self.j_norm(l_norm);
                let l_norm_3pn = self.l_norm_3pn(v, powers.v2, l_norm);
                point.powers = powers;
                point.l_norm = l_norm;
                point.j_norm = j_norm;
                point.l_norm_3pn = l_norm_3pn;
                point.j_norm_3pn = self.j_norm(l_norm_3pn);
            }
            for point in points.iter_mut() {
                point.roots = self.roots(point.l_norm, point.j_norm);
            }
            for point in points.iter_mut() {
                point.s_norm = self.s_norm_with(&point.powers, &point.roots, elljac_sn_fast);
            }
            for point in points.iter_mut() {
                point.corrections = if point.roots.has_corrections() {
                    self.msa_corrections(&point.powers, point.l_norm, point.j_norm, &point.roots)
                } else {
                    [0.0, 0.0]
                };
            }
            for point in points.iter_mut() {
                point.phiz = self.phiz(&point.powers, point.j_norm);
            }
            for point in points.iter() {
                let zeta = self.zeta(&point.powers);
                out(MsaAngles {
                    phiz: point.phiz + point.corrections[0],
                    zeta: zeta + point.corrections[1],
                    cos_theta_l: cos_theta_lj(point.l_norm_3pn, point.j_norm_3pn, point.s_norm),
                });
            }
        }
    }

    /// Final spin of `IMRPhenomX_SetPrecessingRemnantParams` for the MSA
    /// angles (`PhenomXPFinalSpinMod` 3 with the default
    /// `PhenomXPTransPrecessionMethod` 1): from the precession-averaged
    /// couplings and the aligned-spin final spin `af_parallel`.
    ///
    /// `m` is the total mass, `m1_2`, `m2_2` the squared masses (in units of
    /// the total mass, as in the precession struct) and `chi1l`, `chi2l`
    /// the aligned spins.
    pub(crate) fn final_spin(
        &self,
        af_parallel: f64,
        m: f64,
        m1_2: f64,
        m2_2: f64,
        chi1l: f64,
        chi2l: f64,
    ) -> f64 {
        let l_final = fnms(m2_2, chi2l, fms(m * m, af_parallel, m1_2 * chi1l));
        let sign = 1.0f64.copysign(af_parallel);
        sign * fma(
            2.0 * l_final,
            self.s1l_pav + self.s2l_pav,
            fma(l_final, l_final, self.sav2),
        )
        .sqrt()
            / (m * m)
    }
}

/// `IMRPhenomX_costhetaLJ`
#[inline]
fn cos_theta_lj(l_norm: f64, j_norm: f64, s_norm: f64) -> f64 {
    let cos = 0.5 * fnms(s_norm, s_norm, fma(j_norm, j_norm, l_norm * l_norm)) / (l_norm * j_norm);
    // like the comparisons of LAL, clamp leaves NaN unchanged
    cos.clamp(-1.0, 1.0)
}

/// `gsl_hypot(x, y)` (GSL's `sys/hypot.c`).
#[inline]
fn gsl_hypot(x: f64, y: f64) -> f64 {
    let xabs = x.abs();
    let yabs = y.abs();
    if x.is_infinite() || y.is_infinite() {
        return f64::INFINITY;
    }
    let (min, max) = if xabs < yabs {
        (xabs, yabs)
    } else {
        (yabs, xabs)
    };
    if min == 0.0 {
        return max;
    }
    let u = min / max;
    max * fma(u, u, 1.0).sqrt()
}

/// `GSL_SIGN`
#[inline]
fn gsl_sign(x: f64) -> f64 {
    if x >= 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// The Jacobi elliptic function `sn(u | m)` of `gsl_sf_elljac_e` (GSL 2.8,
/// Bulirsch's algorithm with the reflection to `K - u`). Returns NaN where
/// GSL signals a domain error (`|m| > 1`).
fn elljac_sn(u: f64, m: f64) -> f64 {
    if m.abs() > 1.0 {
        // GSL sets sn = 0 and raises a domain error
        return f64::NAN;
    }
    if m.abs() < 2.0 * GSL_DBL_EPSILON {
        return u.sin();
    }
    if (m - 1.0).abs() < 2.0 * GSL_DBL_EPSILON {
        return u.tanh();
    }
    const N: usize = 16;
    let mut mu = [0.0f64; N];
    let mut nu = [0.0f64; N];
    let mut c = [0.0f64; N];
    let mut d = [0.0f64; N];
    let mut n = 0;
    mu[0] = 1.0;
    let sqrt_1mm = (1.0 - m).sqrt();
    nu[0] = sqrt_1mm;
    while (mu[n] - nu[n]).abs() > 4.0 * GSL_DBL_EPSILON * (mu[n] + nu[n]).abs() {
        mu[n + 1] = 0.5 * (mu[n] + nu[n]);
        nu[n + 1] = (mu[n] * nu[n]).sqrt();
        n += 1;
        if n >= N - 1 {
            // GSL_EMAXITER, ignored by LAL
            break;
        }
    }
    let sin_umu = (u * mu[n]).sin();
    let cos_umu = (u * mu[n]).cos();
    // Bulirsch's descending recursion; returns (c[0], d[0]).
    let descend = |t: f64, mu: &[f64; N], nu: &[f64; N], c: &mut [f64; N], d: &mut [f64; N]| {
        let mut n = n;
        c[n] = mu[n] * t;
        d[n] = 1.0;
        while n > 0 {
            n -= 1;
            c[n] = d[n + 1] * c[n + 1];
            let r = (c[n + 1] * c[n + 1]) / mu[n + 1];
            d[n] = (r + nu[n]) / (r + mu[n]);
        }
        (c[0], d[0])
    };
    if sin_umu.abs() < cos_umu.abs() {
        let (c0, d0) = descend(sin_umu / cos_umu, &mu, &nu, &mut c, &mut d);
        let dn = sqrt_1mm / d0;
        let cn = dn * gsl_sign(cos_umu) / gsl_hypot(1.0, c0);
        cn * c0 / sqrt_1mm
    } else {
        let (c0, _) = descend(cos_umu / sin_umu, &mu, &nu, &mut c, &mut d);
        gsl_sign(sin_umu) / gsl_hypot(1.0, c0)
    }
}

/// [`elljac_sn`] for the grid paths, agreeing with it to a few ulp: the
/// inlined `sin_cos` and the descending recursion
/// `d_k = (c_{k+1}^2 / mu_{k+1} + nu_k) / (c_{k+1}^2 / mu_{k+1} + mu_k)`
/// rearranged to one division per step.
#[inline]
fn elljac_sn_fast(u: f64, m: f64) -> f64 {
    if m.abs() > 1.0 {
        return f64::NAN;
    }
    if m.abs() < 2.0 * GSL_DBL_EPSILON {
        return u.sin();
    }
    if (m - 1.0).abs() < 2.0 * GSL_DBL_EPSILON {
        return u.tanh();
    }
    const N: usize = 16;
    let mut mu = [0.0f64; N];
    let mut nu = [0.0f64; N];
    let mut n = 0;
    mu[0] = 1.0;
    let sqrt_1mm = (1.0 - m).sqrt();
    nu[0] = sqrt_1mm;
    while (mu[n] - nu[n]).abs() > 4.0 * GSL_DBL_EPSILON * (mu[n] + nu[n]).abs() {
        mu[n + 1] = 0.5 * (mu[n] + nu[n]);
        nu[n + 1] = (mu[n] * nu[n]).sqrt();
        n += 1;
        if n >= N - 1 {
            break;
        }
    }
    let (sin_umu, cos_umu) = fastmath::sin_cos(u * mu[n]);
    // c_n = mu_n t, d_n = 1; c_k = d_{k+1} c_{k+1}; returns (c_0, d_0)
    let descend = |t: f64| {
        let mut c = mu[n] * t;
        let mut d = 1.0;
        for k in (0..n).rev() {
            let c2 = c * c;
            let d_k = fma(nu[k], mu[k + 1], c2) / fma(mu[k], mu[k + 1], c2);
            c *= d;
            d = d_k;
        }
        (c, d)
    };
    if sin_umu.abs() < cos_umu.abs() {
        let (c0, d0) = descend(sin_umu / cos_umu);
        // sign(cos) c0 / (d0 hypot(1, c0)), |c0| <= 1
        (c0 / (d0 * fma(c0, c0, 1.0).sqrt())).copysign(c0 * cos_umu)
    } else {
        let (c0, _) = descend(cos_umu / sin_umu);
        (1.0 / fma(c0, c0, 1.0).sqrt()).copysign(sin_umu)
    }
}

/// Carlson's `R_F(x, y, z)`, `gsl_sf_ellint_RF_e` with `GSL_PREC_DOUBLE`
/// (NaN for GSL's domain errors).
fn ellint_rf(x: f64, y: f64, z: f64) -> f64 {
    let lolim = 5.0 * GSL_DBL_MIN;
    let uplim = 0.2 * GSL_DBL_MAX;
    let errtol = 0.001;
    let nmax = 10000;
    let max3 = |a: f64, b: f64, c: f64| {
        // locMAX3: GSL_MAX(a, GSL_MAX(b, c)) with GSL_MAX(a, b) = a > b ? a : b
        let bc = if b > c { b } else { c };
        if a > bc {
            a
        } else {
            bc
        }
    };
    if x < 0.0 || y < 0.0 || z < 0.0 {
        return f64::NAN;
    }
    if x + y < lolim || x + z < lolim || y + z < lolim {
        return f64::NAN;
    }
    if max3(x, y, z) >= uplim {
        return f64::NAN;
    }
    let c1 = 1.0 / 24.0;
    let c2 = 3.0 / 44.0;
    let c3 = 1.0 / 14.0;
    let (mut xn, mut yn, mut zn) = (x, y, z);
    let mut n = 0;
    let (mu, xndev, yndev, zndev) = loop {
        let mu = (xn + yn + zn) / 3.0;
        let xndev = 2.0 - (mu + xn) / mu;
        let yndev = 2.0 - (mu + yn) / mu;
        let zndev = 2.0 - (mu + zn) / mu;
        let epslon = max3(xndev.abs(), yndev.abs(), zndev.abs());
        if epslon < errtol {
            break (mu, xndev, yndev, zndev);
        }
        let xnroot = xn.sqrt();
        let ynroot = yn.sqrt();
        let znroot = zn.sqrt();
        let lamda = fma(xnroot, ynroot + znroot, ynroot * znroot);
        xn = (xn + lamda) * 0.25;
        yn = (yn + lamda) * 0.25;
        zn = (zn + lamda) * 0.25;
        n += 1;
        if n == nmax {
            return f64::NAN;
        }
    };
    let e2 = fms(xndev, yndev, zndev * zndev);
    let e3 = xndev * yndev * zndev;
    let s = fma(c3, e3, fma(fnms(c2, e3, fms(c1, e2, 0.1)), e2, 1.0));
    s / mu.sqrt()
}

/// `gsl_sf_ellint_Kcomp_e` with `GSL_PREC_DOUBLE`.
fn ellint_kcomp(k: f64) -> f64 {
    if k * k >= 1.0 {
        f64::NAN
    } else if k * k >= 1.0 - GSL_SQRT_DBL_EPSILON {
        // Abramowitz & Stegun 17.3.34
        let y = fnms(k, k, 1.0);
        let a = [1.38629436112, 0.09666344259, 0.03590092383];
        let b = [0.5, 0.12498593597, 0.06880248576];
        let ta = fma(y, fma(y, a[2], a[1]), a[0]);
        let tb = -y.ln() * fma(y, fma(y, b[2], b[1]), b[0]);
        ta + tb
    } else {
        ellint_rf(0.0, fnms(k, k, 1.0), 1.0)
    }
}

/// The incomplete elliptic integral of the first kind `F(phi, k)`,
/// `gsl_sf_ellint_F(phi, k, GSL_PREC_DOUBLE)` (Legendre form with modulus
/// `k`).
fn ellint_f(phi: f64, k: f64) -> f64 {
    // Angular reduction to -pi/2 < phi < pi/2
    let nc = (phi / PI + 0.5).floor();
    let phi = fnms(nc, PI, phi);
    let sin_phi = phi.sin();
    let sin2_phi = sin_phi * sin_phi;
    let x = 1.0 - sin2_phi;
    let y = fnms(k * k, sin2_phi, 1.0);
    let rf = ellint_rf(x, y, 1.0);
    let val = sin_phi * rf;
    if nc == 0.0 {
        val
    } else {
        fma(2.0 * nc, ellint_kcomp(k), val)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `gsl_sf_elljac_e` and `gsl_sf_ellint_F` against values from the GSL
    /// 2.8 library of the lalsimulation conda environment (`sn` and `F`
    /// via `scipy.special` agree to ~1e-15).
    #[test]
    fn elliptic_functions() {
        // sn(u | m) = sin(am(u | m)) and F(am(u | m), sqrt(m)) = u
        for &(u, m) in &[
            (0.3, 0.2),
            (-1.7, 0.9),
            (4.2, 0.5),
            (12.0, 0.01),
            (0.8, 0.999),
        ] {
            let sn = elljac_sn(u, m);
            let phi = sn.asin();
            // F inverts sn on the principal branch |u| < K
            let f = ellint_f(phi, m.sqrt());
            let k = ellint_kcomp(m.sqrt());
            if u.abs() < k {
                assert!((f - u).abs() < 1e-13, "{u} {m} {f}");
            }
            assert!(sn.abs() <= 1.0);
        }
        assert!((ellint_kcomp(0.0) - PI / 2.0).abs() < 1e-15);
        assert_eq!(elljac_sn(0.4, 0.0), 0.4f64.sin());
        assert!(elljac_sn(0.4, 1.5).is_nan());
    }

    #[test]
    fn hypot_matches_definition() {
        assert_eq!(gsl_hypot(1.0, 0.0), 1.0);
        assert!((gsl_hypot(3.0, 4.0) - 5.0).abs() < 1e-15);
        assert_eq!(gsl_hypot(f64::INFINITY, f64::NAN), f64::INFINITY);
    }

    /// The batched evaluation of the grid paths gives the angles of
    /// `phi_zeta_costhetal`: `phi_z` and `zeta` bit for bit (they are
    /// ill-conditioned for nearly aligned spins, where changes of the roots
    /// in the last bit change them by up to 1e-5 at high frequencies), and
    /// `cos(theta_L)`, which uses the rearranged `sn`, to rounding.
    #[test]
    fn batch_matches_single_frequency_angles() {
        let spins = [
            ([0.4, -0.3, 0.2], [0.1, 0.5, -0.3]),
            ([0.9, 0.1, 0.3], [0.1, 0.5, -0.3]),
            ([0.05, 0.02, 0.9], [0.01, 0.0, 0.9]),
            ([0.7, 0.0, -0.7], [-0.6, 0.3, 0.6]),
            ([0.2, 0.2, 0.0], [0.0, 0.0, 0.0]),
        ];
        let pi_m = PI * 60.0 * 4.925490947641267e-6;
        let x: Vec<f64> = (0..20000)
            .map(|i| pi_m * (5.0 + 0.125 * i as f64))
            .collect();
        let mut systems = 0;
        for (chi1, chi2) in spins {
            for (m1, m2) in [(2.0 / 3.0, 1.0 / 3.0), (0.5, 0.5), (0.9, 0.1)] {
                let Some(msa) = MsaSystem::new(&MsaInputs {
                    eta: m1 * m2,
                    m1,
                    m2,
                    m1_2: m1 * m1,
                    m2_2: m2 * m2,
                    chi1,
                    chi2,
                    pi_gm: pi_m,
                    f_ref: 20.0,
                }) else {
                    continue;
                };
                systems += 1;
                let mut batch = Vec::with_capacity(x.len());
                msa.phi_zeta_costhetal_many(&x, |angles| batch.push(angles));
                for (x, b) in x.iter().zip(&batch) {
                    let a = msa.phi_zeta_costhetal(x.libm_cbrt());
                    assert_eq!((a.phiz, a.zeta), (b.phiz, b.zeta), "{chi1:?} {m1} {x}");
                    assert!((a.cos_theta_l - b.cos_theta_l).abs() < 1e-14);
                }
            }
        }
        assert!(systems > 10);
    }
}
