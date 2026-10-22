// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Phase of the (2,2) mode: `IMRPhenomXPhaseCoefficients`,
//! `IMRPhenomXGetPhaseCoefficients`, `IMRPhenomX_Phase_22_ConnectionCoefficients`,
//! `IMRPhenomX_TimeShift_22`, `IMRPhenomX_Phase_22` and `IMRPhenomX_dPhase_22`
//! from `LALSimIMRPhenomX_internals.c`.

use std::f64::consts::{LN_2, PI};

use crate::constants::EGAMMA;
use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct, POWERS_OF_LALPI};
use crate::imrphenomx::jet::Jet;
use crate::imrphenomx::linalg::lu_solve;
use crate::imrphenomx::utilities::{linb, psi4_to_strain};
use crate::imrphenomx::{inspiral, intermediate, ringdown};

/// Equivalent of `IMRPhenomXPhaseCoefficients` for the default versions
/// (InsPhaseVersion 104, IntPhaseVersion 105, RDPhaseVersion 105) in GR.
///
/// Besides the LAL fields this caches a few combinations of coefficients
/// used in the per-frequency evaluation (`b3_ov2`, `inv_f_damp`, ...).
#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub(crate) struct PhaseCoefficients {
    // Transition frequencies (dimensionless)
    pub(crate) f_phase_ins_min: f64,
    pub(crate) f_phase_ins_max: f64,
    pub(crate) f_phase_rd_min: f64,
    pub(crate) f_phase_rd_max: f64,
    /// Inspiral-intermediate transition.
    pub(crate) f_phase_match_in: f64,
    /// Intermediate-ringdown transition.
    pub(crate) f_phase_match_im: f64,

    /// Connection coefficients (zero until
    /// [`phase_22_connection_coefficients`] is called).
    pub(crate) c1_int: f64,
    pub(crate) c2_int: f64,
    pub(crate) c1_mrd: f64,
    pub(crate) c2_mrd: f64,

    /// Ringdown coefficients.
    pub(crate) c0: f64,
    pub(crate) c1: f64,
    pub(crate) c2: f64,
    pub(crate) c4: f64,
    pub(crate) c_rd: f64,
    pub(crate) c_l: f64,
    pub(crate) c_l_gr: f64,
    pub(crate) c4_ov3: f64,
    pub(crate) c_l_ov_fda: f64,

    /// Intermediate coefficients.
    pub(crate) b0: f64,
    pub(crate) b1: f64,
    pub(crate) b2: f64,
    pub(crate) b3: f64,
    pub(crate) b4: f64,

    /// Inspiral pseudo-PN coefficients of the phase derivative.
    pub(crate) a0: f64,
    pub(crate) a1: f64,
    pub(crate) a2: f64,
    pub(crate) a3: f64,
    pub(crate) a4: f64,
    /// Inspiral pseudo-PN coefficients of the phase.
    pub(crate) sigma1: f64,
    pub(crate) sigma2: f64,
    pub(crate) sigma3: f64,
    pub(crate) sigma4: f64,
    pub(crate) sigma5: f64,

    /// TaylorF2 phase coefficients normalised by 3 / (128 eta pi^(5/3) f^(5/3)),
    /// including powers of pi.
    pub(crate) phi0: f64,
    pub(crate) phi1: f64,
    pub(crate) phi2: f64,
    pub(crate) phi3: f64,
    pub(crate) phi4: f64,
    pub(crate) phi5: f64,
    pub(crate) phi5l: f64,
    pub(crate) phi6: f64,
    pub(crate) phi6l: f64,
    pub(crate) phi7: f64,
    pub(crate) phi8: f64,
    pub(crate) phi8l: f64,
    pub(crate) phi9: f64,
    pub(crate) phi9l: f64,
    pub(crate) phi_initial: f64,
    pub(crate) phi_norm: f64,

    /// Corresponding coefficients of the phase derivative.
    pub(crate) dphi0: f64,
    pub(crate) dphi1: f64,
    pub(crate) dphi2: f64,
    pub(crate) dphi3: f64,
    pub(crate) dphi4: f64,
    pub(crate) dphi5: f64,
    pub(crate) dphi6: f64,
    pub(crate) dphi6l: f64,
    pub(crate) dphi7: f64,
    pub(crate) dphi8: f64,
    pub(crate) dphi8l: f64,
    pub(crate) dphi9: f64,
    pub(crate) dphi9l: f64,

    pub(crate) collocation_points_phase_rd: [f64; 5],
    pub(crate) collocation_values_phase_rd: [f64; 5],
    pub(crate) collocation_points_phase_int: [f64; 5],
    pub(crate) collocation_values_phase_int: [f64; 5],
    pub(crate) collocation_points_phase_ins: [f64; 4],
    pub(crate) collocation_values_phase_ins: [f64; 4],

    // Cached quantities for the per-frequency evaluation.
    /// Ringdown frequency copied from the waveform struct.
    pub(crate) f_ring: f64,
    /// Damping frequency copied from the waveform struct.
    pub(crate) f_damp: f64,
    /// f_damp^2
    pub(crate) f_damp2: f64,
    /// 1 / f_damp
    pub(crate) inv_f_damp: f64,
    /// 1 / (2 f_damp)
    pub(crate) inv_2_f_damp: f64,
    /// 1.5 c1
    pub(crate) c1_3ov2: f64,
    /// b3 / 2
    pub(crate) b3_ov2: f64,
    /// b4 / 3
    pub(crate) b4_ov3: f64,
    /// 4 cL (GR value), numerator of the intermediate Lorentzian.
    pub(crate) int_lorentzian_numerator: f64,
    /// 4 f_damp^2, width of the intermediate Lorentzian.
    pub(crate) int_lorentzian_width2: f64,
    /// 2 cL / f_damp (GR value), prefactor of the intermediate arctangent.
    pub(crate) int_atan_prefactor: f64,
    /// 5 / (128 pi^(5/3)), normalisation of the inspiral phase derivative.
    pub(crate) dphase_norm: f64,
}

/// Gauss-Chebyshev collocation nodes on [0, 1] for four points.
const GPOINTS4: [f64; 4] = [0.0, 1.0 / 4.0, 3.0 / 4.0, 1.0];

/// Gauss-Chebyshev collocation nodes on [0, 1] for five points.
fn gpoints5() -> [f64; 5] {
    let s = 1.0 / (2.0 * 2.0f64.sqrt());
    [0.0, 1.0 / 2.0 - s, 1.0 / 2.0, 1.0 / 2.0 + s, 1.0]
}

impl PhaseCoefficients {
    /// `IMRPhenomXGetPhaseCoefficients`. The connection coefficients are
    /// left at zero; call [`phase_22_connection_coefficients`] afterwards.
    pub(crate) fn new(wf: &WaveformStruct) -> Self {
        let pi = &*POWERS_OF_LALPI;
        let gpoints5 = gpoints5();
        let f_ring = wf.f_ring;
        let f_damp = wf.f_damp;
        let dphase0 = wf.dphase0;

        // Matching regions, Eqs. 5.10, 5.11 and 7.7 of arXiv:2001.11412.
        let f_im_match = 0.6 * (0.5 * f_ring + wf.f_isco);
        let f_in_match = wf.f_meco;
        let deltaf = (f_im_match - f_in_match) * 0.03;
        let f_phase_match_in = f_in_match - 1.0 * deltaf;
        let f_phase_match_im = f_im_match + 0.5 * deltaf;
        let f_phase_ins_min = 0.0026;
        let f_phase_ins_max = 1.020 * wf.f_meco;
        let f_phase_rd_min = f_im_match;
        let f_phase_rd_max = fma(1.25, f_damp, f_ring);

        let phi_norm = -(3.0 * pi.m_five_thirds) / (128.0);

        let chi1l = wf.chi1l;
        let chi2l = wf.chi2l;
        let chi1l2l = chi1l * chi2l;
        let chi1l2 = wf.chi1l * wf.chi1l;
        let chi1l3 = wf.chi1l * chi1l2;
        let chi2l2 = wf.chi2l * wf.chi2l;
        let chi2l3 = wf.chi2l * chi2l2;
        let eta = wf.eta;
        let eta2 = eta * eta;
        let eta3 = eta * eta2;
        let delta = wf.delta;
        let (s, dchi) = (wf.s_tot_r, wf.dchi);

        /* ---------------------------- Ringdown ---------------------------- */
        let deltax = f_phase_rd_max - f_phase_rd_min;
        let xmin = f_phase_rd_min;
        let mut points_rd = [0.0; 5];
        for (point, g) in points_rd.iter_mut().zip(gpoints5) {
            *point = fma(g, deltax, xmin);
        }
        // Collocation point 4 is set to the ringdown frequency ~ dip in Lorentzian.
        points_rd[3] = f_ring;

        let rd_v4 = ringdown::phase_22_v4(eta, s, dchi, wf.delta);
        let mut values_rd = [
            ringdown::phase_22_d12(eta, s, dchi, wf.delta),
            ringdown::phase_22_d24(eta, s, dchi, wf.delta),
            ringdown::phase_22_d34(eta, s, dchi, wf.delta),
            rd_v4,
            ringdown::phase_22_d54(eta, s, dchi, wf.delta),
        ];
        // v_j = d_{j4} + v4
        values_rd[4] += values_rd[3];
        values_rd[2] += values_rd[3];
        values_rd[1] += values_rd[3];
        values_rd[0] += values_rd[1];
        let phase_rd = values_rd[0];

        // ansatzRD(f) = a0 + a1 f^(-1/3) + a2 f^(-2) + a4 f^(-4) + aRD / (f_damp^2 + (f - f_ring)^2)
        let mut matrix = [[0.0; 5]; 5];
        for (row, &ff) in matrix.iter_mut().zip(&points_rd) {
            let invff = 1.0 / ff;
            let ff1 = invff.cbrt();
            let ff2 = invff * invff;
            let ff3 = ff2 * ff2;
            let ff4 = -dphase0 / fma(f_damp, f_damp, (ff - f_ring) * (ff - f_ring));
            *row = [1.0, ff1, ff2, ff3, ff4];
        }
        let [c0, c1, c2, c4, c_rd] = lu_solve(matrix, values_rd);
        let c_l = -(dphase0 * c_rd);
        // IMRPhenomXPNR co-precessing deviation (zero otherwise)
        let c_l = fma(wf.coprec.dev, wf.coprec.nu4, c_l);

        /* ---------------------------- Inspiral ---------------------------- */
        let deltax = f_phase_ins_max - f_phase_ins_min;
        let xmin = f_phase_ins_min;
        let mut points_ins = [0.0; 4];
        for (point, g) in points_ins.iter_mut().zip(GPOINTS4) {
            *point = fma(g, deltax, xmin);
        }
        let s_ins = wf.chi_pn_hat;
        let v3 = inspiral::phase_22_v3(eta, s_ins, dchi, wf.delta);
        // v_i = d_i3 + v3
        let values_ins = [
            inspiral::phase_22_d13(eta, s_ins, dchi, wf.delta) + v3,
            inspiral::phase_22_d23(eta, s_ins, dchi, wf.delta) + v3,
            v3,
            inspiral::phase_22_d43(eta, s_ins, dchi, wf.delta) + v3,
        ];
        let mut matrix = [[0.0; 4]; 4];
        for (row, &ff) in matrix.iter_mut().zip(&points_ins) {
            let ff1 = ff.cbrt();
            let ff2 = ff1 * ff1;
            *row = [1.0, ff1, ff2, ff];
        }
        let [a0, a1, a2, a3] = lu_solve(matrix, values_ins);
        let a4 = 0.0;

        // The pseudo-PN coefficients are normalised such that (dphase0 / eta) f^{8/3} a_j,
        // so they are rescaled by f^{-8/3} in the PN phasing.
        let sigma1 = (-5.0 / 3.0) * a0;
        let sigma2 = (-5.0 / 4.0) * a1;
        let sigma3 = (-5.0 / 5.0) * a2;
        let sigma4 = (-5.0 / 6.0) * a3;
        let sigma5 = (-5.0 / 7.0) * a4;

        /* --------------------- TaylorF2 PN coefficients --------------------- */
        // Normalised by 3 / (128 eta [pi M f]^{5/3}), including powers of pi.
        let phi0_ns = 1.0;
        let phi2_ns = (3715.0 / 756.0 + 55.0 * eta / 9.0) * pi.two_thirds;
        let phi3_ns = -16.0 * pi.two;
        let phi3_s = fma(
            113.0,
            fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)),
            -(76.0 * (chi1l + chi2l) * eta),
        ) / 6.0
            * pi.itself;
        let phi4_ns =
            (15293365.0 / 508032.0 + 27145.0 * eta / 504.0 + 3085.0 * eta2 / 72.0) * pi.four_thirds;
        let phi4_s =
            -5.0 * fma(
                -(81.0 * chi2l2),
                fma(2.0, eta, -1.0 + delta),
                fma(
                    81.0 * chi1l2,
                    fma(-2.0, eta, 1.0 + delta),
                    316.0 * chi1l2l * eta,
                ),
            ) / 16.0
                * pi.four_thirds;
        let phi5l_ns = 5.0 * fma(-6552.0, eta, 46374.0) * PI / 4536.0 * pi.five_thirds;
        let phi5l_s = fma(
            85680.0 * (chi1l + chi2l),
            eta2,
            fma(
                -732985.0,
                fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)),
                -(560.0 * fma(-1213.0, chi1l + chi2l, 63.0 * (chi1l - chi2l) * delta) * eta),
            ),
        ) / 4536.0
            * pi.five_thirds;
        let phi6_ns = (11583231236531.0 / 4.69421568e9
            - 5.0 * eta * fma(588.0 * eta, fma(102260.0, eta, -45633.0), 3147553127.0)
                / 3.048192e6
            - 6848.0 * EGAMMA / 21.0
            - 640.0 * pi.two / 3.0
            + 2255.0 * eta * pi.two / 12.0
            - 13696.0 * LN_2 / 21.0
            - 6848.0 * pi.log / 63.0)
            * pi.two;
        let mut phi6_s =
            5.0 * fma(
                227.0,
                fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)),
                -(156.0 * (chi1l + chi2l) * eta),
            ) * PI
                / 3.0
                * pi.two;
        phi6_s = fma(
            5.0 * fma(
                -(7.0 * chi1l2),
                fma(
                    9808.0,
                    eta2,
                    fma(
                        -15103.0,
                        1.0 + delta,
                        2.0 * fma(6580.0, delta, 21683.0) * eta,
                    ),
                ),
                fma(
                    20.0 * chi1l2l * eta,
                    fma(12488.0, eta, 11763.0),
                    7.0 * chi2l2
                        * fma(
                            -9808.0,
                            eta2,
                            fma(
                                -15103.0,
                                -1.0 + delta,
                                2.0 * fma(6580.0, delta, -21683.0) * eta,
                            ),
                        ),
                ),
            ) / 4032.0,
            pi.two,
            phi6_s,
        );
        let phi6l_ns = -6848.0 / 63.0 * pi.two;
        let phi7_ns = 5.0 * fma(168.0 * fma(-29618.0, eta, 75703.0), eta, 15419335.0) * PI
            / 254016.0
            * pi.seven_thirds;
        let mut phi7_s =
            5.0 * fma(
                90514368.0 * (chi1l + chi2l),
                eta3,
                fma(
                    -(1008.0 * fma(208433.0, chi1l + chi2l, 25011.0 * (chi1l - chi2l) * delta)),
                    eta2,
                    fma(
                        -5030016755.0,
                        fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)),
                        4.0 * fma(
                            2113331119.0,
                            chi1l + chi2l,
                            675484362.0 * (chi1l - chi2l) * delta,
                        ) * eta,
                    ),
                ),
            ) / 6.096384e6
                * pi.seven_thirds;
        phi7_s = fma(
            -5.0 * fma(
                -(57.0 * chi2l2),
                fma(2.0, eta, -1.0 + delta),
                fma(
                    57.0 * chi1l2,
                    fma(-2.0, eta, 1.0 + delta),
                    220.0 * chi1l2l * eta,
                ),
            ) * PI,
            pi.seven_thirds,
            phi7_s,
        );
        phi7_s = fma(
            fma(
                40.0 * (chi1l + chi2l) * fma(17.0, chi2l2, fma(17.0, chi1l2, -(14.0 * chi1l2l))),
                eta2,
                fma(
                    14585.0,
                    fma(-chi2l3, -1.0 + delta, chi1l3 * (1.0 + delta)),
                    -(5.0
                        * fma(
                            chi1l3,
                            fma(2985.0, delta, 8819.0),
                            fma(
                                -(8439.0 * chi1l2 * chi2l),
                                1.0 + delta,
                                fma(
                                    chi2l3,
                                    fma(-2985.0, delta, 8819.0),
                                    8439.0 * chi1l * chi2l2 * (-1.0 + delta),
                                ),
                            ),
                        )
                        * eta),
                ),
            ) / 48.0,
            pi.seven_thirds,
            phi7_s,
        );
        let phi8_s =
            -5.0 * fma(
                94344.0 * (chi1l + chi2l),
                eta2,
                fma(
                    1263141.0,
                    fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)),
                    -(2.0 * fma(794075.0, chi1l + chi2l, 178533.0 * (chi1l - chi2l) * delta) * eta),
                ),
            ) * PI
                * (-1.0 + pi.log)
                / 9072.0
                * pi.eight_thirds;
        let phi8l_s =
            -5.0 * fma(
                94344.0 * (chi1l + chi2l),
                eta2,
                fma(
                    1263141.0,
                    fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)),
                    -(2.0 * fma(794075.0, chi1l + chi2l, 178533.0 * (chi1l - chi2l) * delta) * eta),
                ),
            ) * PI
                / 9072.0
                * pi.eight_thirds;

        let phi0 = phi0_ns;
        let phi1 = 0.0;
        let phi2 = phi2_ns;
        let phi3 = phi3_ns + phi3_s;
        let phi4 = phi4_ns + phi4_s;
        let phi5 = 0.0;
        let phi5l = phi5l_ns + phi5l_s;
        let phi6 = phi6_ns + phi6_s;
        let phi6l = phi6l_ns;
        let phi7 = phi7_ns + phi7_s;
        let phi8 = phi8_s;
        let phi8l = phi8l_s;
        let phi9 = 0.0;
        let phi9l = 0.0;

        // Phase derivative coefficients
        let dphi0 = phi0;
        let dphi1 = 4.0 / 5.0 * phi1;
        let dphi2 = 3.0 / 5.0 * phi2;
        let dphi3 = 2.0 / 5.0 * phi3;
        let dphi4 = 1.0 / 5.0 * phi4;
        let dphi5 = -3.0 / 5.0 * phi5l;
        let dphi6 = fma(-(1.0 / 5.0), phi6, -(3.0 / 5.0 * phi6l));
        let dphi6l = -1.0 / 5.0 * phi6l;
        let dphi7 = -2.0 / 5.0 * phi7;
        let dphi8 = fma(-(3.0 / 5.0), phi8, -(3.0 / 5.0 * phi8l));
        let dphi8l = -3.0 / 5.0 * phi8l;
        let dphi9 = fma(-(4.0 / 5.0), phi9, -(3.0 / 5.0 * phi9l));
        let dphi9l = -3.0 / 5.0 * phi9l;

        // Phase derivative at the inspiral-intermediate transition, the first
        // collocation point of the intermediate fit.
        let p = UsefulPowers::new(f_phase_match_in);
        let mut phase_in = dphi0;
        phase_in = fma(dphi1, p.one_third, phase_in);
        phase_in = fma(dphi2, p.two_thirds, phase_in);
        phase_in = fma(dphi3, p.itself, phase_in);
        phase_in = fma(dphi4, p.four_thirds, phase_in);
        phase_in = fma(dphi5, p.five_thirds, phase_in);
        phase_in = fma(dphi6, p.two, phase_in);
        phase_in = fma(dphi6l * p.two, p.log, phase_in);
        phase_in = fma(dphi7, p.seven_thirds, phase_in);
        phase_in = fma(dphi8, p.eight_thirds, phase_in);
        phase_in = fma(dphi8l * p.eight_thirds, p.log, phase_in);
        phase_in = fma(dphi9, p.three, phase_in);
        phase_in = fma(dphi9l * p.three, p.log, phase_in);
        let mut pseudo = fma(a0, p.eight_thirds, a1 * p.three);
        pseudo = fma(a2 * p.eight_thirds, p.two_thirds, pseudo);
        pseudo = fma(a3 * p.eight_thirds, p.itself, pseudo);
        pseudo = fma(a4 * p.eight_thirds, p.four_thirds, pseudo);
        phase_in += pseudo;
        phase_in = phase_in * p.m_eight_thirds * dphase0;

        /* -------------------------- Intermediate -------------------------- */
        let deltax = f_phase_match_im - f_phase_match_in;
        let xmin = f_phase_match_in;
        let mut points_int = [0.0; 5];
        for (point, g) in points_int.iter_mut().zip(gpoints5) {
            *point = fma(g, deltax, xmin);
        }
        let v2_im_m_rd_v4 = intermediate::phase_22_v2m_rdv4(eta, s, dchi, wf.delta);
        let v3_im_m_rd_v4 = intermediate::phase_22_v3m_rdv4(eta, s, dchi, wf.delta);
        let v2_im = intermediate::phase_22_v2(eta, s, dchi, wf.delta);
        let mut values_int = [
            phase_in,
            // Weighted average of the direct and conditioned fits.
            fma(0.75, v2_im_m_rd_v4 + rd_v4, 0.25 * v2_im),
            v3_im_m_rd_v4 + rd_v4,
            intermediate::phase_22_d43(eta, s, dchi, wf.delta),
            phase_rd,
        ];
        // v4 = d43 + v3
        values_int[3] += values_int[2];

        let mut matrix = [[0.0; 5]; 5];
        let mut rhs = [0.0; 5];
        for (i, &ff) in points_int.iter().enumerate() {
            let ff1 = f_ring / ff;
            let ff2 = ff1 * ff1;
            let ff3 = ff1 * ff2;
            let ff4 = ff2 * ff2;
            let ff0 = if i == 0 {
                (4.0 * c_l) / fma(2.0 * f_damp, 2.0 * f_damp, (ff - f_ring) * (ff - f_ring))
            } else {
                (4. * c_l) / fma(4.0 * f_damp, f_damp, (ff - f_ring) * (ff - f_ring))
            };
            matrix[i] = [1.0, ff1, ff2, ff3, ff4];
            rhs[i] = values_int[i] - ff0;
        }
        let x = lu_solve(matrix, rhs);
        let b0 = x[0];
        let b1 = x[1] * f_ring;
        let b2 = x[2] * f_ring * f_ring;
        let b3 = x[3] * f_ring * f_ring * f_ring;
        let b4 = x[4] * f_ring * f_ring * f_ring * f_ring;
        // IMRPhenomXPNR co-precessing deviations (zero otherwise)
        let b1 = fma(wf.coprec.dev, wf.coprec.zeta2, b1);
        let b4 = fma(wf.coprec.dev, wf.coprec.zeta1, b4);

        let c_l_gr = c_l;

        Self {
            f_phase_ins_min,
            f_phase_ins_max,
            f_phase_rd_min,
            f_phase_rd_max,
            f_phase_match_in,
            f_phase_match_im,
            c1_int: 0.0,
            c2_int: 0.0,
            c1_mrd: 0.0,
            c2_mrd: 0.0,
            c0,
            c1,
            c2,
            c4,
            c_rd,
            c_l,
            c_l_gr,
            c4_ov3: c4 / 3.0,
            c_l_ov_fda: c_l / f_damp,
            b0,
            b1,
            b2,
            b3,
            b4,
            a0,
            a1,
            a2,
            a3,
            a4,
            sigma1,
            sigma2,
            sigma3,
            sigma4,
            sigma5,
            phi0,
            phi1,
            phi2,
            phi3,
            phi4,
            phi5,
            phi5l,
            phi6,
            phi6l,
            phi7,
            phi8,
            phi8l,
            phi9,
            phi9l,
            phi_initial: -PI / 4.0,
            phi_norm,
            dphi0,
            dphi1,
            dphi2,
            dphi3,
            dphi4,
            dphi5,
            dphi6,
            dphi6l,
            dphi7,
            dphi8,
            dphi8l,
            dphi9,
            dphi9l,
            collocation_points_phase_rd: points_rd,
            collocation_values_phase_rd: values_rd,
            collocation_points_phase_int: points_int,
            collocation_values_phase_int: values_int,
            collocation_points_phase_ins: points_ins,
            collocation_values_phase_ins: values_ins,
            f_ring,
            f_damp,
            f_damp2: f_damp * f_damp,
            inv_f_damp: 1.0 / f_damp,
            inv_2_f_damp: 1.0 / (2.0 * f_damp),
            c1_3ov2: 1.5 * c1,
            b3_ov2: b3 / 2.0,
            b4_ov3: b4 / 3.0,
            int_lorentzian_numerator: 4.0 * c_l_gr,
            int_lorentzian_width2: 4.0 * f_damp * f_damp,
            int_atan_prefactor: 2.0 * c_l_gr / f_damp,
            dphase_norm: inspiral::inspiral_dphase_norm(),
        }
    }

    /// Phase (without the 1/eta scaling and linear terms) using the region
    /// boundaries of `IMRPhenomXASGenerateFD`: `mf < fPhaseMatchIN` is
    /// inspiral, `mf > fPhaseMatchIM` ringdown, intermediate otherwise.
    /// Long inspirals have phases of order 1e6 rad, so this is evaluated
    /// in nested form with fused multiply-adds.
    #[inline]
    pub(crate) fn phase_xas(&self, mf: f64, powers: &UsefulPowers) -> f64 {
        if mf < self.f_phase_match_in {
            inspiral::phase_22_ansatz_int(mf, powers, self)
        } else if mf > self.f_phase_match_im {
            let phi = ringdown::phase_22_ansatz_int(mf, powers, self);
            fma(self.c2_mrd, mf, phi + self.c1_mrd)
        } else {
            let phi = intermediate::phase_22_ansatz_int(mf, powers, self);
            fma(self.c2_int, mf, phi + self.c1_int)
        }
    }
}

/// `IMRPhenomX_Phase_22_ConnectionCoefficients`: C1 and C2 coefficients
/// making the phase and its derivative continuous at the transitions.
pub(crate) fn phase_22_connection_coefficients(phase: &mut PhaseCoefficients) {
    let f_ins = phase.f_phase_match_in;
    let f_int = phase.f_phase_match_im;

    // phi_Inspiral(f) = phi_Intermediate(f) + C1 + C2 f at f_ins
    let powers_of_f_ins = UsefulPowers::new(f_ins);
    let dphi_ins = inspiral::phase_22_ansatz(f_ins, &powers_of_f_ins, phase);
    let dphi_int = intermediate::phase_22_ansatz(f_ins, &powers_of_f_ins, phase);
    phase.c2_int = dphi_ins - dphi_int;
    let phi_in = inspiral::phase_22_ansatz_int(f_ins, &powers_of_f_ins, phase);
    let phi_im = intermediate::phase_22_ansatz_int(f_ins, &powers_of_f_ins, phase);
    phase.c1_int = fma(-phase.c2_int, f_ins, phi_in - phi_im);

    // phi_Intermediate(f) = phi_Ringdown(f) + C1 + C2 f at f_int
    let powers_of_f_int = UsefulPowers::new(f_int);
    let phi_imc = fma(
        phase.c2_int,
        f_int,
        intermediate::phase_22_ansatz_int(f_int, &powers_of_f_int, phase) + phase.c1_int,
    );
    let phi_rd = ringdown::phase_22_ansatz_int(f_int, &powers_of_f_int, phase);
    let dphi_intc = intermediate::phase_22_ansatz(f_int, &powers_of_f_int, phase) + phase.c2_int;
    let dphi_rd = ringdown::phase_22_ansatz(f_int, &powers_of_f_int, phase);
    phase.c2_mrd = dphi_intc - dphi_rd;
    phase.c1_mrd = fma(-phase.c2_mrd, f_int, phi_imc - phi_rd);
}

/// `IMRPhenomX_Phase_22`: the phase before scaling by 1/eta, with the
/// region boundaries of LAL (`IMRPhenomX_StepFuncBool`, i.e. `mf >= f`).
pub(crate) fn phase_22(mf: f64, powers: &UsefulPowers, phase: &PhaseCoefficients) -> f64 {
    if mf < phase.f_phase_match_in {
        inspiral::phase_22_ansatz_int(mf, powers, phase)
    } else if mf >= phase.f_phase_match_im {
        fma(
            phase.c2_mrd,
            mf,
            ringdown::phase_22_ansatz_int(mf, powers, phase) + phase.c1_mrd,
        )
    } else {
        fma(
            phase.c2_int,
            mf,
            intermediate::phase_22_ansatz_int(mf, powers, phase) + phase.c1_int,
        )
    }
}

/// `IMRPhenomX_dPhase_22`: the phase derivative before scaling by 1/eta.
pub(crate) fn dphase_22(mf: f64, powers: &UsefulPowers, phase: &PhaseCoefficients) -> f64 {
    if mf < phase.f_phase_match_in {
        inspiral::phase_22_ansatz(mf, powers, phase)
    } else if mf >= phase.f_phase_match_im {
        ringdown::phase_22_ansatz(mf, powers, phase) + phase.c2_mrd
    } else {
        intermediate::phase_22_ansatz(mf, powers, phase) + phase.c2_int
    }
}

/// The phase of the 22 ([`phase_22`]) with its first ([`dphase_22`]) and
/// second derivatives, the latter by differentiating the phase-derivative
/// ansatz of the region of `mf`.
pub(crate) fn phase_22_jet(mf: f64, powers: &UsefulPowers, phase: &PhaseCoefficients) -> Jet {
    let x = Jet::variable(mf);
    let p = phase;
    let dphase = if mf < p.f_phase_match_in {
        let pw = |k: f64| Jet::power(mf, k);
        let log = x.ln();
        let pn = pw(2. / 3.) * p.dphi2
            + x * p.dphi3
            + pw(4. / 3.) * p.dphi4
            + pw(5. / 3.) * p.dphi5
            + pw(2.) * p.dphi6
            + pw(2.) * log * p.dphi6l
            + pw(7. / 3.) * p.dphi7
            + pw(8. / 3.) * p.dphi8
            + pw(8. / 3.) * log * p.dphi8l
            + p.dphi0;
        let pseudo = pw(8. / 3.) * p.a0 + pw(3.) * p.a1 + pw(10. / 3.) * p.a2 + pw(11. / 3.) * p.a3;
        (pn + pseudo) * pw(-8. / 3.) * p.dphase_norm
    } else if mf >= p.f_phase_match_im {
        let dfr = x - p.f_ring;
        Jet::power(mf, -1. / 3.) * p.c1
            + Jet::power(mf, -2.) * p.c2
            + Jet::power(mf, -4.) * p.c4
            + p.c_l / (dfr * dfr + p.f_damp * p.f_damp)
            + p.c0
    } else {
        let dfr = x - p.f_ring;
        x.recip() * p.b1
            + Jet::power(mf, -2.) * p.b2
            + Jet::power(mf, -3.) * p.b3
            + Jet::power(mf, -4.) * p.b4
            + (4. * p.c_l_gr) / (dfr * dfr + 4. * p.f_damp * p.f_damp)
            + p.b0
    };
    Jet::new(
        phase_22(mf, powers, phase),
        dphase_22(mf, powers, phase),
        dphase.d1,
    )
}

/// `IMRPhenomX_TimeShift_22`: linear-in-frequency phase coefficient `linb`
/// aligning the peak of the strain near t = 0.
pub(crate) fn time_shift_22(phase: &PhaseCoefficients, wf: &WaveformStruct) -> f64 {
    // linb is a parameter-space fit of dphi22(fring22 - fdamp22).
    let linb = linb(wf.eta, wf.s_tot_r, wf.dchi, wf.delta);
    let f_ref_fit = wf.f_ring - wf.f_damp;
    let powers_of_f_ref_fit = UsefulPowers::new(f_ref_fit);
    let dphi22_ref = 1.0 / wf.eta * dphase_22(f_ref_fit, &powers_of_f_ref_fit, phase);
    // Align the peak of psi4, then correct to align the peak of the strain.
    let psi4_to_strain = psi4_to_strain(wf.eta, wf.s_tot_r, wf.dchi);
    fma(-(2. * PI), 500. + psi4_to_strain, linb - dphi22_ref)
}
