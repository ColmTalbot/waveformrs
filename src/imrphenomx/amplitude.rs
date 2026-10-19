// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! Amplitude of the (2,2) mode: `IMRPhenomXAmpCoefficients`,
//! `IMRPhenomXGetAmplitudeCoefficients` and `IMRPhenomX_Amplitude_22` from
//! `LALSimIMRPhenomX_internals.c`.

use std::f64::consts::PI;

use crate::imrphenomx::internals::{fma, UsefulPowers, WaveformStruct, POWERS_OF_LALPI};
use crate::imrphenomx::{inspiral, intermediate, ringdown};

/// Equivalent of `IMRPhenomXAmpCoefficients` for the default versions
/// (InsAmpVersion 103, IntAmpVersion 104, RDAmpVersion 103).
#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub(crate) struct AmpCoefficients {
    // Transition frequencies (dimensionless)
    pub(crate) f_amp_ins_min: f64,
    pub(crate) f_amp_ins_max: f64,
    /// Inspiral-intermediate transition.
    pub(crate) f_amp_match_in: f64,
    pub(crate) f_amp_int_max: f64,
    /// Peak of the ringdown amplitude, intermediate-ringdown transition.
    pub(crate) f_amp_rd_min: f64,

    // TaylorF2 PN amplitude coefficients (including powers of pi)
    pub(crate) pn_initial: f64,
    pub(crate) pn_one_third: f64,
    pub(crate) pn_two_thirds: f64,
    pub(crate) pn_three_thirds: f64,
    pub(crate) pn_four_thirds: f64,
    pub(crate) pn_five_thirds: f64,
    pub(crate) pn_six_thirds: f64,
    pub(crate) pn_seven_thirds: f64,
    pub(crate) pn_eight_thirds: f64,
    pub(crate) pn_nine_thirds: f64,

    /// Inspiral pseudo-PN coefficients.
    pub(crate) rho1: f64,
    pub(crate) rho2: f64,
    pub(crate) rho3: f64,

    /// Intermediate polynomial coefficients.
    pub(crate) delta0: f64,
    pub(crate) delta1: f64,
    pub(crate) delta2: f64,
    pub(crate) delta3: f64,
    pub(crate) delta4: f64,
    pub(crate) delta5: f64,

    /// Ringdown coefficients.
    pub(crate) v1_rd: f64,
    pub(crate) gamma1: f64,
    pub(crate) gamma2: f64,
    pub(crate) gamma3: f64,
    /// gamma2 / (f_damp gamma3)
    pub(crate) gamma_r: f64,
    /// (gamma3 f_damp)^2
    pub(crate) gamma_d2: f64,
    /// f_damp gamma1 gamma3
    pub(crate) gamma_d13: f64,

    pub(crate) collocation_points_amp_ins: [f64; 4],
    pub(crate) collocation_values_amp_ins: [f64; 4],
    pub(crate) collocation_points_amp_int: [f64; 4],
    pub(crate) collocation_values_amp_int: [f64; 4],
    /// Derivatives of the inverse amplitude at the intermediate boundaries.
    pub(crate) d1: f64,
    pub(crate) d4: f64,

    /// Ringdown frequency copied from the waveform struct.
    pub(crate) f_ring: f64,
    /// Damping frequency copied from the waveform struct.
    pub(crate) f_damp: f64,
}

impl AmpCoefficients {
    /// `IMRPhenomXGetAmplitudeCoefficients`.
    pub(crate) fn new(wf: &WaveformStruct) -> Self {
        let pi = &*POWERS_OF_LALPI;
        let chi1l = wf.chi1l;
        let chi2l = wf.chi2l;
        let chi1l2 = wf.chi1l2;
        let chi1l3 = wf.chi1l3;
        let chi2l2 = wf.chi2l2;
        let delta = wf.delta;
        let eta = wf.eta;
        let eta2 = wf.eta2;
        let eta3 = wf.eta3;
        let f_ring = wf.f_ring;
        let f_damp = wf.f_damp;

        // Phenomenological ringdown coefficients: gamma2 = lambda and
        // gamma3 = sigma in arXiv:2001.11412.
        let gamma2 = ringdown::amp_22_gamma2(eta, wf.s_tot_r, wf.dchi, wf.delta);
        let gamma3 = ringdown::amp_22_gamma3(eta, wf.s_tot_r, wf.dchi, wf.delta);
        // IMRPhenomXPNR co-precessing deviation (zero otherwise)
        let gamma3 = fma(wf.coprec.dev, wf.coprec.mu2, gamma3);

        // Peak ringdown frequency, Eq. 5.14 in arXiv:2001.11412.
        let f_amp_rd_min = ringdown::amp_22_peak_frequency(gamma2, gamma3, f_ring, f_damp);

        // Value of v1RD at fAmpRDMin is used to calculate gamma1.
        let v1_rd = ringdown::amp_22_v1(eta, wf.s_tot_r, wf.dchi, wf.delta);
        let v1_rd = fma(wf.coprec.dev, wf.coprec.mu1, v1_rd);
        let f1 = f_amp_rd_min;

        // Solving linear equation: ansatzRD(F1) == v1
        let gamma1 = (v1_rd / (f_damp * gamma3))
            * fma(
                f_damp * f_damp * gamma3,
                gamma3,
                fma(f_ring, f_ring, fma(f1, f1, -(2.0 * f1 * f_ring))),
            )
            * (((f1 - f_ring) * gamma2) / (f_damp * gamma3)).exp();

        let gamma_r = gamma2 / (f_damp * gamma3);
        let gamma_d2 = (gamma3 * f_damp) * (gamma3 * f_damp);
        let gamma_d13 = f_damp * gamma1 * gamma3;

        // Inspiral-intermediate transition frequencies, Eq. 5.7 in arXiv:2001.11412
        let f_amp_ins_min = 0.0026;
        let f_amp_ins_max = wf.f_meco + 0.25 * (wf.f_isco - wf.f_meco);
        let f_amp_match_in = f_amp_ins_max;
        let f_amp_int_max = f_amp_rd_min;

        // TaylorF2 PN amplitude coefficients from the PN re-expansion of the Hlm.
        let pn_initial = 1.0;
        let pn_one_third = 0.0;
        let pn_two_thirds = fma(1804.0, eta, -969.0) / 672.0 * pi.two_thirds;
        let pn_three_thirds = fma(
            -(44.0 * (chi1l + chi2l)),
            eta,
            fma(
                -(81.0 * chi2l),
                delta,
                fma(81.0, chi1l + chi2l, 81.0 * chi1l * delta),
            ),
        ) / 48.0
            * pi.itself;
        let pn_four_thirds = fma(
            24.0,
            fma(
                1473794.0,
                eta2,
                fma(
                    428652.0 * chi2l2,
                    -1.0 + delta,
                    fma(
                        10584.0,
                        fma(81.0, chi2l2, fma(81.0, chi1l2, -(94.0 * chi1l * chi2l))),
                        -1975055.0,
                    ) * eta,
                ),
            ),
            fma(-(10287648.0 * chi1l2), 1.0 + delta, -27312085.0),
        ) / 8.128512e6
            * pi.one_third
            * pi.itself;
        let pn_five_thirds = fma(
            42840.0 * fma(4.0, eta, -1.0),
            PI,
            fma(
                chi1l,
                fma(
                    287213.0,
                    1.0 + delta,
                    -(4.0 * eta * fma(8908.0, eta, fma(2083.0, delta, 93414.0))),
                ),
                fma(
                    -6048.0 * chi1l3,
                    fma(3.0 + delta, eta, -1.0 - delta),
                    chi2l
                        * fma(
                            -35632.0,
                            eta2,
                            fma(
                                -fma(6048.0, chi2l2, 287213.0),
                                -1.0 + delta,
                                4.0 * fma(
                                    2083.0,
                                    delta,
                                    fma(1512.0 * chi2l2, -3.0 + delta, -93414.0),
                                ) * eta,
                            ),
                        ),
                ),
            ),
        ) / 32256.0
            * pi.five_thirds;
        let pn_six_thirds = fma(
            12.0,
            fma(
                eta,
                fma(
                    -77616.0,
                    fma(
                        275520.0,
                        pi.two,
                        fma(
                            -(107520.0 * (chi1l + chi2l)),
                            PI,
                            fma(
                                -chi2l2,
                                fma(282718.0, delta, 158633.0),
                                fma(
                                    638642.0 * chi1l,
                                    chi2l,
                                    chi1l2 * fma(282718.0, delta, -158633.0),
                                ),
                            ),
                        ),
                    ),
                    545384828789.0,
                ),
                fma(
                    -116424.0,
                    fma(
                        147117.0,
                        fma(-chi2l2, -1.0 + delta, chi1l2 * (1.0 + delta)),
                        60928.0 * fma(-chi2l, delta, fma(chi1l, delta, chi1l + chi2l)) * PI,
                    ),
                    fma(
                        28.0 * fma(
                            11088.0,
                            fma(
                                163199.0,
                                chi2l2,
                                fma(163199.0, chi1l2, -(266498.0 * chi1l * chi2l)),
                            ),
                            -3248849057.0,
                        ),
                        eta2,
                        27026893936.0 * eta3,
                    ),
                ),
            ),
            -1242641879927.0,
        ) / 6.0085960704e10
            * pi.two;

        // Pseudo-PN coefficients from the collocation points (version 103).
        let f_ins = [
            0.25 * f_amp_match_in,
            0.50 * f_amp_match_in,
            0.75 * f_amp_match_in,
            1.00 * f_amp_match_in,
        ];
        let v_ins = [
            0.0,
            inspiral::amp_22_v2(eta, wf.chi_pn_hat, wf.dchi, wf.delta),
            inspiral::amp_22_v3(eta, wf.chi_pn_hat, wf.dchi, wf.delta),
            inspiral::amp_22_v4(eta, wf.chi_pn_hat, wf.dchi, wf.delta),
        ];
        let [rho1, rho2, rho3] = inspiral::amp_22_rho(
            [v_ins[1], v_ins[2], v_ins[3]],
            [f_ins[1], f_ins[2], f_ins[3]],
        );

        let mut amp = Self {
            f_amp_ins_min,
            f_amp_ins_max,
            f_amp_match_in,
            f_amp_int_max,
            f_amp_rd_min,
            pn_initial,
            pn_one_third,
            pn_two_thirds,
            pn_three_thirds,
            pn_four_thirds,
            pn_five_thirds,
            pn_six_thirds,
            pn_seven_thirds: rho1,
            pn_eight_thirds: rho2,
            pn_nine_thirds: rho3,
            rho1,
            rho2,
            rho3,
            delta0: 0.0,
            delta1: 0.0,
            delta2: 0.0,
            delta3: 0.0,
            delta4: 0.0,
            delta5: 0.0,
            v1_rd,
            gamma1,
            gamma2,
            gamma3,
            gamma_r,
            gamma_d2,
            gamma_d13,
            collocation_points_amp_ins: f_ins,
            collocation_values_amp_ins: v_ins,
            collocation_points_amp_int: [0.0; 4],
            collocation_values_amp_int: [0.0; 4],
            d1: 0.0,
            d4: 0.0,
            f_ring,
            f_damp,
        };

        // Intermediate region from the derivatives of the inspiral and
        // ringdown fits at the boundaries.
        let f1 = amp.f_amp_match_in;
        let f4 = amp.f_amp_rd_min;
        let powers_of_f1 = UsefulPowers::new(f1);
        let powers_of_f4 = UsefulPowers::new(f4);

        let d1 = inspiral::amp_22_dansatz(f1, wf, &amp);
        let d4 = ringdown::amp_22_dansatz(f4, f_ring, f_damp, &amp);
        let insp_f1 = inspiral::amp_22_ansatz_exact(f1, &powers_of_f1, &amp);
        let rd_f4 = ringdown::amp_22_ansatz_exact(f4, f_ring, f_damp, &amp);

        // D[f^{7/6} Ah^{-1}[f], f] = (7/6) f^{1/6} / Ah[f] - f^{7/6} Ah'[f] / Ah[f]^2
        let d1 = ((7.0 / 6.0) * powers_of_f1.one_sixth / insp_f1)
            - (powers_of_f1.seven_sixths * d1 / (insp_f1 * insp_f1));
        let d4 = ((7.0 / 6.0) * powers_of_f4.one_sixth / rd_f4)
            - (powers_of_f4.seven_sixths * d4 / (rd_f4 * rd_f4));

        // Fourth order polynomial (version 104).
        let f2 = f1 + (1.0 / 2.0) * (f4 - f1);
        let f3 = 0.0;
        let v1 = 1.0 / (powers_of_f1.m_seven_sixths * insp_f1);
        let v2 = 1.0 / intermediate::amp_22_va(eta, wf.s_tot_r, wf.dchi, wf.delta);
        let v2 = fma(wf.coprec.dev, wf.coprec.mu3, v2);
        let v3 = 0.0;
        let v4 = 1.0 / (powers_of_f4.m_seven_sixths * rd_f4);

        let [delta0, delta1, delta2, delta3, delta4, delta5] =
            intermediate::amp_22_deltas(d1, d4, v1, v2, v4, f1, f2, f4);
        amp.delta0 = delta0;
        amp.delta1 = delta1;
        amp.delta2 = delta2;
        amp.delta3 = delta3;
        amp.delta4 = delta4;
        amp.delta5 = delta5;
        amp.collocation_points_amp_int = [f1, f2, f3, f4];
        amp.collocation_values_amp_int = [v1, v2, v3, v4];
        amp.d1 = d1;
        amp.d4 = d4;
        amp
    }

    /// The amplitude ansatz of the region containing `mf`, without the
    /// f^(-7/6) prefactor, using the region boundaries of
    /// `IMRPhenomXASGenerateFD` (`mf < fAmpMatchIN` is inspiral,
    /// `mf > fAmpRDMin` ringdown).
    #[inline]
    pub(crate) fn ansatz(&self, mf: f64, powers: &UsefulPowers) -> f64 {
        if mf < self.f_amp_match_in {
            inspiral::amp_22_ansatz(mf, powers, self)
        } else if mf > self.f_amp_rd_min {
            ringdown::amp_22_ansatz(mf, self)
        } else {
            intermediate::amp_22_ansatz(mf, powers, self)
        }
    }
}

/// `IMRPhenomX_Amplitude_22`: the (2,2) amplitude including the
/// `ampNorm f^(-7/6)` prefactor (multiply by `amp0` for strain).
#[allow(dead_code)]
pub(crate) fn amplitude_22(
    mf: f64,
    powers: &UsefulPowers,
    amp: &AmpCoefficients,
    wf: &WaveformStruct,
) -> f64 {
    let prefactor = wf.amp_norm / powers.seven_sixths;
    if mf < amp.f_amp_match_in {
        prefactor * inspiral::amp_22_ansatz(mf, powers, amp)
    } else if mf >= amp.f_amp_rd_min {
        prefactor * ringdown::amp_22_ansatz(mf, amp)
    } else {
        prefactor * intermediate::amp_22_ansatz(mf, powers, amp)
    }
}
