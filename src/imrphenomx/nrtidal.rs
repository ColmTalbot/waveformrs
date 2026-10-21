// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! The NRTidal corrections of IMRPhenomXAS_NRTidalv2/v3 and
//! IMRPhenomXP_NRTidalv2/v3: `IMRPhenomXGetTidalPhaseCoefficients`,
//! `IMRPhenomX_TidalPhase` and `IMRPhenomX_TidalPhaseDerivative` in
//! `LALSimIMRPhenomX_internals.c`, and the tidal terms of
//! `IMRPhenomXASGenerateFD` and `IMRPhenomXPGenerateFD` in
//! `LALSimIMRPhenomX.c`.

use std::f64::consts::PI;

use crate::constants::MTSUN_SI;
use crate::imrphenomx::internals::{UsefulPowers, WaveformStruct, POWERS_OF_LALPI};
use crate::nrtidal::{self, NRTidalVersion, TidalSeries, NRTIDALV2_COEFFS};

/// `XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff`.
fn qm2so_4pn(m_by_m: f64) -> f64 {
    -720. / 9.6 * m_by_m * m_by_m
}

/// `XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff`.
fn qm2s_4pn(m_by_m: f64) -> f64 {
    240. / 9.6 * m_by_m * m_by_m
}

/// `XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff`.
fn qm2s_6pn(m_by_m: f64) -> f64 {
    (4703.5 / 8.4 + 2935. / 6. * m_by_m - 120. * m_by_m * m_by_m) * m_by_m * m_by_m
}

/// The planck taper of the NRTidalv3 phase in `IMRPhenomX_TidalPhase`
/// (Eq. 46 of arXiv:2311.07456) and, for the derivative, its derivative.
fn v3_taper(mf: f64, mf_merger: f64) -> (f64, f64) {
    let start = 1.15 * mf_merger;
    let end = 1.35 * mf_merger;
    if mf <= start {
        (0.0, 0.0)
    } else if mf >= end {
        (1.0, 0.0)
    } else {
        let exp_arg = (end - start) / (mf - start) + (end - start) / (mf - end);
        let exp_fn = exp_arg.cosh() + exp_arg.sinh();
        let taper = 1. / (exp_fn + 1.0);
        let dfactor = -exp_fn
            * (-(end - start) / ((mf - start) * (mf - start))
                - (end - start) / ((mf - end) * (mf - end)));
        (taper, dfactor * (taper * taper))
    }
}

/// The NRTidal state of the (2,2) mode: the parameters of the waveform
/// struct with `PhenomXTidalFlag` 2 or 3 and the coefficients of
/// `IMRPhenomXGetTidalPhaseCoefficients`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Tidal {
    pub(crate) version: NRTidalVersion,
    m1_si: f64,
    m2_si: f64,
    m_tot: f64,
    lambda1: f64,
    lambda2: f64,
    chi1l: f64,
    chi2l: f64,
    x_a: f64,
    x_b: f64,
    kappa2t: f64,
    /// `3 / (128 X_A X_B)`.
    pfa_n: f64,
    /// `c2PN_tidal`, `c3PN_tidal`, `c3p5PN_tidal`.
    c2pn: f64,
    c3pn: f64,
    c3p5pn: f64,
    /// The merger frequency (Hz) of the version.
    f_merger: f64,
    /// NRTidalv3: the 7.5PN coefficients and those of the phase.
    pn: [f64; 10],
    v3: [f64; 20],
}

impl Tidal {
    /// The tidal state of a waveform struct with NRTidal corrections.
    pub(crate) fn new(wf: &WaveformStruct, version: NRTidalVersion) -> Self {
        let (x_a, x_b) = (wf.m1, wf.m2);
        let (quadparam1, quadparam2) = (wf.quadparam1, wf.quadparam2);
        let c2pn = (qm2so_4pn(x_a) + qm2s_4pn(x_a)) * (quadparam1 - 1.) * wf.chi1l2
            + (qm2so_4pn(x_b) + qm2s_4pn(x_b)) * (quadparam2 - 1.) * wf.chi2l2;
        let c3pn = qm2s_6pn(x_a) * (quadparam1 - 1.) * wf.chi1l2
            + qm2s_6pn(x_b) * (quadparam2 - 1.) * wf.chi2l2;
        let (ss, sss) =
            nrtidal::ho_spin_terms(x_a, x_b, wf.chi1l, wf.chi2l, quadparam1, quadparam2);
        let f_merger = match version {
            NRTidalVersion::V2 => nrtidal::merger_frequency(wf.m_tot, wf.kappa2t, wf.q),
            NRTidalVersion::V3 => nrtidal::merger_frequency_v3(
                wf.m_tot, wf.lambda1, wf.lambda2, wf.q, wf.chi1l, wf.chi2l,
            ),
        };
        let pn = nrtidal::pn_coeffs(x_a);
        let v3 = match version {
            NRTidalVersion::V2 => [0.0; 20],
            NRTidalVersion::V3 => nrtidal::v3_coeffs(x_a, wf.m_tot, wf.lambda1, wf.lambda2, &pn),
        };
        Self {
            version,
            m1_si: wf.m1_si,
            m2_si: wf.m2_si,
            m_tot: wf.m_tot,
            lambda1: wf.lambda1,
            lambda2: wf.lambda2,
            chi1l: wf.chi1l,
            chi2l: wf.chi2l,
            x_a,
            x_b,
            kappa2t: wf.kappa2t,
            pfa_n: 3. / (128. * x_a * x_b),
            c2pn,
            c3pn,
            c3p5pn: ss + sss,
            f_merger,
            pn,
            v3,
        }
    }

    /// The merger frequency (Hz) of the version, where the time shift is
    /// set (unless the grid ends below it).
    pub(crate) fn f_merger(&self) -> f64 {
        self.f_merger
    }

    /// `XLALSimNRTunedTidesFDTidalPhaseFrequencySeries` at `frequencies`
    /// (Hz) for this binary.
    pub(crate) fn series(&self, frequencies: &[f64]) -> TidalSeries {
        nrtidal::tidal_series(
            frequencies,
            self.m1_si,
            self.m2_si,
            self.lambda1,
            self.lambda2,
            self.chi1l,
            self.chi2l,
            self.version,
        )
    }

    /// The spin-induced quadrupole (and, for 3.5PN, octupole) phase that
    /// the main loops add to the tidal phase of the series.
    #[inline]
    pub(crate) fn multipole_phase(&self, powers: &UsefulPowers) -> f64 {
        let pi = &*POWERS_OF_LALPI;
        let mut phase = self.pfa_n * self.c2pn * pi.m_one_third * powers.m_one_third;
        phase += self.pfa_n * self.c3pn * pi.one_third * powers.one_third;
        phase += self.pfa_n * self.c3p5pn * pi.two_thirds * powers.two_thirds;
        phase
    }

    /// The (2,2) mode with the tidal corrections, `amp0 (ampNorm Mf^(-7/6)
    /// A(Mf) + 2 sqrt(pi / 5) A_T) exp(i (phi - phi_T)) window`, for the
    /// IMRPhenomXAS amplitude ansatz `amp` and phase `phi` at `mf` and
    /// the point `index` of `series`.
    #[inline]
    pub(crate) fn h22(
        &self,
        wf: &WaveformStruct,
        powers: &UsefulPowers,
        amp: f64,
        phi: f64,
        series: &TidalSeries,
        index: usize,
    ) -> num_complex::Complex<f64> {
        let pi = &*POWERS_OF_LALPI;
        let phase_tidal = series.phase[index] + self.multipole_phase(powers);
        let amplitude = wf.amp0
            * (wf.amp_norm * powers.m_seven_sixths * amp
                + 2. * (1. / 5.0f64).sqrt() * pi.sqrt * series.amplitude[index]);
        let (sin, cos) = (phi - phase_tidal).sin_cos();
        num_complex::Complex::new(amplitude * cos, amplitude * sin) * series.taper[index]
    }

    /// The amplitude and phase of [`Tidal::h22`], `h = A exp(i phi)`.
    #[inline]
    pub(crate) fn amplitude_phase(
        &self,
        wf: &WaveformStruct,
        powers: &UsefulPowers,
        amp: f64,
        phi: f64,
        series: &TidalSeries,
        index: usize,
    ) -> (f64, f64) {
        let pi = &*POWERS_OF_LALPI;
        let phase_tidal = series.phase[index] + self.multipole_phase(powers);
        let amplitude = wf.amp0
            * (wf.amp_norm * powers.m_seven_sixths * amp
                + 2. * (1. / 5.0f64).sqrt() * pi.sqrt * series.amplitude[index]);
        (amplitude * series.taper[index], phi - phase_tidal)
    }

    /// `IMRPhenomX_TidalPhase` at `mf`.
    pub(crate) fn phase(&self, mf: f64, powers: &UsefulPowers) -> f64 {
        let pi = &*POWERS_OF_LALPI;
        let (x_a, x_b) = (self.x_a, self.x_b);
        let pfa_n = self.pfa_n;
        let mut phase_tidal = 0.;
        phase_tidal += pfa_n * self.c2pn * pi.m_one_third * powers.m_one_third;
        phase_tidal += pfa_n * self.c3pn * pi.one_third * powers.one_third;
        let nr_phase = match self.version {
            NRTidalVersion::V2 => {
                let [c_newt, n1, n3over2, n2, n5over2, n3, d1, d3over2, d2] = NRTIDALV2_COEFFS;
                -((c_newt
                    * self.kappa2t
                    * powers.five_thirds
                    * pi.five_thirds
                    * (1.
                        + powers.two_thirds * n1 * pi.two_thirds
                        + mf * n3over2 * PI
                        + powers.four_thirds * n2 * pi.four_thirds
                        + powers.five_thirds * n5over2 * pi.five_thirds
                        + mf.powf(2.) * n3 * pi.two))
                    / ((1.
                        + d1 * powers.two_thirds * pi.two_thirds
                        + d3over2 * mf * PI
                        + d2 * powers.four_thirds * pi.four_thirds)
                        * x_a
                        * x_b))
            }
            NRTidalVersion::V3 => {
                let terms = self.v3_terms(mf, powers);
                let mf_merger = self.f_merger * (self.m_tot * MTSUN_SI);
                let (taper, _) = v3_taper(mf, mf_merger);
                // Eq. (45) in arXiv:2311.07456
                terms.nr * (1. - taper) + terms.pn * taper
            }
        };
        phase_tidal += nr_phase;
        // the PN SS-tail and SSS terms
        phase_tidal += pfa_n * self.c3p5pn * pi.two_thirds * powers.two_thirds;
        phase_tidal
    }

    /// `IMRPhenomX_TidalPhaseDerivative` at `mf`.
    pub(crate) fn phase_derivative(&self, mf: f64, powers: &UsefulPowers) -> f64 {
        let pi = &*POWERS_OF_LALPI;
        let (x_a, x_b) = (self.x_a, self.x_b);
        let mf_ten_thirds = powers.eight_thirds * powers.two_thirds;
        let pi_ten_thirds = pi.eight_thirds * pi.two_thirds;
        let pfa_n = self.pfa_n;
        let (c2pn, c3pn, c3p5pn) = (self.c2pn, self.c3pn, self.c3p5pn);
        let three_pn = (pfa_n * (-c2pn + c3pn * powers.two_thirds * pi.two_thirds))
            / (3. * powers.four_thirds * pi.one_third);
        let mut dphase = three_pn;
        dphase += (2. * c3p5pn * pfa_n * pi.two_thirds) / (3. * powers.one_third);
        let nr_dphase = match self.version {
            NRTidalVersion::V2 => {
                let [c_newt, n1, n3over2, n2, n5over2, n3, d1, d3over2, d2] = NRTIDALV2_COEFFS;
                (c_newt
                    * self.kappa2t
                    * powers.two_thirds
                    * pi.five_thirds
                    * (-5.
                        - powers.two_thirds * (3. * d1 + 7. * n1) * pi.two_thirds
                        - 2. * mf * (d3over2 + 4. * n3over2) * PI
                        - powers.four_thirds * (d2 + 5. * d1 * n1 + 9. * n2) * pi.four_thirds
                        - 2. * powers.five_thirds
                            * (2. * d3over2 * n1 + 3. * d1 * n3over2 + 5. * n5over2)
                            * pi.five_thirds
                        - mf.powf(2.)
                            * (3. * d2 * n1 + 7. * d1 * n2 + 11. * n3 + 5. * d3over2 * n3over2)
                            * pi.two
                        - 2. * powers.seven_thirds
                            * (3. * d3over2 * n2 + 2. * d2 * n3over2 + 4. * d1 * n5over2)
                            * pi.seven_thirds
                        - powers.eight_thirds
                            * (5. * d2 * n2 + 9. * d1 * n3 + 7. * d3over2 * n5over2)
                            * pi.eight_thirds
                        - 2. * mf.powf(3.) * (4. * d3over2 * n3 + 3. * d2 * n5over2) * pi.three
                        - 7. * d2 * mf_ten_thirds * n3 * pi_ten_thirds))
                    / (3.
                        * (1.
                            + d1 * powers.two_thirds * pi.two_thirds
                            + d3over2 * mf * PI
                            + d2 * powers.four_thirds * pi.four_thirds)
                            .powf(2.)
                        * x_a
                        * x_b)
            }
            NRTidalVersion::V3 => {
                let terms = self.v3_terms(mf, powers);
                let mf_merger = self.f_merger * (self.m_tot * MTSUN_SI);
                let (taper, dtaper) = v3_taper(mf, mf_merger);
                terms.dnr * (1.0 - taper) + terms.dpn * taper - (terms.nr - terms.pn) * dtaper
            }
        };
        dphase += nr_dphase;
        dphase
    }

    /// The NRTidalv3 phase (`NRphaseNRT`), the 7.5PN tidal phase
    /// (`PNtidalphase`) and their derivatives, as written out in
    /// `IMRPhenomX_TidalPhase` and `IMRPhenomX_TidalPhaseDerivative`.
    fn v3_terms(&self, mf: f64, powers: &UsefulPowers) -> V3Terms {
        let pi = &*POWERS_OF_LALPI;
        let c = &self.v3;
        let pn = &self.pn;
        let (s1, s2, s3) = (c[0], c[1], c[2]);
        let s2s3 = s2 * s3;
        let s2mf = -s2 * mf * PI * 2.0;
        let exps2s3 = s2s3.cosh() + s2s3.sinh();
        let exps2mf = s2mf.cosh() + s2mf.sinh();
        // Eq. (27) of arXiv:2311.07456
        let dynk2bar = 1.0 + (s1 - 1.) * (1.0 / (1.0 + exps2mf * exps2s3))
            - ((s1 - 1.0) / (1.0 + exps2s3))
            - 2.0 * (mf * PI) * (s1 - 1.) * s2 * exps2s3 / ((1.0 + exps2s3) * (1.0 + exps2s3));
        let (kappa_a, kappa_b) = (c[4], c[5]);

        // per star: (c_Newt, c_1, c_3/2, c_2, c_5/2), the Padé n_5/2, n_3,
        // d_1 and the constrained n_1, n_3/2, n_2, d_3/2
        let stars = [
            (
                kappa_a,
                &pn[0..5],
                c[6],
                c[7],
                c[8],
                c[12],
                c[13],
                c[14],
                c[15],
            ),
            (
                kappa_b,
                &pn[5..10],
                c[9],
                c[10],
                c[11],
                c[16],
                c[17],
                c[18],
                c[19],
            ),
        ];
        let numerator = |n_1: f64, n_3over2: f64, n_2: f64, n_5over2: f64, n_3: f64| {
            1. + powers.two_thirds * n_1 * pi.two_thirds
                + mf * n_3over2 * PI
                + powers.four_thirds * n_2 * pi.four_thirds
                + powers.five_thirds * n_5over2 * pi.five_thirds
                + mf.powf(2.) * n_3 * pi.two
        };
        let denominator = |d_1: f64, d_3over2: f64| {
            1. + d_1 * powers.two_thirds * pi.two_thirds + d_3over2 * mf * PI
        };

        // the derivative of the enhancement factor
        let s2s3b = -s2 * s3;
        let s2mfb = s2 * mf * PI * 2.0;
        let exps2s3b = s2s3b.cosh() + s2s3b.sinh();
        let exps2mfb = s2mfb.cosh() + s2mfb.sinh();
        let dynk2bar_deriv = (s1 - 1.0)
            * (2.0 * PI)
            * s2
            * exps2mfb
            * exps2s3b
            * (1. / ((1. + exps2mfb * exps2s3b) * (1. + exps2mfb * exps2s3b)))
            - 2.0 * PI * (s1 - 1.0) * s2 * exps2s3 / ((1.0 + exps2s3) * (1.0 + exps2s3));

        let mut terms = V3Terms::default();
        for (kappa, pn, n_5over2, n_3, d_1, n_1, n_3over2, n_2, d_3over2) in stars {
            let [c_newt, c_1, c_3over2, c_2, c_5over2] = [pn[0], pn[1], pn[2], pn[3], pn[4]];
            let ratio = numerator(n_1, n_3over2, n_2, n_5over2, n_3) / denominator(d_1, d_3over2);
            // Eqs. (30), (32) in terms of Mf
            terms.nr +=
                -((c_newt * (kappa * dynk2bar) * powers.five_thirds * pi.five_thirds) * ratio);
            let pn_poly = 1.
                + powers.two_thirds * c_1 * pi.two_thirds
                + mf * c_3over2 * PI
                + powers.four_thirds * c_2 * pi.four_thirds
                + powers.five_thirds * c_5over2 * pi.five_thirds;
            terms.pn += -(c_newt * kappa * powers.five_thirds * pi.five_thirds * pn_poly);

            let dnr1 =
                -((c_newt * kappa * dynk2bar_deriv * powers.five_thirds * pi.five_thirds) * ratio);
            let dnr2 = (c_newt
                * kappa
                * dynk2bar
                * powers.two_thirds
                * pi.five_thirds
                * (-5.
                    - powers.two_thirds * (3. * d_1 + 7. * n_1) * pi.two_thirds
                    - 2. * mf * (d_3over2 + 4. * n_3over2) * PI
                    - powers.four_thirds * (5. * d_1 * n_1 + 9. * n_2) * pi.four_thirds
                    - 2. * powers.five_thirds
                        * (2. * d_3over2 * n_1 + 3. * d_1 * n_3over2 + 5. * n_5over2)
                        * pi.five_thirds
                    - mf.powf(2.)
                        * (7. * d_1 * n_2 + 11. * n_3 + 5. * d_3over2 * n_3over2)
                        * pi.two
                    - 2. * powers.seven_thirds
                        * (3. * d_3over2 * n_2 + 4. * d_1 * n_5over2)
                        * pi.seven_thirds
                    - powers.eight_thirds
                        * (9. * d_1 * n_3 + 7. * d_3over2 * n_5over2)
                        * pi.eight_thirds
                    - 2. * mf.powf(3.) * (4. * d_3over2 * n_3) * pi.three))
                / (3. * denominator(d_1, d_3over2).powf(2.));
            terms.dnr += dnr1 + dnr2;

            // LAL writes (4/3) here, an integer division that evaluates to
            // 1: kept to reproduce LAL (it multiplies a term that vanishes
            // where LAL evaluates the derivative, see
            // docs/lal-nrtidalv3-integer-division.md).
            let four_thirds_in_lal = 1.0;
            let dpn_poly = (2. / 3.) * mf.powf(-1.0 / 3.0) * c_1 * pi.two_thirds
                + c_3over2 * PI
                + four_thirds_in_lal * mf.powf(1.0 / 3.0) * c_2 * pi.four_thirds
                + (5. / 3.) * powers.two_thirds * c_5over2 * pi.five_thirds;
            let factor = -c_newt * kappa;
            terms.dpn += factor
                * ((5. / 3.) * powers.two_thirds * pi.five_thirds * pn_poly
                    + powers.five_thirds * pi.five_thirds * dpn_poly);
        }
        terms
    }
}

/// See [`Tidal::v3_terms`].
#[derive(Default)]
struct V3Terms {
    nr: f64,
    pn: f64,
    dnr: f64,
    dpn: f64,
}
