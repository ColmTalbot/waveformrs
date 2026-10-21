// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Orbit-averaged SpinTaylorT4 equations of motion.
//!
//! Ports of `XLALSimSpinTaylorEnergySpinDerivativeSetup`,
//! `XLALSimInspiralSpinTaylorT4Setup`, `XLALSimInspiralSpinDerivativesAvg`,
//! `XLALSimInspiralSpinTaylorT4DerivativesAvg`,
//! `XLALSimInspiralSetEnergyPNTermsAvg` and
//! `XLALSimInspiralSpinTaylorStoppingTest` from `LALSimInspiralSpinTaylor.c`
//! (lalsimulation 6.2.0), restricted to the orbit-averaged dynamics with
//! `phenomtp = 0` (the only configuration reachable from
//! `XLALSimInspiralSpinTaylorPNEvolveOrbit`).
//!
//! The state vector is LAL's: `[phi, M omega, LNhat (3), S1 (3), S2 (3), E1 (3)]`
//! where the spins are normalised by the total mass, `S_i = (m_i / M)^2 chi_i`,
//! and time is measured in units of the total mass.

use super::coefficients::*;
use super::fp::{cross, dot, fma, fnms};
use super::integrator::OdeSystem;
use super::SpinTaylorError;
use crate::constants::{MSUN_SI, MTSUN_SI, REAL4_EPS};
use std::f64::consts::PI;

/// Number of dynamical variables.
pub const NUM_VARIABLES: usize = 14;

/// LAL termination and failure codes (`LALSIMINSPIRAL_ST_*`).
pub(crate) const ST_TEST_ENERGY: i32 = 1025;
pub(crate) const ST_TEST_OMEGADOUBLEDOT: i32 = 1026;
pub(crate) const ST_TEST_OMEGANAN: i32 = 1028;
pub(crate) const ST_TEST_FREQBOUND: i32 = 1029;
pub(crate) const ST_DERIVATIVE_OMEGANONPOS: i32 = 1030;
pub(crate) const ST_TEST_LARGEV: i32 = 1031;

/// LAL's switch statements over PN orders fall through from the highest
/// order, with `-1` meaning "all orders".
#[inline(always)]
fn includes(order: i32, twice_pn: i32) -> bool {
    order < 0 || order >= twice_pn
}

/// The non-dynamical coefficients of the SpinTaylorT4 equations, the port of
/// the relevant parts of `XLALSimInspiralSpinTaylorTxCoeffs`.
///
/// Sums of coefficients that LAL adds at every derivative evaluation (e.g.
/// `wdot6S1S1Avg + wdot6QMS1S1Avg`) are stored pre-added; this is bit-for-bit
/// identical since the same two numbers are added.
#[derive(Clone, Debug)]
pub struct SpinTaylorT4Coefficients {
    /// Total mass in solar masses.
    pub total_mass: f64,
    pub eta: f64,
    pub m1_by_m: f64,
    pub m2_by_m: f64,
    pub f_start: f64,
    pub f_end: f64,
    pub phase_order: i32,
    pub spin_order: i32,
    pub tidal_order: i32,
    pub lscorr: bool,

    wdotnewt: f64,
    wdotcoeff: [f64; 8],
    wdotlogcoeff: f64,
    wdot3s1o: f64,
    wdot3s2o: f64,
    wdot4s1s2_avg: f64,
    wdot4s1os2o_avg: f64,
    /// `wdot4S1S1Avg + wdot4QMS1S1Avg`
    wdot4s1s1_sum: f64,
    wdot4s2s2_sum: f64,
    wdot4s1os1o_sum: f64,
    wdot4s2os2o_sum: f64,
    wdot5s1o: f64,
    wdot5s2o: f64,
    wdot6s1o: f64,
    wdot6s2o: f64,
    wdot6s1s2_avg: f64,
    wdot6s1os2o_avg: f64,
    wdot6s1s1_sum: f64,
    wdot6s2s2_sum: f64,
    wdot6s1os1o_sum: f64,
    wdot6s2os2o_sum: f64,
    wdottidal10: f64,
    wdottidal12: f64,

    ecoeff: [f64; 8],
    e3s1o: f64,
    e3s2o: f64,
    e4s1s2_avg: f64,
    e4s1os2o_avg: f64,
    e4qms1s1_avg: f64,
    e4qms1os1o_avg: f64,
    e4qms2s2_avg: f64,
    e4qms2os2o_avg: f64,
    e5s1o: f64,
    e5s2o: f64,
    e6s1s2_avg: f64,
    e6s1os2o_avg: f64,
    e6s1s1_sum: f64,
    e6s2s2_sum: f64,
    e6s1os1o_sum: f64,
    e6s2os2o_sum: f64,
    etidal10: f64,
    etidal12: f64,

    s1dot3: f64,
    s2dot3: f64,
    s1dot4s2_avg: f64,
    s1dot4s2o_avg: f64,
    s1dot4qms1o_avg: f64,
    s2dot4qms2o_avg: f64,
    s1dot5: f64,
    s2dot5: f64,
    s1dot6s2_avg: f64,
    s1dot6s2o_avg: f64,
    s1dot6s1o_avg: f64,
    s1dot6qms1o_avg: f64,
    s2dot6s1_avg: f64,
    s2dot6s1o_avg: f64,
    s2dot6s2o_avg: f64,
    s2dot6qms2o_avg: f64,
    omegashift_s1: f64,
    omegashift_s2: f64,

    // Mass-dependent constants LAL recomputes in every call of
    // `XLALSimInspiralSpinDerivativesAvg` / the stopping test.
    c_s1: f64,
    c_s1l: f64,
    c_s2: f64,
    c_s2l: f64,
    l1pn: f64,
    omega_start: f64,
    omega_end: f64,
}

impl SpinTaylorT4Coefficients {
    /// Port of `XLALSimInspiralSpinTaylorT4Setup` with `phenomtp = 0`.
    ///
    /// Masses are in kg, frequencies in Hz.  The orders use LAL's
    /// conventions (twice the PN order, `-1` for all available terms).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        m1_si: f64,
        m2_si: f64,
        f_start: f64,
        f_end: f64,
        lambda1: f64,
        lambda2: f64,
        quadparam1: f64,
        quadparam2: f64,
        spin_order: i32,
        tidal_order: i32,
        phase_order: i32,
        lscorr: bool,
    ) -> Result<Self, SpinTaylorError> {
        if !(-1..=8).contains(&phase_order) {
            return Err(SpinTaylorError::InvalidPhaseOrder(phase_order));
        }
        // spinO = 7 is rejected by XLALSimInspiralSpinTaylorPNEvolveOrbit for
        // the orbit-averaged dynamics.
        if !(-1..=6).contains(&spin_order) {
            return Err(SpinTaylorError::InvalidSpinOrder(spin_order));
        }
        if !matches!(tidal_order, -1 | 0 | 10 | 12) {
            return Err(SpinTaylorError::InvalidTidalOrder(tidal_order));
        }

        let m1 = m1_si / MSUN_SI;
        let m2 = m2_si / MSUN_SI;
        let total_mass = m1 + m2;
        let m1m = m1 / total_mass;
        let m2m = m2 / total_mass;
        let eta = m1m * m2m;

        let mut ecoeff = [0.0; 8];
        if includes(phase_order, 6) {
            ecoeff[6] = pn_energy_6pn(eta);
        }
        if includes(phase_order, 4) {
            ecoeff[4] = pn_energy_4pn(eta);
        }
        if includes(phase_order, 2) {
            ecoeff[2] = pn_energy_2pn(eta);
        }
        ecoeff[0] = 1.;

        let mut wdotcoeff = [0.0; 8];
        let mut wdotlogcoeff = 0.0;
        if includes(phase_order, 7) {
            wdotcoeff[7] = taylor_t4_wdot_7pn(eta);
        }
        if includes(phase_order, 6) {
            wdotcoeff[6] = taylor_t4_wdot_6pn(eta);
            wdotlogcoeff = taylor_t4_wdot_6pn_log(eta);
        }
        if includes(phase_order, 5) {
            wdotcoeff[5] = taylor_t4_wdot_5pn(eta);
        }
        if includes(phase_order, 4) {
            wdotcoeff[4] = taylor_t4_wdot_4pn(eta);
        }
        if includes(phase_order, 3) {
            wdotcoeff[3] = taylor_t4_wdot_3pn(eta);
        }
        if includes(phase_order, 2) {
            wdotcoeff[2] = taylor_t4_wdot_2pn(eta);
        }
        wdotcoeff[0] = 1.;

        let mut new = Self {
            total_mass,
            eta,
            m1_by_m: m1m,
            m2_by_m: m2m,
            f_start,
            f_end,
            phase_order,
            spin_order,
            tidal_order,
            lscorr,
            wdotnewt: taylor_t4_wdot_0pn(eta),
            wdotcoeff,
            wdotlogcoeff,
            wdot3s1o: 0.0,
            wdot3s2o: 0.0,
            wdot4s1s2_avg: 0.0,
            wdot4s1os2o_avg: 0.0,
            wdot4s1s1_sum: 0.0,
            wdot4s2s2_sum: 0.0,
            wdot4s1os1o_sum: 0.0,
            wdot4s2os2o_sum: 0.0,
            wdot5s1o: 0.0,
            wdot5s2o: 0.0,
            wdot6s1o: 0.0,
            wdot6s2o: 0.0,
            wdot6s1s2_avg: 0.0,
            wdot6s1os2o_avg: 0.0,
            wdot6s1s1_sum: 0.0,
            wdot6s2s2_sum: 0.0,
            wdot6s1os1o_sum: 0.0,
            wdot6s2os2o_sum: 0.0,
            wdottidal10: 0.0,
            wdottidal12: 0.0,
            ecoeff,
            e3s1o: 0.0,
            e3s2o: 0.0,
            e4s1s2_avg: 0.0,
            e4s1os2o_avg: 0.0,
            e4qms1s1_avg: 0.0,
            e4qms1os1o_avg: 0.0,
            e4qms2s2_avg: 0.0,
            e4qms2os2o_avg: 0.0,
            e5s1o: 0.0,
            e5s2o: 0.0,
            e6s1s2_avg: 0.0,
            e6s1os2o_avg: 0.0,
            e6s1s1_sum: 0.0,
            e6s2s2_sum: 0.0,
            e6s1os1o_sum: 0.0,
            e6s2os2o_sum: 0.0,
            etidal10: 0.0,
            etidal12: 0.0,
            s1dot3: 0.0,
            s2dot3: 0.0,
            s1dot4s2_avg: 0.0,
            s1dot4s2o_avg: 0.0,
            s1dot4qms1o_avg: 0.0,
            s2dot4qms2o_avg: 0.0,
            s1dot5: 0.0,
            s2dot5: 0.0,
            s1dot6s2_avg: 0.0,
            s1dot6s2o_avg: 0.0,
            s1dot6s1o_avg: 0.0,
            s1dot6qms1o_avg: 0.0,
            s2dot6s1_avg: 0.0,
            s2dot6s1o_avg: 0.0,
            s2dot6s2o_avg: 0.0,
            s2dot6qms2o_avg: 0.0,
            omegashift_s1: 0.0,
            omegashift_s2: 0.0,
            c_s1: l_3pn_si_coeff_avg(m1m),
            c_s1l: l_3pn_sil_coeff_avg(m1m),
            c_s2: l_3pn_si_coeff_avg(m2m),
            c_s2l: l_3pn_sil_coeff_avg(m2m),
            l1pn: l_2pn(eta),
            // omega = PI G M f_GW / c^3, as in the stopping test
            omega_start: PI * total_mass * MTSUN_SI * f_start,
            omega_end: PI * total_mass * MTSUN_SI * f_end,
        };

        // XLALSimSpinTaylorEnergySpinDerivativeSetup: spin-dependent energy
        // and spin-precession coefficients.
        if includes(spin_order, 6) {
            new.e6s1s2_avg = pn_energy_6pn_s1s2_avg(eta);
            new.e6s1os2o_avg = pn_energy_6pn_s1os2o_avg(eta);
            new.e6s1s1_sum =
                pn_energy_6pn_s1s1_avg(m1m) + quadparam1 * pn_energy_6pn_qms1s1_avg(m1m);
            new.e6s2s2_sum =
                pn_energy_6pn_s1s1_avg(m2m) + quadparam2 * pn_energy_6pn_qms1s1_avg(m2m);
            new.e6s1os1o_sum =
                pn_energy_6pn_s1os1o_avg(m1m) + quadparam1 * pn_energy_6pn_qms1os1o_avg(m1m);
            new.e6s2os2o_sum =
                pn_energy_6pn_s1os1o_avg(m2m) + quadparam2 * pn_energy_6pn_qms1os1o_avg(m2m);
            new.s1dot6s2_avg = spin_dot_6pn_s2_avg(m1m);
            new.s2dot6s1_avg = spin_dot_6pn_s2_avg(m2m);
            new.s1dot6s2o_avg = spin_dot_6pn_s2o_avg(m1m);
            new.s1dot6s1o_avg = spin_dot_6pn_s1o_avg(m1m);
            new.s2dot6s1o_avg = spin_dot_6pn_s2o_avg(m2m);
            new.s2dot6s2o_avg = spin_dot_6pn_s1o_avg(m2m);
            new.s1dot6qms1o_avg = quadparam1 * spin_dot_6pn_qmso_avg(m1m);
            new.s2dot6qms2o_avg = quadparam2 * spin_dot_6pn_qmso_avg(m2m);
            new.omegashift_s1 = ldot_3pn_so(m1m);
            new.omegashift_s2 = ldot_3pn_so(m2m);
        }
        if includes(spin_order, 5) {
            new.e5s1o = pn_energy_5pn_so(m1m);
            new.e5s2o = pn_energy_5pn_so(m2m);
            new.s1dot5 = spin_dot_5pn(m1m);
            new.s2dot5 = spin_dot_5pn(m2m);
        }
        if includes(spin_order, 4) {
            new.e4s1s2_avg = pn_energy_4pn_s1s2_avg(eta);
            new.e4s1os2o_avg = pn_energy_4pn_s1os2o_avg(eta);
            new.e4qms1s1_avg = quadparam1 * pn_energy_4pn_qms1s1_avg(m1m);
            new.e4qms1os1o_avg = quadparam1 * pn_energy_4pn_qms1os1o_avg(m1m);
            new.e4qms2s2_avg = quadparam2 * pn_energy_4pn_qms1s1_avg(m2m);
            new.e4qms2os2o_avg = quadparam2 * pn_energy_4pn_qms1os1o_avg(m2m);
            new.s1dot4s2_avg = SPIN_DOT_4PN_S2_AVG;
            new.s1dot4s2o_avg = SPIN_DOT_4PN_S2O_AVG;
            new.s1dot4qms1o_avg = quadparam1 * spin_dot_4pn_qmso_avg(m1m);
            new.s2dot4qms2o_avg = quadparam2 * spin_dot_4pn_qmso_avg(m2m);
        }
        if includes(spin_order, 3) {
            new.e3s1o = pn_energy_3pn_so(m1m);
            new.e3s2o = pn_energy_3pn_so(m2m);
            new.s1dot3 = spin_dot_3pn(m1m);
            new.s2dot3 = spin_dot_3pn(m2m);
        }
        if includes(tidal_order, 12) {
            new.etidal12 = fma(
                lambda1,
                pn_energy_12pn_tidal(m1m),
                lambda2 * pn_energy_12pn_tidal(m2m),
            );
        }
        if includes(tidal_order, 10) {
            new.etidal10 = fma(
                lambda1,
                pn_energy_10pn_tidal(m1m),
                lambda2 * pn_energy_10pn_tidal(m2m),
            );
        }

        // XLALSimInspiralSpinTaylorT4Setup: spin corrections to domega/dt.
        if includes(spin_order, 6) {
            new.wdot6s1o = taylor_t4_wdot_6pn_so(m1m);
            new.wdot6s2o = taylor_t4_wdot_6pn_so(m2m);
            new.wdot6s1s2_avg = taylor_t4_wdot_6pn_s1s2_avg(eta);
            new.wdot6s1os2o_avg = taylor_t4_wdot_6pn_s1os2o_avg(eta);
            new.wdot6s1s1_sum =
                taylor_t4_wdot_6pn_s1s1_avg(m1m) + quadparam1 * taylor_t4_wdot_6pn_qms1s1_avg(m1m);
            new.wdot6s2s2_sum =
                taylor_t4_wdot_6pn_s1s1_avg(m2m) + quadparam2 * taylor_t4_wdot_6pn_qms1s1_avg(m2m);
            new.wdot6s1os1o_sum = taylor_t4_wdot_6pn_s1os1o_avg(m1m)
                + quadparam1 * taylor_t4_wdot_6pn_qms1os1o_avg(m1m);
            new.wdot6s2os2o_sum = taylor_t4_wdot_6pn_s1os1o_avg(m2m)
                + quadparam2 * taylor_t4_wdot_6pn_qms1os1o_avg(m2m);
        }
        if includes(spin_order, 5) {
            new.wdot5s1o = taylor_t4_wdot_5pn_so(m1m);
            new.wdot5s2o = taylor_t4_wdot_5pn_so(m2m);
        }
        if includes(spin_order, 4) {
            new.wdot4s1s2_avg = taylor_t4_wdot_4pn_s1s2_avg(eta);
            new.wdot4s1os2o_avg = taylor_t4_wdot_4pn_s1os2o_avg(eta);
            new.wdot4s1s1_sum =
                taylor_t4_wdot_4pn_s1s1_avg(m1m) + quadparam1 * taylor_t4_wdot_4pn_qms1s1_avg(m1m);
            new.wdot4s1os1o_sum = taylor_t4_wdot_4pn_s1os1o_avg(m1m)
                + quadparam1 * taylor_t4_wdot_4pn_qms1os1o_avg(m1m);
            new.wdot4s2s2_sum =
                taylor_t4_wdot_4pn_s1s1_avg(m2m) + quadparam2 * taylor_t4_wdot_4pn_qms1s1_avg(m2m);
            new.wdot4s2os2o_sum = taylor_t4_wdot_4pn_s1os1o_avg(m2m)
                + quadparam2 * taylor_t4_wdot_4pn_qms1os1o_avg(m2m);
        }
        if includes(spin_order, 3) {
            new.wdot3s1o = taylor_t4_wdot_3pn_so(m1m);
            new.wdot3s2o = taylor_t4_wdot_3pn_so(m2m);
        }
        if includes(tidal_order, 12) {
            new.wdottidal12 = fma(
                lambda1,
                taylor_t4_wdot_12pn_tidal(m1m),
                lambda2 * taylor_t4_wdot_12pn_tidal(m2m),
            );
        }
        if includes(tidal_order, 10) {
            new.wdottidal10 = fma(
                lambda1,
                taylor_t4_wdot_10pn_tidal(m1m),
                lambda2 * taylor_t4_wdot_10pn_tidal(m2m),
            );
        }
        Ok(new)
    }

    /// Port of `XLALSimInspiralSpinDerivativesAvg`: derivatives of
    /// `LNhat`, `E1`, `S1` and `S2`, written into `dvalues[2..14]`.
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    fn spin_derivatives(
        &self,
        v: f64,
        lnh: &[f64; 3],
        e1: &[f64; 3],
        s1: &[f64; 3],
        s2: &[f64; 3],
        lnh_dot_s1: f64,
        lnh_dot_s2: f64,
        dvalues: &mut [f64; NUM_VARIABLES],
    ) {
        let spin_order = self.spin_order;
        let eta = self.eta;
        let ln0mag = eta / v;
        let mut lnmag = ln0mag;

        let mut ds1 = [0.0; 3];
        let mut ds2 = [0.0; 3];
        let mut dlnhat = [0.0; 3];

        if includes(spin_order, 3) {
            let v2 = v * v;
            let omega = v2 * v;
            let v5 = omega * v2;

            let lnh_x_s1 = cross(lnh, s1);
            let lnh_x_s2 = cross(lnh, s2);
            let mut ds1_l = [0.0; 3];
            let mut ds2_l = [0.0; 3];
            let mut dl_l = [0.0; 3];
            for i in 0..3 {
                ds1_l[i] = self.s1dot3 * v5 * lnh_x_s1[i];
                ds2_l[i] = self.s2dot3 * v5 * lnh_x_s2[i];
                ds1[i] = ds1_l[i];
                ds2[i] = ds2_l[i];
                dl_l[i] = -(ds1_l[i] + ds2_l[i]);
                dlnhat[i] = dl_l[i];
            }

            if includes(spin_order, 4) {
                let omega2 = omega * omega;
                let s1_x_s2 = cross(s1, s2);
                let mut ds1_nl = [0.0; 3];
                let mut ds2_nl = [0.0; 3];
                for i in 0..3 {
                    ds1_nl[i] = omega2
                        * fma(
                            -self.s1dot4s2_avg,
                            s1_x_s2[i],
                            self.s1dot4s2o_avg * lnh_dot_s2 * lnh_x_s1[i],
                        );
                    ds2_nl[i] = omega2
                        * fma(
                            self.s1dot4s2_avg,
                            s1_x_s2[i],
                            self.s1dot4s2o_avg * lnh_dot_s1 * lnh_x_s2[i],
                        );
                }
                for i in 0..3 {
                    ds1_nl[i] = fma(
                        omega2 * self.s1dot4qms1o_avg * lnh_dot_s1,
                        lnh_x_s1[i],
                        ds1_nl[i],
                    );
                }
                for i in 0..3 {
                    ds2_nl[i] = fma(
                        omega2 * self.s2dot4qms2o_avg * lnh_dot_s2,
                        lnh_x_s2[i],
                        ds2_nl[i],
                    );
                }
                for i in 0..3 {
                    ds1[i] += ds1_nl[i];
                    ds2[i] += ds2_nl[i];
                    let dl_nl = -(ds1_nl[i] + ds2_nl[i]);
                    dlnhat[i] += dl_nl;
                }

                if includes(spin_order, 5) {
                    let v7 = omega2 * v;
                    lnmag = fma(ln0mag * v2, self.l1pn, lnmag);

                    for i in 0..3 {
                        let ds1_nnl = self.s1dot5 * v7 * lnh_x_s1[i];
                        let ds2_nnl = self.s2dot5 * v7 * lnh_x_s2[i];
                        ds1[i] += ds1_nnl;
                        ds2[i] += ds2_nnl;
                        dlnhat[i] -= ds1_nnl + ds2_nnl;
                    }

                    if self.lscorr {
                        for i in 0..3 {
                            dlnhat[i] = fnms(
                                eta * v2,
                                fma(self.c_s1, ds1_l[i], self.c_s2 * ds2_l[i]),
                                dlnhat[i],
                            );
                        }
                    }

                    if includes(spin_order, 6) {
                        let v8 = omega2 * v2;
                        let q1 = fma(
                            self.s1dot6s1o_avg,
                            lnh_dot_s1,
                            self.s1dot6s2o_avg * lnh_dot_s2,
                        );
                        let q2 = fma(
                            self.s2dot6s1o_avg,
                            lnh_dot_s1,
                            self.s2dot6s2o_avg * lnh_dot_s2,
                        );
                        let mut ds1_n3l = [0.0; 3];
                        let mut ds2_n3l = [0.0; 3];
                        for i in 0..3 {
                            ds1_n3l[i] = v8 * fma(-self.s1dot6s2_avg, s1_x_s2[i], q1 * lnh_x_s1[i]);
                            ds2_n3l[i] = v8 * fma(self.s2dot6s1_avg, s1_x_s2[i], q2 * lnh_x_s2[i]);
                        }
                        for i in 0..3 {
                            ds1_n3l[i] = fma(
                                v8,
                                self.s1dot6qms1o_avg * lnh_dot_s1 * lnh_x_s1[i],
                                ds1_n3l[i],
                            );
                        }
                        for i in 0..3 {
                            ds2_n3l[i] = fma(
                                v8,
                                self.s2dot6qms2o_avg * lnh_dot_s2 * lnh_x_s2[i],
                                ds2_n3l[i],
                            );
                        }
                        for i in 0..3 {
                            ds1[i] += ds1_n3l[i];
                            ds2[i] += ds2_n3l[i];
                            dlnhat[i] -= ds1_n3l[i] + ds2_n3l[i];
                        }

                        if self.lscorr {
                            let cl = fma(self.c_s1l, lnh_dot_s1, self.c_s2l * lnh_dot_s2);
                            for i in 0..3 {
                                let inner = fma(self.c_s1, ds1_nl[i], self.c_s2 * ds2_nl[i])
                                    + cl * dl_l[i] / ln0mag;
                                dlnhat[i] = fnms(eta * v2, inner, dlnhat[i]);
                            }
                        }
                    }
                }
            }
        }

        for d in dlnhat.iter_mut() {
            *d /= lnmag;
        }
        let om = cross(lnh, &dlnhat);
        dvalues[2] = fma(-om[2], lnh[1], om[1] * lnh[2]);
        dvalues[3] = fma(-om[0], lnh[2], om[2] * lnh[0]);
        dvalues[4] = fma(-om[1], lnh[0], om[0] * lnh[1]);
        dvalues[5] = ds1[0];
        dvalues[6] = ds1[1];
        dvalues[7] = ds1[2];
        dvalues[8] = ds2[0];
        dvalues[9] = ds2[1];
        dvalues[10] = ds2[2];
        dvalues[11] = fma(-om[2], e1[1], om[1] * e1[2]);
        dvalues[12] = fma(-om[0], e1[2], om[2] * e1[0]);
        dvalues[13] = fma(-om[1], e1[0], om[0] * e1[1]);
    }

    /// Port of `XLALSimInspiralSpinTaylorT4DerivativesAvg`.
    ///
    /// Returns `0` on success or LAL's error code if `omega <= 0`.
    #[inline(always)]
    pub fn derivatives(
        &self,
        values: &[f64; NUM_VARIABLES],
        dvalues: &mut [f64; NUM_VARIABLES],
    ) -> i32 {
        let omega = values[1];
        let lnh = [values[2], values[3], values[4]];
        let s1 = [values[5], values[6], values[7]];
        let s2 = [values[8], values[9], values[10]];
        let e1 = [values[11], values[12], values[13]];

        if omega <= 0.0 {
            return ST_DERIVATIVE_OMEGANONPOS;
        }

        let v = omega.cbrt();
        let v2 = v * v;
        let v11 = omega * omega * omega * v2;

        let lnh_dot_s1 = dot(&lnh, &s1);
        let lnh_dot_s2 = dot(&lnh, &s2);
        let s1_dot_s2 = dot(&s1, &s2);
        let s1sq = dot(&s1, &s1);
        let s2sq = dot(&s2, &s2);

        let spin_order = self.spin_order;
        let mut wspin3 = 0.0;
        let mut wspin4 = 0.0;
        let mut wspin5 = 0.0;
        let mut wspin6 = 0.0;
        if includes(spin_order, 6) {
            let mut w = fma(self.wdot6s1o, lnh_dot_s1, self.wdot6s2o * lnh_dot_s2);
            w = fma(self.wdot6s1os2o_avg * lnh_dot_s1, lnh_dot_s2, w);
            w = fma(self.wdot6s1s2_avg, s1_dot_s2, w);
            w = fma(self.wdot6s1s1_sum, s1sq, w);
            w = fma(self.wdot6s2s2_sum, s2sq, w);
            w = fma(self.wdot6s1os1o_sum * lnh_dot_s1, lnh_dot_s1, w);
            w = fma(self.wdot6s2os2o_sum * lnh_dot_s2, lnh_dot_s2, w);
            wspin6 = w;
        }
        if includes(spin_order, 5) {
            wspin5 = fma(self.wdot5s1o, lnh_dot_s1, self.wdot5s2o * lnh_dot_s2);
        }
        if includes(spin_order, 4) {
            let w = fma(
                self.wdot4s1s2_avg,
                s1_dot_s2,
                self.wdot4s1os2o_avg * lnh_dot_s1 * lnh_dot_s2,
            );
            let mut x = fma(self.wdot4s1s1_sum, s1sq, self.wdot4s2s2_sum * s2sq);
            x = fma(self.wdot4s1os1o_sum * lnh_dot_s1, lnh_dot_s1, x);
            x = fma(self.wdot4s2os2o_sum * lnh_dot_s2, lnh_dot_s2, x);
            wspin4 = w + x;
        }
        if includes(spin_order, 3) {
            wspin3 = fma(self.wdot3s1o, lnh_dot_s1, self.wdot3s2o * lnh_dot_s2);
        }

        let c = &self.wdotcoeff;
        let mut p = fma(v2, self.wdottidal12, self.wdottidal10);
        p = fma(omega, p, c[7]);
        p = fma(v, p, fma(self.wdotlogcoeff, v.ln(), c[6] + wspin6));
        p = fma(v, p, c[5] + wspin5);
        p = fma(v, p, c[4] + wspin4);
        p = fma(v, p, c[3] + wspin3);
        p = fma(v, p, c[2]);
        p = fma(v, p, c[1]);
        p = fma(v, p, c[0]);
        let domega = self.wdotnewt * v11 * p;

        self.spin_derivatives(v, &lnh, &e1, &s1, &s2, lnh_dot_s1, lnh_dot_s2, dvalues);

        let shift = omegashift(
            s1sq,
            s2sq,
            s1_dot_s2,
            lnh_dot_s1,
            lnh_dot_s2,
            self.omegashift_s1,
            self.omegashift_s2,
        );
        dvalues[0] = omega * fma(omega * omega, shift, 1.);
        dvalues[1] = domega;
        0
    }

    /// Port of `XLALSimInspiralSetEnergyPNTermsAvg`, returning
    /// `[Espin3, Espin4, Espin5, Espin6]` (`Espin7` vanishes for
    /// `phenomtp = 0`).
    #[inline]
    fn energy_spin_terms(
        &self,
        lnh_dot_s1: f64,
        lnh_dot_s2: f64,
        s1sq: f64,
        s2sq: f64,
        s1_dot_s2: f64,
    ) -> [f64; 4] {
        let spin_order = self.spin_order;
        let mut espin = [0.0; 4];
        if includes(spin_order, 6) {
            let mut e = fma(
                self.e6s1s2_avg,
                s1_dot_s2,
                self.e6s1os2o_avg * lnh_dot_s1 * lnh_dot_s2,
            );
            e = fma(self.e6s1s1_sum, s1sq, e);
            e = fma(self.e6s2s2_sum, s2sq, e);
            e = fma(self.e6s1os1o_sum * lnh_dot_s1, lnh_dot_s1, e);
            e = fma(self.e6s2os2o_sum * lnh_dot_s2, lnh_dot_s2, e);
            espin[3] += e;
        }
        if includes(spin_order, 5) {
            espin[2] += fma(self.e5s1o, lnh_dot_s1, self.e5s2o * lnh_dot_s2);
        }
        if includes(spin_order, 4) {
            espin[1] += fma(
                self.e4s1s2_avg,
                s1_dot_s2,
                self.e4s1os2o_avg * lnh_dot_s1 * lnh_dot_s2,
            );
            let mut e = fma(self.e4qms1s1_avg, s1sq, self.e4qms2s2_avg * s2sq);
            e = fma(self.e4qms1os1o_avg * lnh_dot_s1, lnh_dot_s1, e);
            e = fma(self.e4qms2os2o_avg * lnh_dot_s2, lnh_dot_s2, e);
            espin[1] += e;
        }
        if includes(spin_order, 3) {
            espin[0] += fma(self.e3s1o, lnh_dot_s1, self.e3s2o * lnh_dot_s2);
        }
        espin
    }

    /// Port of `XLALSimInspiralSpinTaylorStoppingTest`.
    ///
    /// `prev_domega` is the state LAL keeps in the coefficient struct.
    #[inline]
    #[allow(clippy::if_same_then_else)] // mirrors the C branches
    pub fn stopping_test(
        &self,
        values: &[f64; NUM_VARIABLES],
        dvalues: &[f64; NUM_VARIABLES],
        prev_domega: &mut f64,
    ) -> i32 {
        let omega = values[1];
        let v = omega.cbrt();
        let lnh = [values[2], values[3], values[4]];
        let s1 = [values[5], values[6], values[7]];
        let s2 = [values[8], values[9], values[10]];
        let lnh_dot_s1 = dot(&lnh, &s1);
        let lnh_dot_s2 = dot(&lnh, &s2);
        let s1sq = dot(&s1, &s1);
        let s2sq = dot(&s2, &s2);
        let s1_dot_s2 = dot(&s1, &s2);

        let omega_start = self.omega_start;
        let omega_end = self.omega_end;

        let [espin3, espin4, espin5, espin6] =
            self.energy_spin_terms(lnh_dot_s1, lnh_dot_s2, s1sq, s2sq, s1_dot_s2);
        let espin7 = 0.0;
        let e = &self.ecoeff;
        let mut test = fma(12., self.etidal10, v * v * (14. * self.etidal12));
        test = fma(9., e[7] + espin7, v * v * v * test);
        test = fma(8., e[6] + espin6, v * test);
        test = fma(7., e[5] + espin5, v * test);
        test = fma(6., e[4] + espin4, v * test);
        test = fma(5., e[3] + espin3, v * test);
        test = fma(4., e[2], v * test);
        test = fma(v * v, test, 2.);

        let mut ddomega = dvalues[1] - *prev_domega;
        if self.f_end < self.f_start && self.f_end != 0. && *prev_domega != 0. {
            ddomega *= -1.;
        }
        *prev_domega = dvalues[1];

        if omega_end.abs() > REAL4_EPS && omega_end > omega_start && omega > omega_end {
            ST_TEST_FREQBOUND
        } else if omega_end.abs() > REAL4_EPS && omega_end < omega_start && omega < omega_end {
            ST_TEST_FREQBOUND
        } else if test < 0.0 {
            ST_TEST_ENERGY
        } else if omega.is_nan() {
            ST_TEST_OMEGANAN
        } else if v >= 1. {
            ST_TEST_LARGEV
        } else if ddomega <= 0. {
            ST_TEST_OMEGADOUBLEDOT
        } else {
            0
        }
    }
}

/// Port of the static `omegashift` helper: the 3PN shift between the orbital
/// frequency and the phase derivative due to precession.
#[inline(always)]
fn omegashift(
    s1sq: f64,
    s2sq: f64,
    s1s2: f64,
    lnhs1: f64,
    lnhs2: f64,
    oms1: f64,
    oms2: f64,
) -> f64 {
    let a1 = fnms(lnhs1, lnhs1, s1sq);
    let a2 = fnms(lnhs2, lnhs2, s2sq);
    let a3 = fnms(lnhs1, lnhs2, s1s2);
    let s = fma(oms1 * oms1, a1, oms2 * oms2 * a2);
    -0.25 * fma(2. * oms1 * oms2, a3, s)
}

/// The SpinTaylorT4 system as seen by the integrator: the coefficients plus
/// the mutable state of the stopping test.
pub(crate) struct SpinTaylorT4System<'a> {
    pub coefficients: &'a SpinTaylorT4Coefficients,
    pub prev_domega: f64,
}

impl OdeSystem<NUM_VARIABLES> for SpinTaylorT4System<'_> {
    #[inline(always)]
    fn derivatives(
        &self,
        _t: f64,
        y: &[f64; NUM_VARIABLES],
        dydt: &mut [f64; NUM_VARIABLES],
    ) -> i32 {
        self.coefficients.derivatives(y, dydt)
    }

    #[inline]
    fn stop(&mut self, _t: f64, y: &[f64; NUM_VARIABLES], dydt: &[f64; NUM_VARIABLES]) -> i32 {
        self.coefficients
            .stopping_test(y, dydt, &mut self.prev_domega)
    }
}

/// Convert masses in solar masses to the geometrised time units used by
/// `XLALSimInspiralSpinTaylorPNEvolveOrbit`.
#[inline]
pub(crate) fn mass_in_seconds(m_si: f64) -> f64 {
    m_si / MSUN_SI * MTSUN_SI
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Case {
        /// m1_SI, m2_SI, fStart, fEnd, quadparam1, quadparam2, lambda1, lambda2
        reals: [f64; 8],
        /// phaseO, spinO, tideO
        orders: [i32; 3],
        lscorr: bool,
        values: [f64; NUM_VARIABLES],
        status: i32,
        /// `XLALSimInspiralSpinTaylorT4DerivativesAvg` (lalsimulation 6.2.0)
        expected: [u64; NUM_VARIABLES],
    }

    const CASES: [Case; 4] = [
        Case {
            reals: [
                7.158275534512983e31,
                5.766388624024347e31,
                20.0,
                400.0,
                1.0,
                1.0,
                0.0,
                0.0,
            ],
            orders: [7, 6, 12],
            lscorr: false,
            values: [
                1.5, 0.05, 0.1, 0.2, 0.97, 0.1, 0.05, 0.12, -0.03, 0.02, 0.04, 0.9, 0.1, 0.2,
            ],
            status: 0,
            expected: [
                0x3fa9996fcb75c1eb,
                0x3f1894a2bbb49b21,
                0x3f30a4a0722476a5,
                0xbf3408597b702052,
                0x3f035136dc50fb0d,
                0xbf217aee3d64ff01,
                0x3f3f25230fdb0155,
                0xbf16c658096f5935,
                0xbf13878cd48d7ca3,
                0xbf2d7fd23e014753,
                0x3f0db4513d2e53b5,
                0x3f0dfb671f842303,
                0xbf21b978a6f21482,
                0xbf28de17affb9d23,
            ],
        },
        Case {
            reals: [
                1.893723686379096e32,
                9.468618431895478e30,
                15.0,
                10.0,
                3.0,
                5.0,
                200.0,
                900.0,
            ],
            orders: [-1, -1, -1],
            lscorr: true,
            values: [
                -3.25, 0.12, -0.3, 0.25, 0.9, 0.6, -0.5, 0.2, 0.001, -0.002, 0.0015, 0.7, 0.6, 0.38,
            ],
            status: 0,
            expected: [
                0x3fbe63630bc4b76f,
                0x3f5fab5beaf294c5,
                0xbf919e23f1ccc5ea,
                0xbf94aa51d6bdad35,
                0xbf20ec66eb5726b9,
                0x3f62163dceb3a632,
                0x3f65ba6e1e2058b8,
                0x3edeb3be6bd66c02,
                0x3f1cff24b377ab25,
                0x3f111360477f694c,
                0x3eeb7de9e82627bf,
                0xbf89c59f0273f8f7,
                0x3f38758c9069e6e5,
                0x3f972246fd0326c5,
            ],
        },
        Case {
            reals: [2.0e31, 1.5e31, 20.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            orders: [5, 4, 0],
            lscorr: false,
            values: [
                0.0, 0.01, 0.0, 0.0, 1.0, 0.2, 0.0, 0.1, 0.0, -0.1, 0.05, 1.0, 0.0, 0.0,
            ],
            status: 0,
            expected: [
                0x3f847ae147ae147b,
                0x3e8c1a421224bdf0,
                0xbf03d185849cc3f8,
                0xbf0f09d39f795c83,
                0x0000000000000000,
                0xbea0c6f7a0b5ed8e,
                0x3f11c5a08ac07df0,
                0x3eb0c6f7a0b5ed8e,
                0x3f06ca4297c3a4ff,
                0xbea0c6f7a0b5ed8e,
                0xbeb0c6f7a0b5ed8e,
                0x8000000000000000,
                0x0000000000000000,
                0x3f03d185849cc3f8,
            ],
        },
        Case {
            // omega <= 0 is rejected
            reals: [2.0e31, 1.5e31, 20.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            orders: [7, 6, 12],
            lscorr: false,
            values: [
                0.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.0, 0.1, 0.0, -0.1, 0.05, 1.0, 0.0, 0.0,
            ],
            status: ST_DERIVATIVE_OMEGANONPOS,
            expected: [0; NUM_VARIABLES],
        },
    ];

    #[test]
    fn derivatives_match_lalsimulation() {
        for case in &CASES {
            let [m1, m2, f_start, f_end, q1, q2, l1, l2] = case.reals;
            let [phase_order, spin_order, tidal_order] = case.orders;
            let coefficients = SpinTaylorT4Coefficients::new(
                m1,
                m2,
                f_start,
                f_end,
                l1,
                l2,
                q1,
                q2,
                spin_order,
                tidal_order,
                phase_order,
                case.lscorr,
            )
            .unwrap();
            let mut dvalues = [0.0; NUM_VARIABLES];
            let status = coefficients.derivatives(&case.values, &mut dvalues);
            assert_eq!(status, case.status);
            for (value, expected) in dvalues.iter().zip(case.expected) {
                let expected = f64::from_bits(expected);
                assert!(
                    (value - expected).abs() <= 1e-12 * expected.abs(),
                    "{value:e} != {expected:e}"
                );
            }
        }
    }

    #[test]
    fn invalid_orders_are_rejected() {
        let new = |phase, spin, tidal| {
            SpinTaylorT4Coefficients::new(
                2e31, 1e31, 20.0, 100.0, 0.0, 0.0, 1.0, 1.0, spin, tidal, phase, false,
            )
        };
        assert!(new(7, 6, 12).is_ok());
        assert!(new(-1, -1, -1).is_ok());
        assert_eq!(
            new(9, 6, 12).unwrap_err(),
            SpinTaylorError::InvalidPhaseOrder(9)
        );
        assert_eq!(
            new(7, 7, 12).unwrap_err(),
            SpinTaylorError::InvalidSpinOrder(7)
        );
        assert_eq!(
            new(7, 6, 14).unwrap_err(),
            SpinTaylorError::InvalidTidalOrder(14)
        );
    }

    #[test]
    fn stopping_conditions() {
        let coefficients = SpinTaylorT4Coefficients::new(
            2e31, 1e31, 20.0, 100.0, 0.0, 0.0, 1.0, 1.0, 6, 12, 7, false,
        )
        .unwrap();
        let omega_end = coefficients.omega_end;
        let mut values = [0.0; NUM_VARIABLES];
        values[4] = 1.0;
        values[11] = 1.0;
        let mut dvalues = [0.0; NUM_VARIABLES];
        let mut prev = 0.0;

        // inspiralling: continue
        values[1] = 0.5 * omega_end;
        dvalues[1] = 1.0;
        assert_eq!(coefficients.stopping_test(&values, &dvalues, &mut prev), 0);
        assert_eq!(prev, 1.0);
        // domega/dt decreasing
        dvalues[1] = 0.5;
        assert_eq!(
            coefficients.stopping_test(&values, &dvalues, &mut prev),
            ST_TEST_OMEGADOUBLEDOT
        );
        // above the end frequency
        values[1] = 1.01 * omega_end;
        dvalues[1] = 2.0;
        assert_eq!(
            coefficients.stopping_test(&values, &dvalues, &mut prev),
            ST_TEST_FREQBOUND
        );
        // NaN frequency
        values[1] = f64::NAN;
        assert_eq!(
            coefficients.stopping_test(&values, &dvalues, &mut prev),
            ST_TEST_OMEGANAN
        );
    }
}
