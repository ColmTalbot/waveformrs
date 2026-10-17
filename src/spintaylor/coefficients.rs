//! Post-Newtonian coefficients used by the SpinTaylor family.
//!
//! Direct ports of the functions in `LALSimInspiralPNCoefficients.c`
//! (lalsimulation 6.2.0) that `XLALSimInspiralSpinTaylorT4Setup` uses.  The
//! names follow LAL with the `XLALSimInspiral` prefix dropped.  Operation
//! order (and fused multiply-adds, see [`super::fp`]) follows the C source.
#![allow(non_snake_case)]

use std::f64::consts::PI;

use super::fp::{fma, fnms};
use crate::constants::EGAMMA;

// Energy, Eq. (3.1) of arXiv:0907.0700 and spin terms.

pub(crate) fn pn_energy_2pn(eta: f64) -> f64 {
    -(0.75 + eta / 12.0)
}

pub(crate) fn pn_energy_4pn(eta: f64) -> f64 {
    -fma(1. / 24.0 * eta, eta, fnms(19.0 / 8.0, eta, 27.0 / 8.0))
}

pub(crate) fn pn_energy_6pn(eta: f64) -> f64 {
    const K: f64 = 344.45 / 5.76 - 20.5 / 9.6 * PI * PI;
    let a = fnms(K, eta, 67.5 / 6.4);
    let b = fma(15.5 / 9.6 * eta, eta, a);
    -fma(3.5 / 518.4 * eta * eta, eta, b)
}

pub(crate) fn pn_energy_3pn_so(m_by_m: f64) -> f64 {
    2. / 3. + 2. / m_by_m
}

pub(crate) fn pn_energy_4pn_s1s2_avg(eta: f64) -> f64 {
    1. / eta
}

pub(crate) fn pn_energy_4pn_s1os2o_avg(eta: f64) -> f64 {
    -3. / eta
}

pub(crate) fn pn_energy_4pn_qms1s1_avg(m_by_m: f64) -> f64 {
    0.5 / m_by_m / m_by_m
}

pub(crate) fn pn_energy_4pn_qms1os1o_avg(m_by_m: f64) -> f64 {
    -1.5 / m_by_m / m_by_m
}

pub(crate) fn pn_energy_5pn_so(m_by_m: f64) -> f64 {
    5. / 3. + 3. / m_by_m + 29. * m_by_m / 9. + m_by_m * m_by_m / 9.
}

pub(crate) fn pn_energy_6pn_s1s2_avg(eta: f64) -> f64 {
    2. / eta - 11. / 6.
}

pub(crate) fn pn_energy_6pn_s1os2o_avg(eta: f64) -> f64 {
    -11. / 3. / eta + 2.3 / 1.8
}

pub(crate) fn pn_energy_6pn_s1s1_avg(m_by_m: f64) -> f64 {
    -1. / (m_by_m * m_by_m) - 1. / 6. / m_by_m - 0.5
}

pub(crate) fn pn_energy_6pn_s1os1o_avg(m_by_m: f64) -> f64 {
    6. / (m_by_m * m_by_m) - 1.5 / m_by_m - 1.1 / 1.8
}

pub(crate) fn pn_energy_6pn_qms1s1_avg(m_by_m: f64) -> f64 {
    1.25 / m_by_m / m_by_m + 1.25 / m_by_m + 5. / 12.
}

pub(crate) fn pn_energy_6pn_qms1os1o_avg(m_by_m: f64) -> f64 {
    -3.75 / m_by_m / m_by_m - 3.75 / m_by_m - 1.25
}

pub(crate) fn pn_energy_10pn_tidal(m_by_m: f64) -> f64 {
    -9.0 * m_by_m * m_by_m * m_by_m * m_by_m * (1. - m_by_m)
}

pub(crate) fn pn_energy_12pn_tidal(m_by_m: f64) -> f64 {
    let a = fma(11. / 2., m_by_m, -33. / 2.);
    let b = fnms(11. / 2. * m_by_m, m_by_m, a);
    let c = fma(33. / 2. * m_by_m * m_by_m, m_by_m, b);
    c * m_by_m * m_by_m * m_by_m * m_by_m
}

// TaylorT4 domega/dt, Eq. (3.6) of arXiv:0907.0700 and spin terms.

pub(crate) fn taylor_t4_wdot_0pn(eta: f64) -> f64 {
    96. / 5. * eta
}

pub(crate) fn taylor_t4_wdot_2pn(eta: f64) -> f64 {
    -(1.0 / 336.0) * fma(924.0, eta, 743.0)
}

pub(crate) fn taylor_t4_wdot_3pn(_eta: f64) -> f64 {
    4.0 * PI
}

pub(crate) fn taylor_t4_wdot_3pn_so(m_by_m: f64) -> f64 {
    -19. / 6. - 25. / 4. / m_by_m
}

pub(crate) fn taylor_t4_wdot_4pn(eta: f64) -> f64 {
    fma(59472. * eta, eta, fma(122949., eta, 34103.)) / 18144.
}

pub(crate) fn taylor_t4_wdot_4pn_s1s2_avg(eta: f64) -> f64 {
    -247. / 48. / eta
}

pub(crate) fn taylor_t4_wdot_4pn_s1os2o_avg(eta: f64) -> f64 {
    721. / 48. / eta
}

pub(crate) fn taylor_t4_wdot_4pn_s1s1_avg(m_by_m: f64) -> f64 {
    7. / 96. / m_by_m / m_by_m
}

pub(crate) fn taylor_t4_wdot_4pn_s1os1o_avg(m_by_m: f64) -> f64 {
    -1. / 96. / m_by_m / m_by_m
}

pub(crate) fn taylor_t4_wdot_4pn_qms1s1_avg(m_by_m: f64) -> f64 {
    -2.5 / m_by_m / m_by_m
}

pub(crate) fn taylor_t4_wdot_4pn_qms1os1o_avg(m_by_m: f64) -> f64 {
    7.5 / m_by_m / m_by_m
}

pub(crate) fn taylor_t4_wdot_5pn(eta: f64) -> f64 {
    -(1.0 / 672.0) * PI * fma(15876.0, eta, 4159.0)
}

pub(crate) fn taylor_t4_wdot_5pn_so(m_by_m: f64) -> f64 {
    -809. / (84. * m_by_m) + 13.795 / 1.008 - 527. * m_by_m / 24. - 79. * m_by_m * m_by_m / 6.
}

pub(crate) fn taylor_t4_wdot_6pn(eta: f64) -> f64 {
    // The leading constants and `856/105 log(16)` are folded at compile time
    // in C, the remaining products are contracted with their sums.
    const K0: f64 = 16447.322263 / 139.7088 - 1712. / 105. * EGAMMA;
    let t1 = fnms(561.98689 / 2.17728, eta, K0);
    let t2 = fma(PI * PI, fma(451. / 48., eta, 16. / 3.), t1);
    let t3 = fma(541. / 896. * eta, eta, t2);
    let t4 = fnms(5605. / 2592. * eta * eta, eta, t3);
    t4 - 856. / 105. * LN_16
}

/// `log(16.)`, folded at compile time by the C compiler.
const LN_16: f64 = 2.772588722239781;

pub(crate) fn taylor_t4_wdot_6pn_log(_eta: f64) -> f64 {
    -(1712.0 / 105.0)
}

pub(crate) fn taylor_t4_wdot_6pn_so(m_by_m: f64) -> f64 {
    PI * (-37. / 3. - 151. / 6. / m_by_m)
}

pub(crate) fn taylor_t4_wdot_6pn_s1s2_avg(eta: f64) -> f64 {
    108.79 / (6.72 * eta) + 75.25 / 2.88
}

pub(crate) fn taylor_t4_wdot_6pn_s1os2o_avg(eta: f64) -> f64 {
    162.25 / (2.24 * eta) - 129.31 / 2.88
}

pub(crate) fn taylor_t4_wdot_6pn_s1s1_avg(m_by_m: f64) -> f64 {
    101.9 / (6.4 * m_by_m * m_by_m) + 2.51 / (5.76 * m_by_m) + 13.33 / 5.76
}

pub(crate) fn taylor_t4_wdot_6pn_s1os1o_avg(m_by_m: f64) -> f64 {
    -49.3 / (6.4 * m_by_m * m_by_m) + 197.47 / (5.76 * m_by_m) + 56.45 / 5.76
}

pub(crate) fn taylor_t4_wdot_6pn_qms1s1_avg(m_by_m: f64) -> f64 {
    -6.59 / (2.24 * m_by_m * m_by_m) + 7.3 / (4.8 * m_by_m) - 43. / 4.
}

pub(crate) fn taylor_t4_wdot_6pn_qms1os1o_avg(m_by_m: f64) -> f64 {
    19.77 / (2.24 * m_by_m * m_by_m) - 7.3 / (1.6 * m_by_m) + 129. / 4.
}

pub(crate) fn taylor_t4_wdot_7pn(eta: f64) -> f64 {
    (PI / 12096.0) * fma(731960.0 * eta, eta, fma(717350.0, eta, -13245.0))
}

pub(crate) fn taylor_t4_wdot_10pn_tidal(m_by_m: f64) -> f64 {
    6. * m_by_m * m_by_m * m_by_m * m_by_m * fnms(11., m_by_m, 12.)
}

pub(crate) fn taylor_t4_wdot_12pn_tidal(m_by_m: f64) -> f64 {
    let a = fnms(12263. / 56., m_by_m, 4421. / 56.);
    let b = fma(1893. / 4. * m_by_m, m_by_m, a);
    let c = fnms(661. / 2. * m_by_m * m_by_m, m_by_m, b);
    m_by_m * m_by_m * m_by_m * m_by_m * c
}

// Orbital angular momentum, eq. 4.7 of arXiv:1212.5520 and T1500554.

pub(crate) fn l_2pn(eta: f64) -> f64 {
    1.5 + eta / 6.
}

pub(crate) fn l_3pn_si_coeff_avg(m_by_m: f64) -> f64 {
    -0.75 - 0.25 / m_by_m
}

pub(crate) fn l_3pn_sil_coeff_avg(m_by_m: f64) -> f64 {
    -(1. / 3. + 9. / m_by_m) / 4.
}

pub(crate) fn ldot_3pn_so(m_by_m: f64) -> f64 {
    0.5 + 1.5 / m_by_m
}

// Spin precession.

pub(crate) fn spin_dot_3pn(m_by_m: f64) -> f64 {
    3. / 2. - m_by_m - m_by_m * m_by_m / 2.
}

pub(crate) const SPIN_DOT_4PN_S2_AVG: f64 = 0.5;
pub(crate) const SPIN_DOT_4PN_S2O_AVG: f64 = -1.5;

pub(crate) fn spin_dot_4pn_qmso_avg(m_by_m: f64) -> f64 {
    1.5 * (1. - 1. / m_by_m)
}

pub(crate) fn spin_dot_5pn(m_by_m: f64) -> f64 {
    9. / 8. - m_by_m / 2. + 7. * m_by_m * m_by_m / 12.
        - 7. * m_by_m * m_by_m * m_by_m / 6.
        - m_by_m * m_by_m * m_by_m * m_by_m / 24.
}

fn spin_dot_6pn_s2(m_by_m: f64) -> f64 {
    -1.5 - m_by_m
}

fn spin_dot_6pn_s1n(m_by_m: f64) -> f64 {
    fnms(0.5 * m_by_m, m_by_m, 3.5 - 3. / m_by_m)
}

fn spin_dot_6pn_s2n(m_by_m: f64) -> f64 {
    fma(m_by_m, m_by_m, fma(2., m_by_m, 1.5))
}

fn spin_dot_6pn_s1v(m_by_m: f64) -> f64 {
    fnms(1.5, m_by_m, 3.) - 1.5 / m_by_m
}

fn spin_dot_6pn_s2v(m_by_m: f64) -> f64 {
    1.5 + m_by_m
}

pub(crate) fn spin_dot_6pn_s2_avg(m_by_m: f64) -> f64 {
    fma(
        0.5,
        spin_dot_6pn_s2n(m_by_m) + spin_dot_6pn_s2v(m_by_m),
        spin_dot_6pn_s2(m_by_m),
    )
}

pub(crate) fn spin_dot_6pn_s1o_avg(m_by_m: f64) -> f64 {
    -0.5 * (spin_dot_6pn_s1n(m_by_m) + spin_dot_6pn_s1v(m_by_m))
}

pub(crate) fn spin_dot_6pn_s2o_avg(m_by_m: f64) -> f64 {
    -0.5 * (spin_dot_6pn_s2n(m_by_m) + spin_dot_6pn_s2v(m_by_m))
}

fn spin_dot_6pn_qmsn(m_by_m: f64) -> f64 {
    3. * fnms(0.5 * m_by_m, m_by_m, 0.5 / m_by_m + 1. - m_by_m)
}

fn spin_dot_6pn_qmsv(m_by_m: f64) -> f64 {
    3. * (1. / m_by_m - 1.)
}

pub(crate) fn spin_dot_6pn_qmso_avg(m_by_m: f64) -> f64 {
    -0.5 * (spin_dot_6pn_qmsn(m_by_m) + spin_dot_6pn_qmsv(m_by_m))
}

#[cfg(test)]
mod tests {
    #[test]
    fn folded_log_16() {
        assert_eq!(super::LN_16, 16f64.ln());
        assert_eq!(super::LN_16, 4.0 * std::f64::consts::LN_2);
    }
}
