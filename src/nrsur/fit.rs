//! Parameter-space fits of the precessing NR surrogates.
//!
//! A fit is `sum_i c_i prod_j B_j(k_ij; x_j)` where `B_j` are monomials of
//! order `k_ij` in the (transformed) parameters `x_j`. LALSimulation
//! (`NRSur7dq4_eval_fit`) tabulates the 22 powers `x_j^k` (`k <= 3` for the
//! mass ratio, `k <= 2` for the spins) for every fit; here the table is
//! built once per parameter point ([`PowerTable`]) and shared by all fits
//! evaluated there, and the basis-function orders are turned into indices
//! into that table when the data are loaded. The arithmetic (the product
//! order and the fused `res += c * prod` of clang's arm64 build) is the
//! same, so the fits agree with LAL bit for bit.

use crate::spintaylor::fp::{fma, fnms};

/// Number of entries of the power table: 4 powers of the mass-ratio
/// parameter, 3 of each of the 6 spin parameters.
pub(crate) const N_POWERS: usize = 22;

/// A scalar fit with its basis functions as indices into a [`PowerTable`].
#[derive(Clone, Debug, Default)]
pub(crate) struct Fit {
    coefs: Vec<f64>,
    /// `7 * order + j` for the 7 parameters `j` of each term.
    powers: Vec<[u8; 7]>,
}

impl Fit {
    /// A fit from the `coefs` and the `n x 7` (row-major) basis-function
    /// orders of the data file. Returns `None` if an order exceeds the
    /// power table (LAL would read out of bounds).
    pub(crate) fn new(coefs: Vec<f64>, orders: &[i64]) -> Option<Self> {
        if orders.len() != 7 * coefs.len() {
            return None;
        }
        let powers = orders
            .chunks_exact(7)
            .map(|row| {
                let mut out = [0u8; 7];
                for (j, (&order, out)) in row.iter().zip(&mut out).enumerate() {
                    let index = 7 * order + j as i64;
                    if order < 0 || index >= N_POWERS as i64 {
                        return None;
                    }
                    *out = index as u8;
                }
                Some(out)
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self { coefs, powers })
    }

    pub(crate) fn len(&self) -> usize {
        self.coefs.len()
    }

    /// `NRSur7dq4_eval_fit` with the powers of the parameters tabulated.
    #[inline]
    pub(crate) fn eval(&self, table: &PowerTable) -> f64 {
        let x = &table.0;
        let mut res = 0.0;
        for (&c, p) in self.coefs.iter().zip(&self.powers) {
            let mut prod = x[p[0] as usize];
            prod *= x[p[1] as usize];
            prod *= x[p[2] as usize];
            prod *= x[p[3] as usize];
            prod *= x[p[4] as usize];
            prod *= x[p[5] as usize];
            prod *= x[p[6] as usize];
            res = fma(c, prod, res);
        }
        res
    }
}

/// The powers `x_powers[i] = ipow(p[i % 7], i / 7)` of the 7 fit parameters
/// (`x_powers` in `NRSur7dq4_eval_fit`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct PowerTable(pub(crate) [f64; N_POWERS]);

/// The parts of the NRSur7dq4 fit parameters that depend only on the mass
/// ratio: `q_fit` and its powers, and the effective-spin coefficients.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MassRatioParameters {
    q: f64,
    one_plus_q: f64,
    /// `38 eta / 113`.
    chi_hat_slope: f64,
    /// `1 - 76 eta / 113`.
    chi_hat_denominator: f64,
    q_fit: [f64; 4],
}

// These rescale log(q) from [-0.01, log(4.01)] to [-1, 1].
const NRSUR7DQ4_Q_FIT_SLOPE: f64 = 1.4298059216576398;
const NRSUR7DQ4_Q_FIT_OFFSET: f64 = -0.9857019407834238;

impl MassRatioParameters {
    /// NRSur7dq4: `q_fit = offset + slope log(q)` and the `eta`-dependent
    /// coefficients of `NRSur7dq4_effective_spins`.
    pub(crate) fn nrsur7dq4(q: f64) -> Self {
        let q_fit = fma(NRSUR7DQ4_Q_FIT_SLOPE, q.ln(), NRSUR7DQ4_Q_FIT_OFFSET);
        let one_plus_q = 1.0 + q;
        let eta = q / one_plus_q / one_plus_q;
        Self {
            q,
            one_plus_q,
            chi_hat_slope: 38.0 * eta / 113.0,
            chi_hat_denominator: 1.0 - 76.0 * eta / 113.0,
            q_fit: [1.0, q_fit, q_fit * q_fit, q_fit * q_fit * q_fit],
        }
    }

    /// The power table at the point `x = [q, chi1x, chi1y, chi1z, chi2x,
    /// chi2y, chi2z]` (`x[0]` is ignored: the mass ratio is fixed). The fit
    /// parameters are `[log q, chi1x, chi1y, chiHat, chi2x, chi2y, chi_a]`.
    #[inline]
    pub(crate) fn table(&self, x: &[f64; 7]) -> PowerTable {
        let (chi1z, chi2z) = (x[3], x[6]);
        let chi_wt_avg = fma(self.q, chi1z, chi2z) / self.one_plus_q;
        let chi_hat =
            fnms(self.chi_hat_slope, chi1z + chi2z, chi_wt_avg) / self.chi_hat_denominator;
        let chi_a = (chi1z - chi2z) / 2.0;
        let p = [x[1], x[2], chi_hat, x[4], x[5], chi_a];
        let mut t = [0.0; N_POWERS];
        t[0] = self.q_fit[0];
        t[7] = self.q_fit[1];
        t[14] = self.q_fit[2];
        t[21] = self.q_fit[3];
        for (j, &p) in p.iter().enumerate() {
            t[1 + j] = 1.0;
            t[8 + j] = p;
            t[15 + j] = p * p;
        }
        PowerTable(t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Direct transcription of `NRSur7dq4_eval_fit` (recomputing everything).
    fn lal_eval_fit(coefs: &[f64], orders: &[i64], x: &[f64; 7]) -> f64 {
        fn ipow(base: f64, exponent: i64) -> f64 {
            if exponent == 0 {
                return 1.0;
            }
            let mut res = base;
            for _ in 1..exponent {
                res *= base;
            }
            res
        }
        let q = x[0];
        let eta = q / (1. + q) / (1. + q);
        let chi_wt_avg = fma(q, x[3], x[6]) / (1.0 + q);
        let chi_hat = fnms(38. * eta / 113., x[3] + x[6], chi_wt_avg) / (1. - 76. * eta / 113.);
        let chi_a = (x[3] - x[6]) / 2.;
        let fit_params = [q.ln(), x[1], x[2], chi_hat, x[4], x[5], chi_a];
        let q_fit = fma(NRSUR7DQ4_Q_FIT_SLOPE, fit_params[0], NRSUR7DQ4_Q_FIT_OFFSET);
        let mut x_powers = [0.0; 22];
        for (i, xp) in x_powers.iter_mut().enumerate() {
            let base = if i % 7 == 0 { q_fit } else { fit_params[i % 7] };
            *xp = ipow(base, (i / 7) as i64);
        }
        let mut res = 0.0;
        for (i, &c) in coefs.iter().enumerate() {
            let mut prod = x_powers[7 * orders[7 * i] as usize];
            for j in 1..7 {
                prod *= x_powers[7 * orders[7 * i + j] as usize + j];
            }
            res = fma(c, prod, res);
        }
        res
    }

    #[test]
    fn table_matches_lal_power_loop() {
        let coefs = vec![0.3, -1.2, 2.5, 0.7, -0.01];
        let orders = vec![
            0, 0, 0, 0, 0, 0, 0, //
            1, 0, 2, 1, 0, 0, 1, //
            3, 2, 2, 2, 2, 2, 2, //
            2, 1, 0, 0, 1, 2, 0, //
            1, 1, 1, 1, 1, 1, 1,
        ];
        let fit = Fit::new(coefs.clone(), &orders).unwrap();
        for x in [
            [1.0, 0.1, -0.2, 0.3, 0.4, -0.5, 0.6],
            [3.7, -0.8, 0.05, -0.33, 0.0, 0.71, -0.2],
            [6.0, 0.9, 0.1, 0.2, -0.3, 0.4, 0.5],
        ] {
            let table = MassRatioParameters::nrsur7dq4(x[0]).table(&x);
            assert_eq!(
                fit.eval(&table).to_bits(),
                lal_eval_fit(&coefs, &orders, &x).to_bits()
            );
        }
    }

    #[test]
    fn rejects_orders_beyond_the_table() {
        assert!(Fit::new(vec![1.0], &[4, 0, 0, 0, 0, 0, 0]).is_none());
        assert!(Fit::new(vec![1.0], &[0, 3, 0, 0, 0, 0, 0]).is_none());
        assert!(Fit::new(vec![1.0], &[3, 2, 2, 2, 2, 2, 2]).is_some());
    }
}
