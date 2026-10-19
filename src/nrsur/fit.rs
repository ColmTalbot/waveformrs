// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
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

/// Number of factors per term kept by the compact representation.
const COMPACT: usize = 4;
/// Index of an entry of the power table that is exactly one (`x_1^0`).
const ONE: u8 = 1;

/// A scalar fit with its basis functions as indices into a [`PowerTable`].
///
/// Order-zero factors are exactly one (`ipow(x, 0) = 1.0`) and multiplying
/// by them is exact, so each term keeps only its other factors, in LAL's
/// order, padded with ones. Terms have at most four such factors except for
/// two terms of the NRSur7dq4 data; fits with longer terms keep all seven.
#[derive(Clone, Debug, Default)]
pub(crate) struct Fit {
    coefs: Vec<f64>,
    /// The non-trivial factors of each term (`7 * order + j`).
    compact: Vec<[u8; COMPACT]>,
    /// All factors, if some term has more than [`COMPACT`] non-trivial ones.
    full: Vec<[u8; 7]>,
}

impl Fit {
    /// A fit from the `coefs` and the `n x 7` (row-major) basis-function
    /// orders of the data file. Returns `None` if an order exceeds the
    /// power table (LAL would read out of bounds).
    pub(crate) fn new(coefs: Vec<f64>, orders: &[i64]) -> Option<Self> {
        if orders.len() != 7 * coefs.len() {
            return None;
        }
        let mut full = Vec::with_capacity(coefs.len());
        let mut compact = Vec::with_capacity(coefs.len());
        let mut fits_compact = true;
        for row in orders.chunks_exact(7) {
            let mut all = [0u8; 7];
            let mut short = [ONE; COMPACT];
            let mut n = 0;
            for (j, (&order, out)) in row.iter().zip(&mut all).enumerate() {
                let index = 7 * order + j as i64;
                if order < 0 || index >= N_POWERS as i64 {
                    return None;
                }
                *out = index as u8;
                if order > 0 {
                    if n < COMPACT {
                        short[n] = index as u8;
                    }
                    n += 1;
                }
            }
            fits_compact &= n <= COMPACT;
            full.push(all);
            compact.push(short);
        }
        if fits_compact {
            full = Vec::new();
        } else {
            compact = Vec::new();
        }
        Some(Self {
            coefs,
            compact,
            full,
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.coefs.len()
    }

    /// Continue the sum of a compact fit from term `start` with `res`.
    #[inline(always)]
    fn eval_compact_from(&self, start: usize, mut res: f64, x: &[f64; TABLE_SIZE]) -> f64 {
        for (&c, p) in self.coefs[start..].iter().zip(&self.compact[start..]) {
            res = fma(c, compact_product(x, p), res);
        }
        res
    }

    /// `NRSur7dq4_eval_fit` with the powers of the parameters tabulated.
    #[inline]
    pub(crate) fn eval(&self, table: &PowerTable) -> f64 {
        let x = &table.0;
        if self.full.is_empty() {
            return self.eval_compact_from(0, 0.0, x);
        }
        let mut res = 0.0;
        for (&c, p) in self.coefs.iter().zip(&self.full) {
            let mut prod = x[p[0] as usize & INDEX_MASK];
            for &i in &p[1..] {
                prod *= x[i as usize & INDEX_MASK];
            }
            res = fma(c, prod, res);
        }
        res
    }
}

/// The product of the factors of a compact term, in order.
#[inline(always)]
fn compact_product(x: &[f64; TABLE_SIZE], p: &[u8; COMPACT]) -> f64 {
    let mut prod = x[p[0] as usize & INDEX_MASK];
    prod *= x[p[1] as usize & INDEX_MASK];
    prod *= x[p[2] as usize & INDEX_MASK];
    prod *= x[p[3] as usize & INDEX_MASK];
    prod
}

/// Size of the power table in memory: [`N_POWERS`] rounded up to a power
/// of two, so that masked indices need no bounds checks.
const TABLE_SIZE: usize = 32;
const INDEX_MASK: usize = TABLE_SIZE - 1;

/// The powers `x_powers[i] = ipow(p[i % 7], i / 7)` of the 7 fit parameters
/// (`x_powers` in `NRSur7dq4_eval_fit`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct PowerTable(pub(crate) [f64; TABLE_SIZE]);

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
        let mut t = [0.0; TABLE_SIZE];
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
        // with all seven factors (long terms) and with the compact factors
        let fit = Fit::new(coefs.clone(), &orders).unwrap();
        assert!(fit.compact.is_empty());
        let short_coefs = vec![0.3, -1.2, 0.7, 1.1];
        let short_orders = vec![
            0, 0, 0, 0, 0, 0, 0, //
            1, 0, 2, 1, 0, 0, 1, //
            2, 1, 0, 0, 1, 2, 0, //
            0, 0, 0, 2, 0, 0, 0,
        ];
        let short = Fit::new(short_coefs.clone(), &short_orders).unwrap();
        assert!(short.full.is_empty());
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
            assert_eq!(
                short.eval(&table).to_bits(),
                lal_eval_fit(&short_coefs, &short_orders, &x).to_bits()
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
