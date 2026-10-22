// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Roots of a real polynomial as GSL's `gsl_poly_complex_solve`
//! (`poly/zsolve.c`): the eigenvalues of the balanced companion matrix from
//! the Francis double-shift QR algorithm. The roots are returned in the
//! order and packed layout of GSL (`[re0, im0, re1, im1, ...]`), which
//! IMRPhenomXHM's `CrossZeroP5` depends on.
//!
//! The arithmetic follows GSL operation by operation, with the fused
//! multiply-adds of a clang build (see [`fma`]).

use crate::imrphenomx::internals::fma;

/// `GSL_DBL_EPSILON`
const DBL_EPSILON: f64 = f64::EPSILON;

/// `RADIX` of `balance_companion_matrix`.
const RADIX: f64 = 2.0;
const RADIX2: f64 = RADIX * RADIX;

/// Roots of `a[0] + a[1] x + ... + a[n-1] x^(n-1)` (`n = N + 1`, packed as
/// GSL does), or `None` where `gsl_poly_complex_solve` fails (leading
/// coefficient zero or no convergence).
pub(crate) fn complex_solve<const N: usize>(a: &[f64; N]) -> Option<Vec<f64>> {
    let n = N;
    if n < 2 || a[n - 1] == 0.0 {
        return None;
    }
    let nc = n - 1;
    let mut m = companion_matrix(a, nc);
    balance_companion_matrix(&mut m, nc);
    let mut z = vec![0.0; 2 * nc];
    qr_companion(&mut m, nc, &mut z)?;
    Some(z)
}

/// `set_companion_matrix`: row-major `nc x nc` companion matrix.
fn companion_matrix(a: &[f64], nc: usize) -> Vec<f64> {
    let mut m = vec![0.0; nc * nc];
    for i in 1..nc {
        m[i * nc + i - 1] = 1.0;
    }
    for i in 0..nc {
        m[i * nc + nc - 1] = -a[i] / a[nc];
    }
    m
}

/// `balance_companion_matrix`
fn balance_companion_matrix(m: &mut [f64], nc: usize) {
    let at = |i: usize, j: usize| i * nc + j;
    let mut not_converged = true;
    while not_converged {
        not_converged = false;
        for i in 0..nc {
            // column norm, excluding the diagonal
            let mut col_norm = if i != nc - 1 {
                m[at(i + 1, i)].abs()
            } else {
                let mut col_norm = 0.0;
                for j in 0..nc - 1 {
                    col_norm += m[at(j, nc - 1)].abs();
                }
                col_norm
            };
            // row norm, excluding the diagonal
            let row_norm = if i == 0 {
                m[at(0, nc - 1)].abs()
            } else if i == nc - 1 {
                m[at(i, i - 1)].abs()
            } else {
                m[at(i, i - 1)].abs() + m[at(i, nc - 1)].abs()
            };
            if col_norm == 0.0 || row_norm == 0.0 {
                continue;
            }
            let mut g = row_norm / RADIX;
            let mut f = 1.0;
            let s = col_norm + row_norm;
            while col_norm < g {
                f *= RADIX;
                col_norm *= RADIX2;
            }
            g = row_norm * RADIX;
            while col_norm > g {
                f /= RADIX;
                col_norm /= RADIX2;
            }
            if (row_norm + col_norm) < 0.95 * s * f {
                not_converged = true;
                g = 1.0 / f;
                if i == 0 {
                    m[at(0, nc - 1)] *= g;
                } else {
                    m[at(i, i - 1)] *= g;
                    m[at(i, nc - 1)] *= g;
                }
                if i == nc - 1 {
                    for j in 0..nc {
                        m[at(j, i)] *= f;
                    }
                } else {
                    m[at(i + 1, i)] *= f;
                }
            }
        }
    }
}

/// `qr_companion`: eigenvalues of the upper Hessenberg matrix `h`
/// (1-based indices as in GSL's `FMAT`), stored at `zroot[2 (n - 1)]`
/// as they converge at the bottom of the active block.
fn qr_companion(h: &mut [f64], nc: usize, zroot: &mut [f64]) -> Option<()> {
    let idx = |i: usize, j: usize| (i - 1) * nc + (j - 1);
    let mut t = 0.0;
    let mut n = nc;
    let set = |zroot: &mut [f64], k: usize, re: f64, im: f64| {
        zroot[2 * k] = re;
        zroot[2 * k + 1] = im;
    };

    'next_root: loop {
        if n == 0 {
            return Some(());
        }
        let mut iterations = 0;
        loop {
            // next_iteration: look for a small subdiagonal element
            let mut e = n;
            while e >= 2 {
                let a1 = h[idx(e, e - 1)].abs();
                let a2 = h[idx(e - 1, e - 1)].abs();
                let a3 = h[idx(e, e)].abs();
                if a1 <= DBL_EPSILON * (a2 + a3) {
                    break;
                }
                e -= 1;
            }
            let mut x = h[idx(n, n)];
            if e == n {
                // one real root
                set(zroot, n - 1, x + t, 0.0);
                n -= 1;
                continue 'next_root;
            }
            let mut y = h[idx(n - 1, n - 1)];
            let mut w = h[idx(n - 1, n)] * h[idx(n, n - 1)];
            if e == n - 1 {
                let p = (y - x) / 2.0;
                let q = fma(p, p, w);
                y = q.abs().sqrt();
                x += t;
                if q > 0.0 {
                    // two real roots
                    if p < 0.0 {
                        y = -y;
                    }
                    y += p;
                    set(zroot, n - 1, x - w / y, 0.0);
                    set(zroot, n - 2, x + y, 0.0);
                } else {
                    set(zroot, n - 1, x + p, -y);
                    set(zroot, n - 2, x + p, y);
                }
                n -= 2;
                continue 'next_root;
            }
            // no more roots found yet, do another iteration
            if iterations == 120 {
                return None;
            }
            if iterations % 10 == 0 && iterations > 0 {
                // use an exceptional shift
                t += x;
                for i in 1..=n {
                    h[idx(i, i)] -= x;
                }
                let s = h[idx(n, n - 1)].abs() + h[idx(n - 1, n - 2)].abs();
                y = 0.75 * s;
                x = y;
                w = -0.4375 * s * s;
            }
            iterations += 1;

            // GSL keeps p, q, r of the last step of this loop for the first
            // step of the double QR sweep.
            let (mut p, mut q, mut r);
            let mut m = n - 2;
            loop {
                let z = h[idx(m, m)];
                r = x - z;
                let s0 = y - z;
                p = h[idx(m, m + 1)] + fma(r, s0, -w) / h[idx(m + 1, m)];
                q = h[idx(m + 1, m + 1)] - z - r - s0;
                r = h[idx(m + 2, m + 1)];
                let s = p.abs() + q.abs() + r.abs();
                p /= s;
                q /= s;
                r /= s;
                if m == e {
                    break;
                }
                let a1 = h[idx(m, m - 1)].abs();
                let a2 = h[idx(m - 1, m - 1)].abs();
                let a3 = h[idx(m + 1, m + 1)].abs();
                if a1 * (q.abs() + r.abs()) <= DBL_EPSILON * p.abs() * (a2 + a3) {
                    break;
                }
                m -= 1;
            }

            for i in (m + 2)..=n {
                h[idx(i, i - 2)] = 0.0;
            }
            for i in (m + 3)..=n {
                h[idx(i, i - 3)] = 0.0;
            }

            // double QR step
            for k in m..n {
                let notlast = k != n - 1;
                if k != m {
                    p = h[idx(k, k - 1)];
                    q = h[idx(k + 1, k - 1)];
                    r = if notlast { h[idx(k + 2, k - 1)] } else { 0.0 };
                    x = p.abs() + q.abs() + r.abs();
                    if x == 0.0 {
                        continue;
                    }
                    p /= x;
                    q /= x;
                    r /= x;
                }
                let mut s = fma(r, r, fma(p, p, q * q)).sqrt();
                if p < 0.0 {
                    s = -s;
                }
                if k != m {
                    h[idx(k, k - 1)] = -s * x;
                } else if e != m {
                    h[idx(k, k - 1)] *= -1.0;
                }
                p += s;
                x = p / s;
                y = q / s;
                let z = r / s;
                q /= p;
                r /= p;

                // row modifications
                for j in k..=n {
                    p = fma(q, h[idx(k + 1, j)], h[idx(k, j)]);
                    if notlast {
                        p = fma(r, h[idx(k + 2, j)], p);
                        h[idx(k + 2, j)] = fma(-p, z, h[idx(k + 2, j)]);
                    }
                    h[idx(k + 1, j)] = fma(-p, y, h[idx(k + 1, j)]);
                    h[idx(k, j)] = fma(-p, x, h[idx(k, j)]);
                }

                // column modifications
                let jmax = if k + 3 < n { k + 3 } else { n };
                for i in e..=jmax {
                    p = fma(x, h[idx(i, k)], y * h[idx(i, k + 1)]);
                    if notlast {
                        p = fma(z, h[idx(i, k + 2)], p);
                        h[idx(i, k + 2)] = fma(-p, r, h[idx(i, k + 2)]);
                    }
                    h[idx(i, k + 1)] = fma(-p, q, h[idx(i, k + 1)]);
                    h[idx(i, k)] -= p;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::complex_solve;

    #[test]
    fn roots_of_a_quintic() {
        // (x - 1)(x - 2)(x + 3)(x^2 + 1) = x^5 - 6 x^3 + 6 x^2 - 7 x + 6
        let a = [6.0, -7.0, 6.0, -6.0, 0.0, 1.0];
        let z = complex_solve(&a).unwrap();
        let mut roots: Vec<(f64, f64)> = z.chunks(2).map(|c| (c[0], c[1])).collect();
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let expected = [(-3.0, 0.0), (0.0, -1.0), (0.0, 1.0), (1.0, 0.0), (2.0, 0.0)];
        for ((re, im), (ere, eim)) in roots.iter().zip(expected) {
            assert!(
                (re - ere).abs() < 1e-12 && (im - eim).abs() < 1e-12,
                "{roots:?}"
            );
        }
    }
}
