// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Dense linear solve for the small collocation systems of IMRPhenomX.
//!
//! LALSimulation solves these with `gsl_linalg_LU_decomp` followed by
//! `gsl_linalg_LU_solve`. For matrices this small GSL (>= 2.7) uses the
//! unblocked, partially pivoted algorithm `LU_decomp_L2` (LAPACK `DGETF2`),
//! and the solve applies the row permutation followed by unit lower and
//! upper triangular substitution. The operations below follow the GSL
//! reference implementation (including scaling the pivot column by the
//! reciprocal of the pivot).
//!
//! GSL delegates the level-2 BLAS operations to the CBLAS it is linked
//! against (OpenBLAS in the conda-forge builds). Its `dger` rank-one update
//! accumulates `a_ik += (-l_ij) u_jk` with fused multiply-adds, and its
//! `dtrsv` subtracts a dot product that is accumulated left to right with
//! fused multiply-adds starting from zero. Reproducing these (with
//! [`fma`], fused on the arm64 targets whose LAL builds fuse) makes the
//! solution agree with LALSimulation to the last bit (verified on random
//! 4x4 to 8x8 systems against `gsl_linalg_LU_solve`). On x86-64 OpenBLAS
//! picks its kernels at run time: the unfused operations here match its
//! kernels for CPUs without FMA (verified with the Nehalem kernels that it
//! selects under Rosetta), while on Haswell and later its kernels fuse and
//! the solution agrees to rounding.

use crate::spintaylor::fp::fma;

/// Solve `a x = b` for an `N x N` system using LU decomposition with
/// partial pivoting, replicating `gsl_linalg_LU_decomp` + `gsl_linalg_LU_solve`.
pub(crate) fn lu_solve<const N: usize>(mut a: [[f64; N]; N], mut b: [f64; N]) -> [f64; N] {
    let mut pivots = [0usize; N];

    // LU_decomp_L2
    for j in 0..N {
        // gsl_blas_idamax: first index of the maximum absolute value.
        let mut j_pivot = j;
        let mut max = a[j][j].abs();
        for (i, row) in a.iter().enumerate().skip(j + 1) {
            if row[j].abs() > max {
                max = row[j].abs();
                j_pivot = i;
            }
        }
        pivots[j] = j_pivot;

        if a[j_pivot][j] != 0.0 {
            if j_pivot != j {
                a.swap(j, j_pivot);
            }
            let ajj = a[j][j];
            if ajj.abs() >= f64::MIN_POSITIVE {
                let scale = 1.0 / ajj;
                for row in a.iter_mut().skip(j + 1) {
                    row[j] *= scale;
                }
            } else {
                for row in a.iter_mut().skip(j + 1) {
                    row[j] /= ajj;
                }
            }
        }

        // Rank-one update of the trailing block, gsl_blas_dger with alpha = -1.
        let pivot_row = a[j];
        for row in a.iter_mut().skip(j + 1) {
            let tmp = -row[j];
            for k in (j + 1)..N {
                row[k] = fma(tmp, pivot_row[k], row[k]);
            }
        }
    }

    // gsl_permute_vector with the permutation built from the pivots is
    // equivalent to applying the row swaps in order.
    for (j, &pivot) in pivots.iter().enumerate() {
        b.swap(j, pivot);
    }

    // Forward substitution with the unit lower triangle (cblas_dtrsv).
    for i in 1..N {
        let dot = a[i][..i]
            .iter()
            .zip(&b[..i])
            .fold(0.0, |acc, (aij, bj)| fma(*aij, *bj, acc));
        b[i] -= dot;
    }

    // Back substitution with the upper triangle (cblas_dtrsv).
    b[N - 1] /= a[N - 1][N - 1];
    for i in (0..N - 1).rev() {
        let dot = a[i][i + 1..]
            .iter()
            .zip(&b[i + 1..])
            .fold(0.0, |acc, (aij, bj)| fma(*aij, *bj, acc));
        b[i] = (b[i] - dot) / a[i][i];
    }
    b
}

#[cfg(test)]
mod tests {
    use super::lu_solve;

    #[test]
    fn solves_small_system() {
        let a = [
            [2.0, 1.0, -1.0, 0.5],
            [-3.0, -1.0, 2.0, 1.0],
            [-2.0, 1.0, 2.0, -4.0],
            [1.0, 7.0, 3.0, 2.0],
        ];
        let x = [1.0, -2.0, 3.0, 0.25];
        let b = a.map(|row| row.iter().zip(x).map(|(aij, xj)| aij * xj).sum());
        let solution = lu_solve(a, b);
        for (s, x) in solution.iter().zip(x) {
            assert!((s - x).abs() < 1e-14, "{s} != {x}");
        }
    }

    #[test]
    fn matches_gsl_openblas() {
        // gsl_linalg_LU_decomp + gsl_linalg_LU_solve (GSL 2.8, OpenBLAS
        // 0.3.34) for a collocation-like system.
        let a = [
            [1.0, 0.5, 0.25, 0.125, 0.0625],
            [1.0, 0.25, 0.0625, 0.015625, 0.00390625],
            [1.0, 2.0, 4.0, 8.0, 16.0],
            [1.0, 3.0, 9.0, 27.0, 81.0],
            [1.0, -1.5, 2.25, -3.375, 5.0625],
        ];
        let b = [0.3, -1.7, 2.9, 4.1, -0.6];
        let x = lu_solve(a, b);
        // (osx-arm64; the solution differs between OpenBLAS's kernels by
        // ~1e-14, a condition number of ~100 times rounding)
        let expected = [
            -3.873704390847249,
            8.711915069057929,
            0.6889095031952177,
            -3.2225520511234795,
            0.7734157905586476,
        ];
        for (x, expected) in x.iter().zip(expected) {
            assert!(
                (x - expected).abs() <= 1e-14 * expected.abs(),
                "{x} != {expected}"
            );
        }
    }
}
