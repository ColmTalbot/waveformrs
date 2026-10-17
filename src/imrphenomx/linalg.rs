//! Dense linear solve for the small collocation systems of IMRPhenomX.
//!
//! LALSimulation solves these with `gsl_linalg_LU_decomp` followed by
//! `gsl_linalg_LU_solve`. For matrices this small GSL (>= 2.7) uses the
//! unblocked, partially pivoted algorithm `LU_decomp_L2` (LAPACK `DGETF2`),
//! and the solve applies the row permutation followed by unit lower and
//! upper triangular substitution. The operations below follow the GSL
//! reference implementation (including scaling the pivot column by the
//! reciprocal of the pivot) so that results agree to rounding.

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
                row[k] += pivot_row[k] * tmp;
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
        let mut tmp = b[i];
        for (aij, bj) in a[i][..i].iter().zip(&b[..i]) {
            tmp -= aij * bj;
        }
        b[i] = tmp;
    }

    // Back substitution with the upper triangle (cblas_dtrsv).
    b[N - 1] /= a[N - 1][N - 1];
    for i in (0..N - 1).rev() {
        let mut tmp = b[i];
        for (aij, bj) in a[i][i + 1..].iter().zip(&b[i + 1..]) {
            tmp -= aij * bj;
        }
        b[i] = tmp / a[i][i];
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
}
