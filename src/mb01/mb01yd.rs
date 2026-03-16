//! MB01YD — C := alpha*op(A)*op(A)' + beta*C, C symmetric, A with L codiagonals (SLICOT MB01YD)

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01YdUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01YdTrans {
    NoTrans,
    Trans,
}

/// Symmetric rank-k update: C := alpha*op(A)*op(A)' + beta*C.
/// A has L nonzero subdiagonals (UPLO='U') or superdiagonals (UPLO='L'); only that band is used.
pub fn mb01yd(
    uplo: Mb01YdUplo,
    trans: Mb01YdTrans,
    n: usize,
    k: usize,
    l: usize,
    alpha: f64,
    beta: f64,
    a: &DMatrix<f64>,
    c: &mut DMatrix<f64>,
) -> i32 {
    if a.nrows() != match trans {
        Mb01YdTrans::NoTrans => n,
        Mb01YdTrans::Trans => k,
    } || a.ncols() != match trans {
        Mb01YdTrans::NoTrans => k,
        Mb01YdTrans::Trans => n,
    } {
        return -8;
    }
    if c.nrows() != n || c.ncols() != n {
        return -10;
    }
    if n == 0 {
        return 0;
    }

    // C := beta*C (only the stored triangle)
    if beta == 0.0 {
        for j in 0..n {
            for i in 0..n {
                if (uplo == Mb01YdUplo::Upper && i <= j) || (uplo == Mb01YdUplo::Lower && i >= j) {
                    c[(i, j)] = 0.0;
                }
            }
        }
    } else if beta != 1.0 {
        for j in 0..n {
            for i in 0..n {
                if (uplo == Mb01YdUplo::Upper && i <= j) || (uplo == Mb01YdUplo::Lower && i >= j) {
                    c[(i, j)] *= beta;
                }
            }
        }
    }

    if alpha == 0.0 || k == 0 {
        return 0;
    }

    // Band mask: for UPLO='U', A has upper triangle + L subdiagonals => (i,j) used iff i <= j+L
    // For UPLO='L', A has lower triangle + L superdiagonals => (i,j) used iff i >= j.saturating_sub(L)
    let in_band = |i: usize, j: usize, _nr: usize, _nc: usize| -> bool {
        match uplo {
            Mb01YdUplo::Upper => i <= j + l,
            Mb01YdUplo::Lower => i + l >= j,
        }
    };

    match trans {
        Mb01YdTrans::NoTrans => {
            // C += alpha * A * A'; A is n×k
            for p in 0..n {
                for i in 0..n {
                    if (uplo == Mb01YdUplo::Upper && i > p) || (uplo == Mb01YdUplo::Lower && i < p) {
                        continue;
                    }
                    let mut sum = 0.0;
                    for j in 0..k {
                        if in_band(i, j, n, k) && in_band(p, j, n, k) {
                            sum += a[(i, j)] * a[(p, j)];
                        }
                    }
                    c[(i, p)] += alpha * sum;
                }
            }
        }
        Mb01YdTrans::Trans => {
            // C += alpha * A' * A; A is k×n
            for p in 0..n {
                for i in 0..n {
                    if (uplo == Mb01YdUplo::Upper && i > p) || (uplo == Mb01YdUplo::Lower && i < p) {
                        continue;
                    }
                    let mut sum = 0.0;
                    for j in 0..k {
                        if in_band(j, i, k, n) && in_band(j, p, k, n) {
                            sum += a[(j, i)] * a[(j, p)];
                        }
                    }
                    c[(i, p)] += alpha * sum;
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01yd_upper_notrans_full_band() {
        // C := 1*A*A' + 0*C, n=2, k=2, L=1 (full 2x2), upper C
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
        let mut c = DMatrix::zeros(2, 2);
        assert_eq!(
            mb01yd(
                Mb01YdUplo::Upper,
                Mb01YdTrans::NoTrans,
                2,
                2,
                1,
                1.0,
                0.0,
                &a,
                &mut c,
            ),
            0
        );
        // A*A' = [1 2; 0 3] * [1 0; 2 3] = [5 6; 6 9]
        assert!((c[(0, 0)] - 5.0).abs() < 1e-10);
        assert!((c[(0, 1)] - 6.0).abs() < 1e-10);
        assert!((c[(1, 1)] - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_mb01yd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut c = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(
            mb01yd(
                Mb01YdUplo::Upper,
                Mb01YdTrans::NoTrans,
                0,
                0,
                0,
                1.0,
                0.0,
                &a,
                &mut c,
            ),
            0
        );
    }
}
