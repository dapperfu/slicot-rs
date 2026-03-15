//! DLACPY_SLC — Copy all or part of a matrix (SLICOT/LAPACK auxiliary)
//
// Copies A into B for the specified triangular part or full matrix.
// Real version; UPLO 'U' = upper, 'L' = lower, else full.

use nalgebra::DMatrix;

/// Which part of the matrix to copy.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum DlacpyUplo {
    /// Upper triangular part.
    Upper,
    /// Lower triangular part.
    Lower,
    /// Full matrix.
    All,
}

/// Copies the specified part of A into B. B must be M×N (same shape as A).
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn dlacpy_slc(uplo: DlacpyUplo, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    let m = a.nrows();
    let n = a.ncols();
    if b.nrows() != m || b.ncols() != n {
        return -6;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    match uplo {
        DlacpyUplo::Upper => {
            for j in 0..n {
                for i in 0..m.min(j + 1) {
                    b[(i, j)] = a[(i, j)];
                }
            }
        }
        DlacpyUplo::Lower => {
            for j in 0..n {
                for i in j..m {
                    b[(i, j)] = a[(i, j)];
                }
            }
        }
        DlacpyUplo::All => {
            for j in 0..n {
                for i in 0..m {
                    b[(i, j)] = a[(i, j)];
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
    fn test_dlacpy_slc_all() {
        let a = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let mut b = DMatrix::zeros(2, 3);
        assert_eq!(dlacpy_slc(DlacpyUplo::All, &a, &mut b), 0);
        assert_eq!(a, b);
    }

    #[test]
    fn test_dlacpy_slc_upper() {
        let a = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let mut b = DMatrix::zeros(3, 3);
        assert_eq!(dlacpy_slc(DlacpyUplo::Upper, &a, &mut b), 0);
        assert_eq!(b[(0, 0)], 1.0);
        assert_eq!(b[(0, 1)], 2.0);
        assert_eq!(b[(1, 1)], 5.0);
        assert_eq!(b[(0, 2)], 3.0);
        assert_eq!(b[(1, 2)], 6.0);
        assert_eq!(b[(2, 2)], 9.0);
    }
}
