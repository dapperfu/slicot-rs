//! MA02AD — Transpose all or part of a matrix into another (SLICOT MA02AD)

use nalgebra::DMatrix;

/// Part of the matrix to transpose.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02AdJob {
    /// Upper triangular part.
    Upper,
    /// Lower triangular part.
    Lower,
    /// Full matrix.
    All,
}

/// Transposes the specified part of A into B. B must be N×M (so B = A' in the specified part).
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn ma02ad(job: Ma02AdJob, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    let m = a.nrows();
    let n = a.ncols();
    if b.nrows() != n || b.ncols() != m {
        return -6;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    match job {
        Ma02AdJob::Upper => {
            for j in 0..n {
                for i in 0..m.min(j + 1) {
                    b[(j, i)] = a[(i, j)];
                }
            }
        }
        Ma02AdJob::Lower => {
            for j in 0..n {
                for i in j..m {
                    b[(j, i)] = a[(i, j)];
                }
            }
        }
        Ma02AdJob::All => {
            for j in 0..n {
                for i in 0..m {
                    b[(j, i)] = a[(i, j)];
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
    fn test_ma02ad_all() {
        let a = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let mut b = DMatrix::zeros(3, 2);
        assert_eq!(ma02ad(Ma02AdJob::All, &a, &mut b), 0);
        assert_eq!(b[(0, 0)], 1.0);
        assert_eq!(b[(1, 0)], 2.0);
        assert_eq!(b[(0, 1)], 4.0);
        assert_eq!(b[(2, 1)], 6.0);
    }
}
