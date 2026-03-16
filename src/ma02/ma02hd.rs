//! MA02HD — Check if A = DIAG*I (identity-like: ones on diagonal, zeros elsewhere) (SLICOT MA02HD)

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02HdJob {
    Upper,
    Lower,
    All,
}

/// Returns true if A = DIAG*I (I has ones on diagonal, zeros elsewhere); false otherwise.
/// If min(M,N) = 0, returns false.
pub fn ma02hd(job: Ma02HdJob, a: &DMatrix<f64>, diag: f64) -> bool {
    let m = a.nrows();
    let n = a.ncols();
    if m == 0 || n == 0 {
        return false;
    }
    let eps = 1e-15_f64 * (1.0_f64).max(diag.abs()).max(1.0);
    match job {
        Ma02HdJob::Upper => {
            for j in 0..n {
                for i in 0..m {
                    if i <= j {
                        let expected = if i == j { diag } else { 0.0 };
                        if (a[(i, j)] - expected).abs() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        Ma02HdJob::Lower => {
            for j in 0..n {
                for i in 0..m {
                    if i >= j {
                        let expected = if i == j { diag } else { 0.0 };
                        if (a[(i, j)] - expected).abs() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        Ma02HdJob::All => {
            for i in 0..m {
                for j in 0..n {
                    let expected = if i == j { diag } else { 0.0 };
                    if (a[(i, j)] - expected).abs() > eps {
                        return false;
                    }
                }
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02hd_identity() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        assert!(ma02hd(Ma02HdJob::All, &a, 1.0));
        assert!(!ma02hd(Ma02HdJob::All, &a, 2.0));
    }

    #[test]
    fn test_ma02hd_scalar_multiple() {
        let a = DMatrix::from_row_slice(2, 2, &[3.0, 0.0, 0.0, 3.0]);
        assert!(ma02hd(Ma02HdJob::All, &a, 3.0));
    }

    #[test]
    fn test_ma02hd_false_off_diagonal() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 0.0, 1.0]);
        assert!(!ma02hd(Ma02HdJob::All, &a, 1.0));
    }

    #[test]
    fn test_ma02hd_zero_dim() {
        let a = DMatrix::<f64>::zeros(0, 0);
        assert!(!ma02hd(Ma02HdJob::All, &a, 1.0));
    }
}
