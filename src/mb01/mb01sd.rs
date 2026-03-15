//! MB01SD — Scale matrix by row and/or column factors (SLICOT MB01SD)

use nalgebra::DMatrix;

/// Scaling mode.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01SdJobs {
    /// A := diag(r) * A (row scaling).
    Row,
    /// A := A * diag(c) (column scaling).
    Column,
    /// A := diag(r) * A * diag(c).
    Both,
}

/// Scales A in place. r has length m, c has length n.
pub fn mb01sd(jobs: Mb01SdJobs, a: &mut DMatrix<f64>, r: &[f64], c: &[f64]) -> i32 {
    let m = a.nrows();
    let n = a.ncols();
    if m == 0 || n == 0 {
        return 0;
    }
    if (jobs == Mb01SdJobs::Row || jobs == Mb01SdJobs::Both) && r.len() < m {
        return -5;
    }
    if (jobs == Mb01SdJobs::Column || jobs == Mb01SdJobs::Both) && c.len() < n {
        return -6;
    }
    match jobs {
        Mb01SdJobs::Column => {
            for j in 0..n {
                let cj = c[j];
                for i in 0..m {
                    a[(i, j)] *= cj;
                }
            }
        }
        Mb01SdJobs::Row => {
            for i in 0..m {
                let ri = r[i];
                for j in 0..n {
                    a[(i, j)] *= ri;
                }
            }
        }
        Mb01SdJobs::Both => {
            for j in 0..n {
                let cj = c[j];
                for i in 0..m {
                    a[(i, j)] *= r[i] * cj;
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
    fn test_mb01sd_row() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let r = [2.0, 0.5];
        let c = [1.0, 1.0];
        assert_eq!(mb01sd(Mb01SdJobs::Row, &mut a, &r, &c), 0);
        assert_eq!(a[(0, 0)], 2.0);
        assert_eq!(a[(1, 1)], 2.0);
    }
}
