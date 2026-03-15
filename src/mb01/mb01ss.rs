//! MB01SS — Scale symmetric matrix by diag(D) or inv(diag(D)) (SLICOT MB01SS)

use nalgebra::DMatrix;

/// Scaling mode.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01SsJobs {
    /// A := diag(D)*A*diag(D).
    Scale,
    /// A := inv(diag(D))*A*inv(diag(D)).
    Inverse,
}

/// Which triangle is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01SsUplo {
    Upper,
    Lower,
}

/// Scales the symmetric matrix A (only the specified triangle) in place. D has length n.
pub fn mb01ss(jobs: Mb01SsJobs, uplo: Mb01SsUplo, a: &mut DMatrix<f64>, d: &[f64]) -> i32 {
    let n = a.nrows();
    if a.ncols() != n || d.len() < n {
        return -1;
    }
    if n == 0 {
        return 0;
    }
    match (jobs, uplo) {
        (Mb01SsJobs::Scale, Mb01SsUplo::Upper) => {
            for j in 0..n {
                let dj = d[j];
                for i in 0..=j {
                    a[(i, j)] *= dj * d[i];
                }
            }
        }
        (Mb01SsJobs::Scale, Mb01SsUplo::Lower) => {
            for j in 0..n {
                let dj = d[j];
                for i in j..n {
                    a[(i, j)] *= dj * d[i];
                }
            }
        }
        (Mb01SsJobs::Inverse, Mb01SsUplo::Upper) => {
            for j in 0..n {
                for i in 0..=j {
                    let denom = d[i] * d[j];
                    if denom != 0.0 {
                        a[(i, j)] /= denom;
                    }
                }
            }
        }
        (Mb01SsJobs::Inverse, Mb01SsUplo::Lower) => {
            for j in 0..n {
                for i in j..n {
                    let denom = d[i] * d[j];
                    if denom != 0.0 {
                        a[(i, j)] /= denom;
                    }
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
    fn test_mb01ss_scale_upper() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
        let d = [2.0, 3.0];
        assert_eq!(mb01ss(Mb01SsJobs::Scale, Mb01SsUplo::Upper, &mut a, &d), 0);
        assert_eq!(a[(0, 0)], 4.0);
        assert_eq!(a[(0, 1)], 12.0);
        assert_eq!(a[(1, 1)], 27.0);
    }
}
