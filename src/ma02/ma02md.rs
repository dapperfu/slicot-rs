//! MA02MD — Compute norms of a real skew-symmetric matrix (SLICOT MA02MD)
//
// One norm, Frobenius, infinity, or max absolute value. For skew-symmetric, infinity = one norm.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02MdNorm {
    One,
    Frobenius,
    Infinity,
    MaxAbs,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02MdUplo {
    Upper,
    Lower,
}

/// Returns the specified norm of the skew-symmetric matrix A (only one triangle stored).
pub fn ma02md(norm: Ma02MdNorm, uplo: Ma02MdUplo, a: &DMatrix<f64>, dwork: &mut [f64]) -> f64 {
    let n = a.nrows();
    if n == 0 || a.ncols() != n {
        return 0.0;
    }

    match norm {
        Ma02MdNorm::MaxAbs => {
            let mut m = 0.0_f64;
            match uplo {
                Ma02MdUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            m = m.max(a[(i, j)].abs());
                        }
                    }
                }
                Ma02MdUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            m = m.max(a[(i, j)].abs());
                        }
                    }
                }
            }
            m
        }
        Ma02MdNorm::Frobenius => {
            let mut sum = 0.0;
            match uplo {
                Ma02MdUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            let x = a[(i, j)];
                            sum += x * x;
                        }
                    }
                }
                Ma02MdUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            let x = a[(i, j)];
                            sum += x * x;
                        }
                    }
                }
            }
            (2.0 * sum).sqrt()
        }
        Ma02MdNorm::One | Ma02MdNorm::Infinity => {
            if dwork.len() < n {
                return 0.0;
            }
            for i in 0..n {
                dwork[i] = 0.0;
            }
            match uplo {
                Ma02MdUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            let x = a[(i, j)].abs();
                            dwork[i] += x;
                            dwork[j] += x;
                        }
                    }
                }
                Ma02MdUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            let x = a[(i, j)].abs();
                            dwork[i] += x;
                            dwork[j] += x;
                        }
                    }
                }
            }
            dwork.iter().fold(0.0_f64, |a, &b| a.max(b))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02md_zero() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut dwork = vec![0.0; 1];
        assert_eq!(ma02md(Ma02MdNorm::One, Ma02MdUplo::Upper, &a, &mut dwork), 0.0);
    }

    #[test]
    fn test_ma02md_frobenius_2x2() {
        // Skew-symmetric: [0 1; -1 0], Frobenius = sqrt(1+1) = sqrt(2)
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, -1.0, 0.0]);
        let mut dwork = vec![0.0; 2];
        let f = ma02md(Ma02MdNorm::Frobenius, Ma02MdUplo::Upper, &a, &mut dwork);
        assert!((f - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_ma02md_max_abs() {
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 3.0, -3.0, 0.0]);
        let mut dwork = vec![0.0; 2];
        assert!((ma02md(Ma02MdNorm::MaxAbs, Ma02MdUplo::Upper, &a, &mut dwork) - 3.0).abs() < 1e-10);
    }
}
