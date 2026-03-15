//! MA02MZ — Compute norms of a complex skew-Hermitian matrix (SLICOT MA02MZ)
//
// One, Frobenius, infinity, or max absolute value. Diagonal real parts assumed zero; only one triangle used.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02MzNorm {
    One,
    Frobenius,
    Infinity,
    MaxAbs,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02MzUplo {
    Upper,
    Lower,
}

fn cabs(re: f64, im: f64) -> f64 {
    (re * re + im * im).sqrt()
}

/// Returns the specified norm. a_re, a_im column-major LDA×N. dwork length >= N for One/Infinity.
pub fn ma02mz(
    norm: Ma02MzNorm,
    uplo: Ma02MzUplo,
    n: usize,
    a_re: &[f64],
    a_im: &[f64],
    lda: usize,
    dwork: &mut [f64],
) -> f64 {
    if n == 0 || lda < n {
        return 0.0;
    }
    match norm {
        Ma02MzNorm::MaxAbs => {
            let mut m = 0.0_f64;
            match uplo {
                Ma02MzUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            m = m.max(cabs(a_re[i + j * lda], a_im[i + j * lda]));
                        }
                    }
                }
                Ma02MzUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            m = m.max(cabs(a_re[i + j * lda], a_im[i + j * lda]));
                        }
                    }
                }
            }
            m
        }
        Ma02MzNorm::Frobenius => {
            let mut sum = 0.0_f64;
            match uplo {
                Ma02MzUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            let c = cabs(a_re[i + j * lda], a_im[i + j * lda]);
                            sum += c * c;
                        }
                    }
                }
                Ma02MzUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            let c = cabs(a_re[i + j * lda], a_im[i + j * lda]);
                            sum += c * c;
                        }
                    }
                }
            }
            (2.0 * sum).sqrt()
        }
        Ma02MzNorm::One | Ma02MzNorm::Infinity => {
            if dwork.len() < n {
                return 0.0;
            }
            for i in 0..n {
                dwork[i] = 0.0;
            }
            match uplo {
                Ma02MzUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            let x = cabs(a_re[i + j * lda], a_im[i + j * lda]);
                            dwork[i] += x;
                            dwork[j] += x;
                        }
                    }
                }
                Ma02MzUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            let x = cabs(a_re[i + j * lda], a_im[i + j * lda]);
                            dwork[i] += x;
                            dwork[j] += x;
                        }
                    }
                }
            }
            dwork[..n].iter().fold(0.0_f64, |a, &b| a.max(b))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02mz_zero() {
        let mut dwork = vec![0.0; 1];
        assert_eq!(
            ma02mz(Ma02MzNorm::One, Ma02MzUplo::Upper, 0, &[], &[], 0, &mut dwork),
            0.0
        );
    }

    #[test]
    fn test_ma02mz_frobenius_2x2() {
        // Skew-Hermitian: [0 i; -i 0], upper triangle (0,1) = i => cabs=1
        let a_re = vec![0.0, 0.0, 0.0, 0.0];
        let a_im = vec![0.0, 1.0, -1.0, 0.0];
        let mut dwork = vec![0.0; 2];
        let f = ma02mz(Ma02MzNorm::Frobenius, Ma02MzUplo::Upper, 2, &a_re, &a_im, 2, &mut dwork);
        assert!((f - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_ma02mz_max_abs() {
        let a_re = vec![0.0, 0.0, 0.0, 0.0];
        let a_im = vec![0.0, 3.0, -3.0, 0.0];
        let mut dwork = vec![0.0; 2];
        assert!(
            (ma02mz(Ma02MzNorm::MaxAbs, Ma02MzUplo::Upper, 2, &a_re, &a_im, 2, &mut dwork) - 3.0)
                .abs()
                < 1e-10
        );
    }
}
