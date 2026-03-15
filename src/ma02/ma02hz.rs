//! MA02HZ — Check if complex A = DIAG*I (identity-like) (SLICOT MA02HZ)

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02HzJob {
    Upper,
    Lower,
    All,
}

/// Returns true if A = DIAG*I (complex), false otherwise. min(M,N)=0 => false.
/// a_re, a_im are column-major LDA×N.
pub fn ma02hz(
    job: Ma02HzJob,
    m: usize,
    n: usize,
    diag_re: f64,
    diag_im: f64,
    a_re: &[f64],
    a_im: &[f64],
    lda: usize,
) -> bool {
    if m == 0 || n == 0 || lda < m {
        return false;
    }
    let eps = 1e-15 * (1.0_f64).max(diag_re.abs()).max(diag_im.abs()).max(1.0);
    match job {
        Ma02HzJob::Upper => {
            for j in 0..n {
                for i in 0..m {
                    if i <= j {
                        let re = a_re[i + j * lda];
                        let im = a_im[i + j * lda];
                        let (exp_re, exp_im) = if i == j {
                            (diag_re, diag_im)
                        } else {
                            (0.0, 0.0)
                        };
                        if (re - exp_re).abs() > eps || (im - exp_im).abs() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        Ma02HzJob::Lower => {
            for j in 0..n {
                for i in 0..m {
                    if i >= j {
                        let re = a_re[i + j * lda];
                        let im = a_im[i + j * lda];
                        let (exp_re, exp_im) = if i == j {
                            (diag_re, diag_im)
                        } else {
                            (0.0, 0.0)
                        };
                        if (re - exp_re).abs() > eps || (im - exp_im).abs() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        Ma02HzJob::All => {
            for j in 0..n {
                for i in 0..m {
                    let re = a_re[i + j * lda];
                    let im = a_im[i + j * lda];
                    let (exp_re, exp_im) = if i == j {
                        (diag_re, diag_im)
                    } else {
                        (0.0, 0.0)
                    };
                    if (re - exp_re).abs() > eps || (im - exp_im).abs() > eps {
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
    fn test_ma02hz_identity() {
        let a_re = vec![1.0, 0.0, 0.0, 1.0];
        let a_im = vec![0.0, 0.0, 0.0, 0.0];
        assert!(ma02hz(
            Ma02HzJob::All,
            2,
            2,
            1.0,
            0.0,
            &a_re,
            &a_im,
            2,
        ));
    }

    #[test]
    fn test_ma02hz_zero_dim() {
        assert!(!ma02hz(Ma02HzJob::All, 0, 0, 1.0, 0.0, &[], &[], 0));
    }
}
