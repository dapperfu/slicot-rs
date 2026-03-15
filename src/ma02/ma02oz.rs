//! MA02OZ — Number of zero rows of a complex (skew-)Hamiltonian matrix (SLICOT MA02OZ)
//
// H = [A D; E ±A']. DE: lower = E, columns 1..M = upper D. Complex; zero = both re and im zero.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02OzSkew {
    Hamiltonian,
    SkewHamiltonian,
}

fn cabs(re: f64, im: f64) -> f64 {
    (re * re + im * im).sqrt()
}

fn get_e_re_im(
    de_re: &[f64],
    de_im: &[f64],
    ldde: usize,
    _m: usize,
    skew: bool,
    i: usize,
    j: usize,
) -> (f64, f64) {
    if i >= j {
        (de_re[i + j * ldde], de_im[i + j * ldde])
    } else if skew {
        (-de_re[j + i * ldde], -de_im[j + i * ldde])
    } else {
        (de_re[j + i * ldde], -de_im[j + i * ldde])
    }
}

fn get_d_re_im(
    de_re: &[f64],
    de_im: &[f64],
    ldde: usize,
    _m: usize,
    skew: bool,
    i: usize,
    j: usize,
) -> (f64, f64) {
    if i <= j {
        (de_re[i + (j + 1) * ldde], de_im[i + (j + 1) * ldde])
    } else if skew {
        (-de_re[j + (i + 1) * ldde], -de_im[j + (i + 1) * ldde])
    } else {
        (de_re[j + (i + 1) * ldde], -de_im[j + (i + 1) * ldde])
    }
}

/// Returns the number of zero rows of H. A and DE are complex, column-major.
pub fn ma02oz(
    skew: Ma02OzSkew,
    m: usize,
    a_re: &[f64],
    a_im: &[f64],
    lda: usize,
    de_re: &[f64],
    de_im: &[f64],
    ldde: usize,
) -> usize {
    if m == 0 || lda < m || ldde < m {
        return 0;
    }
    let skew_bool = skew == Ma02OzSkew::SkewHamiltonian;
    let sign_re = if skew_bool { -1.0 } else { 1.0 };
    let sign_im = 0.0_f64;

    let mut count = 0_usize;

    for i in 0..m {
        let mut row_zero = true;
        for j in 0..m {
            if cabs(a_re[i + j * lda], a_im[i + j * lda]) > 0.0 {
                row_zero = false;
                break;
            }
        }
        if row_zero {
            for j in 0..m {
                let (dr, di) = get_d_re_im(de_re, de_im, ldde, m, skew_bool, i, j);
                if cabs(dr, di) > 0.0 {
                    row_zero = false;
                    break;
                }
            }
        }
        if row_zero {
            count += 1;
        }
    }

    for i in 0..m {
        let mut row_zero = true;
        for j in 0..m {
            let (er, ei) = get_e_re_im(de_re, de_im, ldde, m, skew_bool, i, j);
            if cabs(er, ei) > 0.0 {
                row_zero = false;
                break;
            }
        }
        if row_zero {
            for j in 0..m {
                let ar = sign_re * a_re[j + i * lda] - sign_im * a_im[j + i * lda];
                let ai = sign_re * a_im[j + i * lda] + sign_im * a_re[j + i * lda];
                if cabs(ar, ai) > 0.0 {
                    row_zero = false;
                    break;
                }
            }
        }
        if row_zero {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02oz_zero_dim() {
        assert_eq!(
            ma02oz(Ma02OzSkew::SkewHamiltonian, 0, &[], &[], 0, &[], &[], 0),
            0
        );
    }

    #[test]
    fn test_ma02oz_all_zero() {
        let m = 2;
        let a_re = vec![0.0; 4];
        let a_im = vec![0.0; 4];
        let de_re = vec![0.0; m * (m + 1)];
        let de_im = vec![0.0; m * (m + 1)];
        assert_eq!(
            ma02oz(Ma02OzSkew::SkewHamiltonian, m, &a_re, &a_im, 2, &de_re, &de_im, 2),
            4
        );
    }
}
