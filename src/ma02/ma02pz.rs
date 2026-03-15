//! MA02PZ — Count zero rows and zero columns of a complex matrix (SLICOT MA02PZ)
//
// A is M×N, column-major: a_re[i + j*lda], a_im[i + j*lda].

/// Counts zero rows and zero columns. a_re and a_im are column-major, length at least lda*n.
/// Returns 0; nzr and nzc are set.
pub fn ma02pz(
    m: usize,
    n: usize,
    a_re: &[f64],
    a_im: &[f64],
    lda: usize,
    nzr: &mut i32,
    nzc: &mut i32,
) -> i32 {
    *nzr = 0;
    *nzc = 0;
    if m == 0 || n == 0 {
        return 0;
    }
    if lda < m {
        return -5;
    }
    let need = lda * n;
    if a_re.len() < need || a_im.len() < need {
        return -4;
    }

    for j in 0..n {
        let mut zero = true;
        for i in 0..m {
            let idx = i + j * lda;
            if a_re[idx] != 0.0 || a_im[idx] != 0.0 {
                zero = false;
                break;
            }
        }
        if zero {
            *nzc += 1;
        }
    }
    for i in 0..m {
        let mut zero = true;
        for j in 0..n {
            let idx = i + j * lda;
            if a_re[idx] != 0.0 || a_im[idx] != 0.0 {
                zero = false;
                break;
            }
        }
        if zero {
            *nzr += 1;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02pz() {
        let a_re = [0.0, 0.0, 1.0, 0.0];
        let a_im = [0.0, 0.0, 0.0, 0.0];
        let mut nzr = -1;
        let mut nzc = -1;
        assert_eq!(ma02pz(2, 2, &a_re, &a_im, 2, &mut nzr, &mut nzc), 0);
        assert_eq!(nzr, 1);
        assert_eq!(nzc, 1);
    }
}
