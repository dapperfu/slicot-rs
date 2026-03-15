//! MA02CZ — Pertranspose central band of a complex square matrix (SLICOT MA02CZ)
//
// Same as MA02CD but for complex A. a_re, a_im column-major LDA×N.

/// Pertransposes the central band (KL subdiagonals, diagonal, KU superdiagonals) in place.
/// a_re and a_im are column-major, length at least lda*n.
pub fn ma02cz(
    n: usize,
    kl: usize,
    ku: usize,
    a_re: &mut [f64],
    a_im: &mut [f64],
    lda: usize,
) -> i32 {
    if n <= 1 {
        return 0;
    }
    if lda < n {
        return -5;
    }
    let need = lda * n;
    if a_re.len() < need || a_im.len() < need {
        return -4;
    }

    for i in 1..=kl.min(n.saturating_sub(2)) {
        let i1 = (n - i) / 2;
        if i1 > 0 {
            for j in 0..i1 {
                let r1 = i + j;
                let c1 = j;
                let r2 = n - i1 + j;
                let c2 = n - i1 - i + j;
                if r2 < n && c2 < n {
                    let idx1 = r1 + c1 * lda;
                    let idx2 = r2 + c2 * lda;
                    a_re.swap(idx1, idx2);
                    a_im.swap(idx1, idx2);
                }
            }
        }
    }
    for i in 1..=ku.min(n.saturating_sub(2)) {
        let i1 = (n - i) / 2;
        if i1 > 0 {
            for j in 0..i1 {
                let r1 = j;
                let c1 = i + j;
                let r2 = n - i1 - i + j;
                let c2 = n - i1 + j;
                if r2 < n && c2 < n {
                    let idx1 = r1 + c1 * lda;
                    let idx2 = r2 + c2 * lda;
                    a_re.swap(idx1, idx2);
                    a_im.swap(idx1, idx2);
                }
            }
        }
    }
    let i1 = n / 2;
    if i1 > 0 {
        for j in 0..i1 {
            let r1 = j;
            let c1 = j;
            let r2 = n - i1 + j;
            let c2 = n - i1 + j;
            if r2 < n && c2 < n {
                let idx1 = r1 + c1 * lda;
                let idx2 = r2 + c2 * lda;
                a_re.swap(idx1, idx2);
                a_im.swap(idx1, idx2);
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02cz_diag() {
        let mut a_re = [1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0];
        let mut a_im = [0.0; 9];
        assert_eq!(ma02cz(3, 0, 0, &mut a_re, &mut a_im, 3), 0);
        assert_eq!(a_re[0], 3.0);
        assert_eq!(a_re[4], 2.0);
        assert_eq!(a_re[8], 1.0);
    }
}
