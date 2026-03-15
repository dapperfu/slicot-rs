//! MA02GZ — Column interchanges on a complex matrix (SLICOT MA02GZ)
//
// For each column index j in [k1, k2], swap column j with column ipiv(j).
// ipiv accessed with stride incx (if incx < 0, apply in reverse order k2 down to k1).

/// Applies column interchanges. Indices 0-based. a_re, a_im column-major LDA×M (M >= max column index).
pub fn ma02gz(
    n: usize,
    a_re: &mut [f64],
    a_im: &mut [f64],
    lda: usize,
    k1: usize,
    k2: usize,
    ipiv: &[i32],
    incx: i32,
) -> i32 {
    if incx == 0 || n == 0 {
        return 0;
    }
    if lda < n {
        return -3;
    }
    let m = a_re.len().min(a_im.len()) / lda.max(1);
    if m == 0 {
        return 0;
    }

    let k2_eff = k2.min(m.saturating_sub(1));
    let mut jx: i32 = if incx > 0 {
        k1 as i32
    } else {
        k2_eff as i32
    };

    let (start, end, step): (i32, i32, i32) = if incx > 0 {
        (k1 as i32, k2_eff as i32, 1)
    } else {
        (k2_eff as i32, k1 as i32, -1)
    };

    let mut j = start;
    while (step > 0 && j <= end) || (step < 0 && j >= end) {
        let ju = j as usize;
        let jp = ipiv
            .get(jx.max(0) as usize)
            .copied()
            .unwrap_or(ju as i32) as usize;
        if jp != ju && jp < m {
            for i in 0..n {
                let idx_j = i + ju * lda;
                let idx_jp = i + jp * lda;
                a_re.swap(idx_j, idx_jp);
                a_im.swap(idx_j, idx_jp);
            }
        }
        jx += incx;
        j += step;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02gz_swap_cols() {
        // 2×3 column-major: col0=(1,2), col1=(3,4), col2=(5,6). Swap col1 with col2 (k1=k2=1 only).
        let mut a_re = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut a_im = [0.0; 6];
        let ipiv = [0, 2, 1]; // for j=1: swap col 1 with col 2
        ma02gz(2, &mut a_re, &mut a_im, 2, 1, 1, &ipiv, 1);
        assert_eq!(a_re[2], 5.0);
        assert_eq!(a_re[4], 3.0);
    }
}
