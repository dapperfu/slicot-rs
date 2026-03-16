//! MC03NX — Construction of pencil s*E - A from polynomial matrix (SLICOT MC03NX)
//
// Given P(s) = P(0) + P(1)*s + ... + P(DP)*s^DP, builds the pencil s*E - A as in the spec.

/// P is MP×NP×(DP+1) in increasing powers. A and E are output (DP*MP) × ((DP-1)*MP+NP).
pub fn mc03nx(
    mp: i32,
    np: i32,
    dp: i32,
    p: &[f64],
    ldp1: usize,
    ldp2: usize,
    a: &mut [f64],
    lda: usize,
    e: &mut [f64],
    lde: usize,
) -> i32 {
    let (mp, np, dp) = (mp as usize, np as usize, dp as usize);
    if dp < 1 {
        return -3;
    }
    let nrows = dp * mp;
    let ncols = (dp - 1) * mp + np;
    if lda < nrows.max(1) || lde < nrows.max(1) {
        return -8;
    }
    if a.len() < nrows * lda || e.len() < nrows * lde {
        return -8;
    }
    if ldp1 < mp.max(1) || ldp2 < np.max(1) || p.len() < ldp1 * ldp2 * (dp + 1) {
        return -5;
    }
    let base_p = ldp1 * ldp2;
    for i in 0..nrows {
        for j in 0..ncols {
            a[i + j * lda] = 0.0;
            e[i + j * lde] = 0.0;
        }
    }
    for ib in 0..dp {
        for i in 0..mp {
            let row = ib * mp + i;
            if ib + 1 < dp {
                a[row + ((ib + 1) * mp + i) * lda] = 1.0;
            }
            if ib + 1 < dp {
                e[(ib + 1) * mp + i + (ib * mp + i) * lde] = 1.0;
            }
        }
    }
    for i in 0..mp {
        for j in 0..np {
            a[(dp - 1) * mp + i + ((dp - 1) * mp + j) * lda] =
                p[i + j * ldp1 + 0 * base_p];
        }
    }
    for ib in 0..dp {
        for i in 0..mp {
            for j in 0..np {
                e[ib * mp + i + ((dp - 1) * mp + j) * lde] =
                    -p[i + j * ldp1 + (ib + 1) * base_p];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc03nx_degree1() {
        let mut a = vec![0.0; 4 * 4];
        let mut e = vec![0.0; 4 * 4];
        let p = [1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0];
        assert_eq!(mc03nx(2, 2, 1, &p, 2, 2, &mut a, 2, &mut e, 2), 0);
    }
}
