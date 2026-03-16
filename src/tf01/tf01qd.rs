//! TF01QD — Markov parameters from transfer function matrix G(z) (SLICOT TF01QD)
//!
//! Computes M(1),...,M(N) from AR/MA coefficients per (i,j) element.

/// IORD(i,j) = order r of G_ij. AR and MA stored row-by-row (1,1),(1,2),...,(NC,NB), each with r coefficients (decreasing powers).
/// H(i, (k-1)*NB+j) = (i,j) element of M(k).
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tf01qd(
    nc: usize,
    nb: usize,
    n: usize,
    iord: &[i32],
    ar: &[f64],
    ma: &[f64],
    h: &mut [f64],
    ldh: usize,
) -> i32 {
    if nc == 0 || nb == 0 || n == 0 {
        return 0;
    }
    if iord.len() < nc * nb {
        return -5;
    }
    if h.len() < nc * (n * nb) || ldh < nc {
        return -8;
    }

    let mut ar_offset = 0_usize;
    let mut ma_offset = 0_usize;
    for i in 0..nc {
        for j in 0..nb {
            let r = iord[i * nb + j] as usize;
            if r == 0 {
                continue;
            }
            if ar_offset + r > ar.len() || ma_offset + r > ma.len() {
                return -6;
            }
            for k in 1..=n {
                let idx = i + ((k - 1) * nb + j) * ldh;
                if idx >= h.len() {
                    continue;
                }
                let (ar_ij, ma_ij) = (&ar[ar_offset..ar_offset + r], &ma[ma_offset..ma_offset + r]);
                let mut m_k = if k <= r { ma_ij[k - 1] } else { 0.0 };
                if k <= r {
                    for p in 1..k {
                        m_k -= ar_ij[p - 1] * h[i + ((k - p - 1) * nb + j) * ldh];
                    }
                } else {
                    for p in 0..r {
                        let col = (k - p - 2) * nb + j;
                        m_k -= ar_ij[p] * h[i + col * ldh];
                    }
                }
                h[idx] = m_k;
            }
            ar_offset += r;
            ma_offset += r;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tf01qd_smoke() {
        let iord = [1, 1, 1, 1];
        let ar = [0.5, -0.8, 0.5, -0.8];
        let ma = [1.0, 1.0, 1.0, 1.0];
        let mut h = vec![0.0; 2 * 3 * 2];
        assert_eq!(tf01qd(2, 2, 3, &iord, &ar, &ma, &mut h, 2), 0);
        assert!(h[0].is_finite());
    }
}
