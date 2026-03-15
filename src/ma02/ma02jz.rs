//! MA02JZ — Residual || Q^H Q - I ||_F for complex unitary symplectic Q (SLICOT MA02JZ)
//
// Q = [ op(Q1) op(Q2); -op(Q2) op(Q1) ]. res is workspace (unused; kept for API).

/// Computes residual for complex Q. q1_re, q1_im, q2_re, q2_im are n×n column-major (LDQ×n).
/// (Q^H*Q)(1,1) = Q1^H*Q1 + Q2^H*Q2, (Q^H*Q)(1,2) = Q1^H*Q2 - Q2^H*Q1.
pub fn ma02jz(
    tran1: bool,
    tran2: bool,
    n: usize,
    q1_re: &[f64],
    q1_im: &[f64],
    ldq1: usize,
    q2_re: &[f64],
    q2_im: &[f64],
    ldq2: usize,
    _res: &mut [f64],
) -> f64 {
    if n == 0 || ldq1 < n || ldq2 < n {
        return 0.0;
    }
    // P11 = Q1^H*Q1 + Q2^H*Q2 - I (complex), P12 = Q1^H*Q2 - Q2^H*Q1 (complex)
    let mut sum_sq = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let mut p11_re = -if i == j { 1.0 } else { 0.0 };
            let mut p11_im = 0.0;
            let mut p12_re = 0.0;
            let mut p12_im = 0.0;
            for k in 0..n {
                let (q1_ki_re, q1_ki_im) = if tran1 {
                    (q1_re[i + k * ldq1], -q1_im[i + k * ldq1])
                } else {
                    (q1_re[k + i * ldq1], -q1_im[k + i * ldq1])
                };
                let (q1_kj_re, q1_kj_im) = if tran1 {
                    (q1_re[j + k * ldq1], q1_im[j + k * ldq1])
                } else {
                    (q1_re[k + j * ldq1], q1_im[k + j * ldq1])
                };
                let (q2_ki_re, q2_ki_im) = if tran2 {
                    (q2_re[i + k * ldq2], -q2_im[i + k * ldq2])
                } else {
                    (q2_re[k + i * ldq2], -q2_im[k + i * ldq2])
                };
                let (q2_kj_re, q2_kj_im) = if tran2 {
                    (q2_re[j + k * ldq2], q2_im[j + k * ldq2])
                } else {
                    (q2_re[k + j * ldq2], q2_im[k + j * ldq2])
                };
                // conj(q1_ki)*q1_kj
                p11_re += q1_ki_re * q1_kj_re + q1_ki_im * q1_kj_im;
                p11_im += q1_ki_re * q1_kj_im - q1_ki_im * q1_kj_re;
                // conj(q2_ki)*q2_kj
                p11_re += q2_ki_re * q2_kj_re + q2_ki_im * q2_kj_im;
                p11_im += q2_ki_re * q2_kj_im - q2_ki_im * q2_kj_re;
                // conj(q1_ki)*q2_kj - conj(q2_ki)*q1_kj
                p12_re += q1_ki_re * q2_kj_re + q1_ki_im * q2_kj_im - (q2_ki_re * q1_kj_re + q2_ki_im * q1_kj_im);
                p12_im += q1_ki_re * q2_kj_im - q1_ki_im * q2_kj_re - (q2_ki_re * q1_kj_im - q2_ki_im * q1_kj_re);
            }
            sum_sq += p11_re * p11_re + p11_im * p11_im;
            sum_sq += p12_re * p12_re + p12_im * p12_im;
        }
    }
    // Block (2,1) = -conj(P12), block (2,2) = P11; count each once (P11 twice total, P12 twice)
    sum_sq.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02jz_identity() {
        let n = 2;
        let q1_re = vec![1.0, 0.0, 0.0, 1.0];
        let q1_im = vec![0.0, 0.0, 0.0, 0.0];
        let q2_re = vec![0.0, 0.0, 0.0, 0.0];
        let q2_im = vec![0.0, 0.0, 0.0, 0.0];
        let mut res = vec![0.0; 4];
        let r = ma02jz(false, false, n, &q1_re, &q1_im, 2, &q2_re, &q2_im, 2, &mut res);
        assert!(r < 1e-10, "identity Q should give residual ~0, got {}", r);
    }

    #[test]
    fn test_ma02jz_zero_dim() {
        let mut res = vec![0.0];
        assert_eq!(ma02jz(false, false, 0, &[], &[], 0, &[], &[], 0, &mut res), 0.0);
    }
}
