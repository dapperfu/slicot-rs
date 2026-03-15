//! SB02MV — Normalize N×N matrix so max |A(i,j)| = 1 (SLICOT support).

/// A := A / max_ij |A(i,j)|. Column-major, LDA. No-op if max is zero.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02mv(n: usize, a: &mut [f64], lda: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if lda < n || a.len() < lda * n {
        return -3;
    }
    let mut scale = 0.0_f64;
    for j in 0..n {
        for i in 0..n {
            scale = scale.max(a[i + j * lda].abs());
        }
    }
    if scale > 0.0 {
        for k in 0..(lda * n) {
            a[k] /= scale;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02mv_n0() {
        let mut a = [0.0];
        assert_eq!(sb02mv(0, &mut a, 1), 0);
    }

    #[test]
    fn test_sb02mv_normalize() {
        let mut a = [2.0, 4.0, 6.0, 8.0];
        assert_eq!(sb02mv(2, &mut a, 2), 0);
        assert!((a[0] - 0.25).abs() < 1e-10);
        assert!((a[3] - 1.0).abs() < 1e-10);
    }
}
