//! SB02MU — Negate N×N matrix in-place (SLICOT support).

/// A := -A. Column-major, leading dimension LDA.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02mu(n: usize, a: &mut [f64], lda: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if lda < n || a.len() < lda * n {
        return -3;
    }
    for k in 0..(lda * n) {
        a[k] = -a[k];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02mu_n0() {
        let mut a = [0.0];
        assert_eq!(sb02mu(0, &mut a, 1), 0);
    }

    #[test]
    fn test_sb02mu_negate() {
        let mut a = [1.0, -2.0, 3.0, 4.0];
        assert_eq!(sb02mu(2, &mut a, 2), 0);
        assert!((a[0] + 1.0).abs() < 1e-10);
        assert!((a[1] - 2.0).abs() < 1e-10);
    }
}
