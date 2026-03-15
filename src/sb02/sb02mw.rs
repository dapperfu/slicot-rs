//! SB02MW — Set N×N matrix to identity (SLICOT support).

/// A := I (identity). Column-major, LDA.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02mw(n: usize, a: &mut [f64], lda: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if lda < n || a.len() < lda * n {
        return -3;
    }
    for j in 0..n {
        for i in 0..n {
            a[i + j * lda] = if i == j { 1.0 } else { 0.0 };
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02mw_n0() {
        let mut a = [0.0];
        assert_eq!(sb02mw(0, &mut a, 1), 0);
    }

    #[test]
    fn test_sb02mw_identity() {
        let mut a = [0.0; 4];
        assert_eq!(sb02mw(2, &mut a, 2), 0);
        assert!((a[0] - 1.0).abs() < 1e-10);
        assert!((a[1] - 0.0).abs() < 1e-10);
        assert!((a[2] - 0.0).abs() < 1e-10);
        assert!((a[3] - 1.0).abs() < 1e-10);
    }
}
