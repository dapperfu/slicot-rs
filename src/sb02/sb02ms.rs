//! SB02MS — Symmetrize N×N matrix from lower triangle (SLICOT support).
//!
//! Sets A(i,j) = A(j,i) for i < j so that A becomes symmetric.

/// Symmetrizes the N×N matrix A in-place from lower triangle (column-major, LDA).
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02ms(n: usize, a: &mut [f64], lda: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if lda < n || a.len() < lda * n {
        return -3;
    }
    for j in 0..n {
        for i in (j + 1)..n {
            a[j + i * lda] = a[i + j * lda];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ms_n0() {
        let mut a = [0.0];
        assert_eq!(sb02ms(0, &mut a, 1), 0);
    }

    #[test]
    fn test_sb02ms_symmetrize() {
        let mut a = [1.0, 2.0, 0.5, 3.0]; // lower: [1, 2; 0.5, 3] col-major
        assert_eq!(sb02ms(2, &mut a, 2), 0);
        assert!((a[2] - 2.0).abs() < 1e-10); // (0,1) := (1,0)
    }
}
