//! SB02MR — Symmetrize N×N matrix from upper triangle (SLICOT support).
//!
//! Sets A(i,j) = A(j,i) for i > j so that A becomes symmetric.

/// Symmetrizes the N×N matrix A in-place from upper triangle (column-major, LDA).
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02mr(n: usize, a: &mut [f64], lda: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if lda < n || a.len() < lda * n {
        return -3;
    }
    for j in 0..n {
        for i in (j + 1)..n {
            a[i + j * lda] = a[j + i * lda];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02mr_n0() {
        let mut a = [0.0];
        assert_eq!(sb02mr(0, &mut a, 1), 0);
    }

    #[test]
    fn test_sb02mr_symmetrize() {
        let mut a = [1.0, 0.5, 2.0, 3.0]; // col-major: (0,0)=1,(1,0)=0.5,(0,1)=2,(1,1)=3. Upper triangle: 1,2,3.
        assert_eq!(sb02mr(2, &mut a, 2), 0);
        assert!((a[1] - 2.0).abs() < 1e-10); // (1,0) := (0,1) = 2
    }
}
