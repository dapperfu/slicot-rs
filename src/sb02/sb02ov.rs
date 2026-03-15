//! SB02OV — Extract diagonal of N×N matrix into vector (SLICOT support).

use nalgebra::DMatrix;

/// Copies the diagonal of A into b: b(i) = A(i,i). b must have length at least n.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02ov(n: usize, a: &DMatrix<f64>, b: &mut [f64]) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n || b.len() < n {
        return -3;
    }
    for i in 0..n {
        b[i] = a[(i, i)];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ov_diag() {
        let a = DMatrix::from_row_slice(2, 2, &[5.0, 1.0, 1.0, 6.0]);
        let mut b = [0.0; 2];
        assert_eq!(sb02ov(2, &a, &mut b), 0);
        assert!((b[0] - 5.0).abs() < 1e-10);
        assert!((b[1] - 6.0).abs() < 1e-10);
    }
}
