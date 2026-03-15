//! SB02MX — Matrix-vector product b := A*b (SLICOT support).

use nalgebra::DMatrix;

/// Computes b := A * b (in-place: b overwritten with result).
/// A is n×n, b has length n.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02mx(n: usize, a: &DMatrix<f64>, b: &mut [f64]) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || b.len() < n {
        return -3;
    }
    let b_vec: Vec<f64> = (0..n).map(|i| b[i]).collect();
    for i in 0..n {
        b[i] = (0..n).map(|j| a[(i, j)] * b_vec[j]).sum();
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02mx_n0() {
        let a = DMatrix::zeros(0, 0);
        let mut b = [0.0];
        assert_eq!(sb02mx(0, &a, &mut b), 0);
    }

    #[test]
    fn test_sb02mx_mv() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let mut b = [1.0, 0.0];
        assert_eq!(sb02mx(2, &a, &mut b), 0);
        assert!((b[0] - 1.0).abs() < 1e-10);
        assert!((b[1] - 3.0).abs() < 1e-10);
    }
}
