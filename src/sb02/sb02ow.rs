//! SB02OW — Copy matrix transpose: X := A' (SLICOT support).

use nalgebra::DMatrix;

/// Copies the transpose of A into X: X = A'. A and X are n×n.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02ow(n: usize, a: &DMatrix<f64>, x: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n || x.nrows() < n || x.ncols() < n {
        return -3;
    }
    for i in 0..n {
        for j in 0..n {
            x[(i, j)] = a[(j, i)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ow_transpose() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let mut x = DMatrix::zeros(2, 2);
        assert_eq!(sb02ow(2, &a, &mut x), 0);
        assert!((x[(0, 1)] - 3.0).abs() < 1e-10); // x = A' => (0,1) = A(1,0) = 3
        assert!((x[(1, 0)] - 2.0).abs() < 1e-10); // (1,0) = A(0,1) = 2
    }
}
