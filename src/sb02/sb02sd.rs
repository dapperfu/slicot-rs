//! SB02SD — Copy matrix: X := A (SLICOT support).

use nalgebra::DMatrix;

/// Copies A into X: X = A. A and X are n×n.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02sd(n: usize, a: &DMatrix<f64>, x: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n || x.nrows() < n || x.ncols() < n {
        return -3;
    }
    for i in 0..n {
        for j in 0..n {
            x[(i, j)] = a[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02sd() {
        let a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb02sd(1, &a, &mut x), 0);
        assert!((x[(0, 0)] - 1.0).abs() < 1e-10);
    }
}
