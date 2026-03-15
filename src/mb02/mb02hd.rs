//! MB02HD — Solve A*X = B (SLICOT MB02). Dense fallback.

use nalgebra::DMatrix;

use crate::mb02::common;

/// Solves A*X = B, overwrites B with X. Returns 0, 1 if singular, -1 if invalid.
pub fn mb02hd(n: usize, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    common::solve_ax_b(n, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb02hd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb02hd(0, &a, &mut b), 0);
    }

    #[test]
    fn test_mb02hd_solve() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 0.0, 2.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 2.0]);
        assert_eq!(mb02hd(2, &a, &mut b), 0);
        assert!((b[(0, 0)] - 0.0).abs() < 1e-10);
        assert!((b[(1, 0)] - 1.0).abs() < 1e-10);
    }
}
