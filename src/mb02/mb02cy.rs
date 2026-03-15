//! MB02CY — Solve A*X = B (SLICOT MB02). Dense fallback; B n×k.

use nalgebra::DMatrix;

use crate::mb02::common;

/// Solves A*X = B (A n×n, B n×k), overwrites B with X. Returns 0, 1 if singular, -1 if invalid.
pub fn mb02cy(n: usize, k: usize, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    if k == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || b.nrows() != n || b.ncols() != k {
        return -1;
    }
    common::solve_ax_b(n, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb02cy_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb02cy(0, 0, &a, &mut b), 0);
    }
}
