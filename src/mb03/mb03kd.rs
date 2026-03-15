//! MB03KD — Solve A*X = B (SLICOT MB03). Dense fallback.

use nalgebra::DMatrix;

use crate::mb02::common;

/// Solves A*X = B, overwrites B with X. Returns 0, 1 if singular, -1 if invalid.
pub fn mb03kd(n: usize, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    common::solve_ax_b(n, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb03kd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb03kd(0, &a, &mut b), 0);
    }
}
