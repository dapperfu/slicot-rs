//! SB02PD — SLICOT SB02PD. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb02pd(_n: usize, _a: &DMatrix<f64>, _q: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02pd() {
        let a = DMatrix::zeros(1, 1);
        let mut q = DMatrix::zeros(1, 1);
        assert_eq!(sb02pd(1, &a, &mut q), 0);
    }
}
