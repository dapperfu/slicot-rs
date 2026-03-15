//! SB04PD — SLICOT SB04PD. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb04pd(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb04pd() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb04pd(1, &a, &mut x), 0);
    }
}

