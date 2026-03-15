//! SB04RD — SLICOT SB04RD. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb04rd(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb04rd() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb04rd(1, &a, &mut x), 0);
    }
}

