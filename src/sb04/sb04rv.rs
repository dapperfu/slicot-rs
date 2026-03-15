//! SB04RV — SLICOT SB04RV. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb04rv(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb04rv() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb04rv(1, &a, &mut x), 0);
    }
}

