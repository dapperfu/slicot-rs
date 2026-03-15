//! SB04OW — SLICOT SB04OW. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb04ow(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb04ow() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb04ow(1, &a, &mut x), 0);
    }
}

