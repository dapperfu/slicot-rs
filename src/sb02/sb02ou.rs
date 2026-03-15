//! SB02OU — SLICOT SB02OU. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb02ou(_n: usize, _a: &DMatrix<f64>, _x: &mut [f64]) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ou() {
        let a = DMatrix::zeros(1, 1);
        let mut x = [0.0];
        assert_eq!(sb02ou(1, &a, &mut x), 0);
    }
}
