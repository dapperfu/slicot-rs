//! SB02OX — SLICOT SB02OX. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb02ox(_n: usize, _a: &DMatrix<f64>, _b: &mut [f64]) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ox() {
        let a = DMatrix::zeros(1, 1);
        let mut b = [0.0];
        assert_eq!(sb02ox(1, &a, &mut b), 0);
    }
}
