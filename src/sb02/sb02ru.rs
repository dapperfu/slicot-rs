//! SB02RU — SLICOT SB02RU. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb02ru(_n: usize, _a: &DMatrix<f64>, _x: &mut [f64]) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ru() {
        let a = DMatrix::zeros(1, 1);
        let mut x = [0.0];
        assert_eq!(sb02ru(1, &a, &mut x), 0);
    }
}
