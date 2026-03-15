//! SB08CD — SLICOT SB08CD. Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb08cd(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sb08cd() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb08cd(1, &a, &mut x), 0);
    }
}
