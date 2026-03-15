//! SB03SX — SLICOT SB03SX. Stub.
use nalgebra::DMatrix;

pub fn sb03sx(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sb03sx() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb03sx(1, &a, &mut x), 0);
    }
}
