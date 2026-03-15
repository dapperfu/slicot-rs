//! SB03OV — SLICOT SB03OV. Stub.
use nalgebra::DMatrix;

pub fn sb03ov(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sb03ov() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb03ov(1, &a, &mut x), 0);
    }
}
