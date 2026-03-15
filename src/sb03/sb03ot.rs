//! SB03OT — SLICOT SB03OT. Stub.
use nalgebra::DMatrix;

pub fn sb03ot(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sb03ot() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb03ot(1, &a, &mut x), 0);
    }
}
