//! SB03MX — SLICOT SB03MX. Stub.
use nalgebra::DMatrix;

pub fn sb03mx(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sb03mx() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb03mx(1, &a, &mut x), 0);
    }
}
