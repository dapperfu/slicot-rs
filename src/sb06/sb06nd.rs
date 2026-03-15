//! SB06ND — SLICOT SB06ND (H2-norm). Stub.

use nalgebra::DMatrix;

/// Stub: returns 0.
pub fn sb06nd(_n: usize, _a: &DMatrix<f64>, _b: &DMatrix<f64>, _c: &DMatrix<f64>, _norm: &mut f64) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb06nd() {
        let a = DMatrix::zeros(1, 1);
        let b = DMatrix::zeros(1, 1);
        let c = DMatrix::zeros(1, 1);
        let mut norm = 0.0;
        assert_eq!(sb06nd(1, &a, &b, &c, &mut norm), 0);
    }
}
