//! SB10JD — Positive feedback controller (discrete-time loop shaping).

use nalgebra::DMatrix;

/// Discrete-time loop shaping controller. Placeholder delegates to continuous with same interface.
pub fn sb10jd(
    n: usize,
    m: usize,
    np: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    factor: f64,
    nk: &mut usize,
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
    rcond: &mut [f64],
) -> i32 {
    use crate::sb10::sb10id::sb10id;
    sb10id(n, m, np, a, b, c, d, factor, nk, ak, bk, ck, dk, rcond)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10jd_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 1);
        let c = DMatrix::zeros(1, 0);
        let d = DMatrix::zeros(1, 1);
        let mut nk = 0;
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 2];
        assert_eq!(sb10jd(0, 1, 1, &a, &b, &c, &d, 1.0, &mut nk, &mut ak, &mut bk, &mut ck, &mut dk, &mut rcond), 0);
    }
}
