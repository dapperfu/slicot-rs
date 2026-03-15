//! SB10YD — H2 optimal controller (discrete-time) full routine.

use nalgebra::DMatrix;

use crate::sb10::sb10ed::sb10ed;

/// H2 optimal discrete-time controller (wrapper around SB10ED).
pub fn sb10yd(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &mut DMatrix<f64>,
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    sb10ed(n, m, np, ncon, nmeas, a, b, c, d, ak, bk, ck, dk, rcond, tol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10yd_n0() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let mut d = DMatrix::zeros(2, 2);
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 8];
        assert_eq!(sb10yd(0, 2, 2, 1, 1, &mut a, &mut b, &c, &mut d, &mut ak, &mut bk, &mut ck, &mut dk, &mut rcond, 1e-10), 0);
    }
}
