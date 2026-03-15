//! SB10SD — State feedback and output injection for H-infinity (discrete-time).

use nalgebra::DMatrix;

use crate::sb10::sb10dd::sb10dd;

/// Computes F and H for discrete H-infinity for given gamma (X, Z from SB10DD).
pub fn sb10sd(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    gamma: f64,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    f: &mut DMatrix<f64>,
    h: &mut DMatrix<f64>,
    x: &mut DMatrix<f64>,
    z: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    if n == 0 {
        return 0;
    }
    let mut ak = DMatrix::zeros(n, n);
    let mut bk = DMatrix::zeros(n, nmeas);
    let mut ck = DMatrix::zeros(ncon, n);
    let mut dk = DMatrix::zeros(ncon, nmeas);
    sb10dd(n, m, np, ncon, nmeas, gamma, a, b, c, d, &mut ak, &mut bk, &mut ck, &mut dk, x, z, rcond, tol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10sd_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let d = DMatrix::zeros(2, 2);
        let mut f = DMatrix::zeros(1, 0);
        let mut h = DMatrix::zeros(0, 1);
        let mut x = DMatrix::zeros(0, 0);
        let mut z = DMatrix::zeros(0, 0);
        let mut rcond = [0.0; 8];
        assert_eq!(sb10sd(0, 2, 2, 1, 1, 10.0, &a, &b, &c, &d, &mut f, &mut h, &mut x, &mut z, &mut rcond, 1e-10), 0);
    }
}
