//! SB10ZD — Normalization for H-infinity (discrete-time).

use nalgebra::DMatrix;

use crate::sb10::sb10ud::sb10ud;

/// Normalizes discrete-time plant for H-infinity design.
pub fn sb10zd(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    tu: &mut DMatrix<f64>,
    ty: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    sb10ud(n, m, np, ncon, nmeas, b, c, d, tu, ty, rcond, tol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10zd_n0() {
        let mut b = DMatrix::zeros(0, 2);
        let mut c = DMatrix::zeros(2, 0);
        let mut d = DMatrix::zeros(2, 2);
        let mut tu = DMatrix::zeros(1, 1);
        let mut ty = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 2];
        assert_eq!(sb10zd(0, 2, 2, 1, 1, &mut b, &mut c, &mut d, &mut tu, &mut ty, &mut rcond, 1e-10), 0);
    }
}
