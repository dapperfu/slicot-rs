//! SB10PD — Normalization for H-infinity controller design (continuous-time).

use nalgebra::DMatrix;

/// Normalizes plant for H-infinity design (D12, D21 to unit form with gamma scaling).
pub fn sb10pd(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    gamma: f64,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    if gamma <= 0.0 {
        return -7;
    }
    use crate::sb10::sb10ud::sb10ud;
    let mut tu = DMatrix::zeros(ncon, ncon);
    let mut ty = DMatrix::zeros(nmeas, nmeas);
    sb10ud(n, m, np, ncon, nmeas, b, c, d, &mut tu, &mut ty, rcond, tol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10pd_n0() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 2);
        let mut c = DMatrix::zeros(2, 0);
        let mut d = DMatrix::zeros(2, 2);
        let mut rcond = [0.0; 2];
        assert_eq!(sb10pd(0, 2, 2, 1, 1, 10.0, &mut a, &mut b, &mut c, &mut d, &mut rcond, 1e-10), 0);
    }
}
