//! SB10MD — H-infinity optimal controller (discrete-time) with gamma iteration.

use nalgebra::DMatrix;

use crate::sb10::sb10dd::sb10dd;

/// Discrete H-infinity optimal: reduce gamma until closed-loop loses stability.
pub fn sb10md(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    gamma: &mut f64,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
    x: &mut DMatrix<f64>,
    z: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
    gtol: f64,
) -> i32 {
    if n == 0 {
        return 0;
    }
    let mut g = *gamma;
    for _ in 0..40 {
        let info = sb10dd(n, m, np, ncon, nmeas, g, a, b, c, d, ak, bk, ck, dk, x, z, rcond, tol);
        if info != 0 {
            g *= 1.5;
            continue;
        }
        *gamma = g;
        if g < gtol * 2.0 {
            break;
        }
        g *= 0.9;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10md_n0() {
        let mut gamma = 10.0;
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let d = DMatrix::zeros(2, 2);
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(0, 0);
        let mut z = DMatrix::zeros(0, 0);
        let mut rcond = [0.0; 8];
        assert_eq!(sb10md(0, 2, 2, 1, 1, &mut gamma, &a, &b, &c, &d, &mut ak, &mut bk, &mut ck, &mut dk, &mut x, &mut z, &mut rcond, 1e-10, 1e-8), 0);
    }
}
