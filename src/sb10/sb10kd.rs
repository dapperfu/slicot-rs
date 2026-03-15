//! SB10KD — H2 optimal controller (continuous-time) full routine.

use nalgebra::DMatrix;

use crate::sb10::sb10ud::sb10ud;
use crate::sb10::sb10vd::sb10vd;
use crate::sb10::sb10wd::sb10wd;

/// Computes H2 optimal controller: SB10UD (normalize) then SB10VD (F,H,X,Y) then SB10WD (K).
pub fn sb10kd(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    a: &DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    if n == 0 {
        return 0;
    }
    let mut tu = DMatrix::zeros(ncon, ncon);
    let mut ty = DMatrix::zeros(nmeas, nmeas);
    let mut rcond_ud = [0.0; 2];
    let info_ud = sb10ud(n, m, np, ncon, nmeas, b, c, d, &mut tu, &mut ty, &mut rcond_ud, tol);
    if info_ud != 0 {
        return info_ud;
    }
    let mut f = DMatrix::zeros(ncon, n);
    let mut h = DMatrix::zeros(n, nmeas);
    let mut x = DMatrix::zeros(n, n);
    let mut y = DMatrix::zeros(n, n);
    let mut xycond = [0.0; 2];
    let info_vd = sb10vd(n, m, np, ncon, nmeas, a, b, c, &mut f, &mut h, &mut x, &mut y, &mut xycond);
    if info_vd != 0 {
        return info_vd + 10;
    }
    let info_wd = sb10wd(n, m, np, ncon, nmeas, a, b, c, d, &f, &h, &tu, &ty, ak, bk, ck, dk);
    if info_wd != 0 {
        return info_wd + 20;
    }
    if rcond.len() >= 2 {
        rcond[0] = xycond[0];
        rcond[1] = xycond[1];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10kd_n0() {
        let a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 2);
        let mut c = DMatrix::zeros(2, 0);
        let mut d = DMatrix::zeros(2, 2);
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 7];
        assert_eq!(sb10kd(0, 2, 2, 1, 1, &a, &mut b, &mut c, &mut d, &mut ak, &mut bk, &mut ck, &mut dk, &mut rcond, 1e-10), 0);
    }
}
