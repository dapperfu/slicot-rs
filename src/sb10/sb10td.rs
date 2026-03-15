//! SB10TD — H-infinity controller from F, H (discrete-time).

use nalgebra::DMatrix;

use crate::sb10::sb10wd::sb10wd;

/// Forms discrete H-infinity controller K from F, H, TU, TY.
pub fn sb10td(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    f: &DMatrix<f64>,
    h: &DMatrix<f64>,
    tu: &DMatrix<f64>,
    ty: &DMatrix<f64>,
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
) -> i32 {
    sb10wd(n, m, np, ncon, nmeas, a, b, c, d, f, h, tu, ty, ak, bk, ck, dk)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10td_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let d = DMatrix::zeros(2, 2);
        let f = DMatrix::zeros(1, 0);
        let h = DMatrix::zeros(0, 1);
        let tu = DMatrix::zeros(1, 1);
        let ty = DMatrix::zeros(1, 1);
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        assert_eq!(sb10td(0, 2, 2, 1, 1, &a, &b, &c, &d, &f, &h, &tu, &ty, &mut ak, &mut bk, &mut ck, &mut dk), 0);
    }
}
