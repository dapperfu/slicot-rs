//! SB10ED — H2 optimal state controller (discrete-time).

use nalgebra::DMatrix;

use crate::sb10::dare;
use crate::sb10::sb10ud::sb10ud;

/// H2 optimal discrete-time controller. Uses SB10UD-style normalization then two DAREs.
pub fn sb10ed(
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
    if n == 0 {
        return 0;
    }
    let m1 = m - ncon;
    let np1 = np - nmeas;
    let m2 = ncon;
    let np2 = nmeas;
    let mut b_work = b.clone();
    let mut c_work = c.clone();
    let mut d_work = d.clone();
    let mut tu = DMatrix::zeros(m2, m2);
    let mut ty = DMatrix::zeros(np2, np2);
    let mut rcond_ud = [0.0; 2];
    let info_ud = sb10ud(n, m, np, ncon, nmeas, &mut b_work, &mut c_work, &mut d_work, &mut tu, &mut ty, &mut rcond_ud, tol);
    if info_ud != 0 {
        return info_ud + 2;
    }
    if rcond.len() >= 7 {
        rcond[0] = rcond_ud[0];
        rcond[1] = rcond_ud[1];
    }
    let mut c1 = DMatrix::zeros(np1, n);
    let mut b2 = DMatrix::zeros(n, m2);
    for i in 0..np1 {
        for j in 0..n {
            c1[(i, j)] = c_work[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..m2 {
            b2[(i, j)] = b_work[(i, m1 + j)];
        }
    }
    let qx = c1.transpose() * &c1;
    let rx = DMatrix::identity(m2, m2);
    let mut x = DMatrix::zeros(n, n);
    if !dare::dare(a, &b2, &qx, &rx, &mut x, tol, 80) {
        return 6;
    }
    let mut c2 = DMatrix::zeros(np2, n);
    let mut b1 = DMatrix::zeros(n, m1);
    for i in 0..np2 {
        for j in 0..n {
            c2[(i, j)] = c_work[(np1 + i, j)];
        }
    }
    for i in 0..n {
        for j in 0..m1 {
            b1[(i, j)] = b_work[(i, j)];
        }
    }
    let qz = &b1 * b1.transpose();
    let rz = DMatrix::identity(np2, np2);
    let at = a.transpose();
    let c2t = c2.transpose();
    let mut z = DMatrix::zeros(n, n);
    if !dare::dare(&at, &c2t, &qz, &rz, &mut z, tol, 80) {
        return 8;
    }
    let rpbxb = &rx + b2.transpose() * &x * &b2;
    let rpbxb_inv = match rpbxb.try_inverse() {
        Some(inv) => inv,
        None => return 7,
    };
    let f = -&rpbxb_inv * b2.transpose() * &x * &*a;
    let rpcyc = &rz + &c2 * &z * &c2t;
    let rpcyc_inv = match rpcyc.try_inverse() {
        Some(inv) => inv,
        None => return 9,
    };
    let l = -&*a * &z * c2.transpose() * &rpcyc_inv;
    for i in 0..n {
        for j in 0..n {
            ak[(i, j)] = a[(i, j)] + (b2.clone() * &f)[(i, j)] + (l.clone() * &c2)[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..np2 {
            bk[(i, j)] = l[(i, j)];
        }
    }
    for i in 0..m2 {
        for j in 0..n {
            ck[(i, j)] = f[(i, j)];
        }
    }
    for i in 0..m2 {
        for j in 0..np2 {
            dk[(i, j)] = 0.0;
        }
    }
    b.copy_from(&b_work);
    d.copy_from(&d_work);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10ed_n0() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let mut d = DMatrix::zeros(2, 2);
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 8];
        assert_eq!(sb10ed(0, 2, 2, 1, 1, &mut a, &mut b, &c, &mut d, &mut ak, &mut bk, &mut ck, &mut dk, &mut rcond, 1e-10), 0);
    }
}
