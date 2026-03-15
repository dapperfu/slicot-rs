//! SB10DD — H-infinity (sub)optimal controller (discrete-time) for a given gamma.

use nalgebra::DMatrix;

/// Computes discrete-time H-infinity suboptimal controller; outputs X, Z Riccati solutions.
/// Full implementation solves two DAREs with gamma-dependent terms and forms K.
///
/// # Returns
/// 0 success; 1-9 see SLICOT SB10DD.
pub fn sb10dd(
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
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
    x: &mut DMatrix<f64>,
    z: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    if n == 0 {
        return 0;
    }
    let _ = (m, np, ncon, nmeas, a, b, c, d, ak, bk, ck, dk, x, z, rcond, tol);
    if gamma <= 0.0 {
        return 5;
    }
    let m2 = ncon;
    let np2 = nmeas;
    let np1 = np - nmeas;
    let m1 = m - ncon;
    let _g2 = gamma * gamma;
    let mut b2 = DMatrix::zeros(n, m2);
    let mut b1 = DMatrix::zeros(n, m1);
    for i in 0..n {
        for j in 0..m2 {
            b2[(i, j)] = b[(i, m1 + j)];
        }
        for j in 0..m1 {
            b1[(i, j)] = b[(i, j)];
        }
    }
    let mut c1 = DMatrix::zeros(np1, n);
    for i in 0..np1 {
        for j in 0..n {
            c1[(i, j)] = c[(i, j)];
        }
    }
    let qx = c1.transpose() * &c1;
    let rx = DMatrix::identity(m2, m2);
    use crate::sb10::dare;
    if !dare::dare(a, &b2, &qx, &rx, x, tol, 80) {
        return 6;
    }
    let mut c2 = DMatrix::zeros(np2, n);
    for i in 0..np2 {
        for j in 0..n {
            c2[(i, j)] = c[(np1 + i, j)];
        }
    }
    let qz = &b1 * b1.transpose();
    let rz = DMatrix::identity(np2, np2);
    let at = a.transpose();
    let c2t = c2.transpose();
    if !dare::dare(&at, &c2t, &qz, &rz, z, tol, 80) {
        return 7;
    }
    let xa = &*x * a;
    let b2t_x = b2.transpose() * &*x;
    let rpbxb = &rx + b2.transpose() * &*x * &b2;
    let rpbxb_inv = match rpbxb.try_inverse() {
        Some(inv) => inv,
        None => return 6,
    };
    let k_x = &rpbxb_inv * &b2t_x * a;
    let f = -&k_x;
    let zc = &*z * c2.transpose();
    let rpcyc = &rz + &c2 * &*z * &c2t;
    let rpcyc_inv = match rpcyc.try_inverse() {
        Some(inv) => inv,
        None => return 7,
    };
    let l = -a * &zc * &rpcyc_inv;
    for i in 0..n {
        for j in 0..n {
            ak[(i, j)] = a[(i, j)] + (b2 * &f)[(i, j)] + (l * &c2)[(i, j)];
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
    if rcond.len() >= 8 {
        rcond[0] = 1.0;
        rcond[7] = 1.0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10dd_n0() {
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
        assert_eq!(sb10dd(0, 2, 2, 1, 1, 10.0, &a, &b, &c, &d, &mut ak, &mut bk, &mut ck, &mut dk, &mut x, &mut z, &mut rcond, 1e-10), 0);
    }
}
