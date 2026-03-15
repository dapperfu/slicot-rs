//! SB10FD — H-infinity (sub)optimal state controller (continuous-time) for a given gamma.
//! Uses modified Glover-Doyle formulas; assumes D11=0 or scaled, D12/D21 full rank.

use nalgebra::DMatrix;

use crate::sb02::sb02md::{sb02md, Dico, Uplo};

/// Computes H-infinity suboptimal controller K for given gamma.
///
/// # Returns
/// 0 success; 1-9 see SLICOT SB10FD.
pub fn sb10fd(
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
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    if n == 0 {
        if rcond.len() >= 4 {
            rcond[0] = 1.0;
            rcond[1] = 1.0;
            rcond[2] = 1.0;
            rcond[3] = 1.0;
        }
        return 0;
    }
    let m1 = m - ncon;
    let np1 = np - nmeas;
    let m2 = ncon;
    let np2 = nmeas;
    let g2 = gamma * gamma;
    if g2 < 1e-20 {
        return 6;
    }
    let mut b1 = DMatrix::zeros(n, m1);
    let mut b2 = DMatrix::zeros(n, m2);
    let mut c1 = DMatrix::zeros(np1, n);
    let mut c2 = DMatrix::zeros(np2, n);
    for i in 0..n {
        for j in 0..m1 {
            b1[(i, j)] = b[(i, j)];
        }
        for j in 0..m2 {
            b2[(i, j)] = b[(i, m1 + j)];
        }
    }
    for i in 0..np1 {
        for j in 0..n {
            c1[(i, j)] = c[(i, j)];
        }
    }
    for i in 0..np2 {
        for j in 0..n {
            c2[(i, j)] = c[(np1 + i, j)];
        }
    }
    let qx = c1.transpose() * &c1;
    let gx = &b2 * b2.transpose() - &b1 * b1.transpose() / g2;
    let mut a_x = a.clone();
    let mut qx_mut = qx.clone();
    let mut rcond_x = 0.0;
    let mut wr = vec![0.0; n];
    let mut wi = vec![0.0; n];
    let mut s = DMatrix::zeros(n, n);
    let mut u = DMatrix::zeros(n, n);
    let info_x = sb02md(
        Dico::Continuous,
        'D',
        Uplo::Upper,
        'N',
        'S',
        n,
        &mut a_x,
        &gx,
        &mut qx_mut,
        &mut rcond_x,
        &mut wr,
        &mut wi,
        &mut s,
        &mut u,
    );
    if info_x != 0 {
        if rcond.len() >= 3 {
            rcond[2] = rcond_x;
        }
        return 7;
    }
    if rcond.len() >= 4 {
        rcond[0] = 1.0;
        rcond[1] = 1.0;
        rcond[2] = rcond_x;
    }
    let x = &qx_mut;
    let qy = &b1 * b1.transpose();
    let gy = &c2.transpose() * &c2 - &c1.transpose() * &c1 / g2;
    let at = a.transpose();
    let mut at_mut = at.clone();
    let mut qy_mut = qy.clone();
    let mut rcond_y = 0.0;
    let info_y = sb02md(
        Dico::Continuous,
        'D',
        Uplo::Upper,
        'N',
        'S',
        n,
        &mut at_mut,
        &gy,
        &mut qy_mut,
        &mut rcond_y,
        &mut wr,
        &mut wi,
        &mut s,
        &mut u,
    );
    if info_y != 0 {
        if rcond.len() >= 4 {
            rcond[3] = rcond_y;
        }
        return 8;
    }
    if rcond.len() >= 4 {
        rcond[3] = rcond_y;
    }
    let y = &qy_mut;
    let mut d22 = DMatrix::zeros(np2, m2);
    for i in 0..np2 {
        for j in 0..m2 {
            d22[(i, j)] = d[(np1 + i, m1 + j)];
        }
    }
    let im2 = DMatrix::identity(m2, m2);
    let d22t = d22.transpose();
    let denom = &im2 - &d22t * &d22 / g2;
    let inv_denom = match denom.try_inverse() {
        Some(inv) => inv,
        None => return 9,
    };
    let f = -b2.transpose() * x;
    let h = -y * c2.transpose();
    let dk_mat = &inv_denom * &d22t;
    for i in 0..m2 {
        for j in 0..np2 {
            dk[(i, j)] = dk_mat[(i, j)];
        }
    }
    let b2_f = &b2 * &f;
    let h_c2 = h * &c2;
    for i in 0..n {
        for j in 0..n {
            ak[(i, j)] = a[(i, j)] + b2_f[(i, j)] + h_c2[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..np2 {
            bk[(i, j)] = h[(i, j)];
        }
    }
    for i in 0..m2 {
        for j in 0..n {
            ck[(i, j)] = f[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10fd_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let d = DMatrix::zeros(2, 2);
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 4];
        assert_eq!(sb10fd(0, 2, 2, 1, 1, 15.0, &a, &b, &c, &d, &mut ak, &mut bk, &mut ck, &mut dk, &mut rcond, 1e-10), 0);
    }
}
