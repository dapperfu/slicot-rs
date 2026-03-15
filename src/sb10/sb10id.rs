//! SB10ID — Positive feedback controller (McFarlane/Glover loop shaping), continuous-time.

use nalgebra::DMatrix;

use crate::sb02::sb02md::{sb02md, Dico, Uplo};

/// Computes loop-shaping controller for shaped plant G. FACTOR >= 1 (1 = optimal).
pub fn sb10id(
    n: usize,
    m: usize,
    np: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    factor: f64,
    nk: &mut usize,
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
    rcond: &mut [f64],
) -> i32 {
    if n == 0 {
        *nk = 0;
        return 0;
    }
    if factor < 1.0 {
        return -10;
    }
    let qx = c.transpose() * c;
    let gx = &b * b.transpose() / (factor * factor);
    let mut a_x = a.clone();
    let mut qx_mut = qx.clone();
    let mut rcond_x = 0.0;
    let mut wr = vec![0.0; n];
    let mut wi = vec![0.0; n];
    let mut s = DMatrix::zeros(n, n);
    let mut u = DMatrix::zeros(n, n);
    let info_x = sb02md(Dico::Continuous, 'D', Uplo::Upper, 'N', 'S', n, &mut a_x, &gx, &mut qx_mut, &mut rcond_x, &mut wr, &mut wi, &mut s, &mut u);
    if info_x != 0 {
        if rcond.len() >= 1 {
            rcond[0] = rcond_x;
        }
        return 1;
    }
    let x = qx_mut;
    let qy = &b * b.transpose();
    let gy = &c.transpose() * c / (factor * factor);
    let at = a.transpose();
    let mut at_mut = at.clone();
    let mut qy_mut = qy.clone();
    let mut rcond_y = 0.0;
    let info_y = sb02md(Dico::Continuous, 'D', Uplo::Upper, 'N', 'S', n, &mut at_mut, &gy, &mut qy_mut, &mut rcond_y, &mut wr, &mut wi, &mut s, &mut u);
    if info_y != 0 {
        if rcond.len() >= 2 {
            rcond[1] = rcond_y;
        }
        return 2;
    }
    if rcond.len() >= 2 {
        rcond[0] = rcond_x;
        rcond[1] = rcond_y;
    }
    let y = qy_mut;
    let ip = DMatrix::identity(np, np);
    let im = DMatrix::identity(m, m);
    let dd = d * d.transpose() / (factor * factor);
    let mut ip_minus = ip.clone();
    for i in 0..np {
        for j in 0..np {
            ip_minus[(i, j)] -= dd[(i, j)];
        }
    }
    let inv_ip = match ip_minus.try_inverse() {
        Some(inv) => inv,
        None => return 4,
    };
    let dk_val = (inv_ip * d.transpose()).transpose() / (factor * factor);
    let im_minus_dk_d = &im - &dk_val * d;
    if im_minus_dk_d.determinant().abs() < 1e-14 {
        return 5;
    }
    let f = -b.transpose() * &x / (factor * factor);
    let h = -y * c.transpose() / (factor * factor);
    for i in 0..n {
        for j in 0..n {
            ak[(i, j)] = a[(i, j)] + (b * &f)[(i, j)] + (h * c)[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..np {
            bk[(i, j)] = h[(i, j)];
        }
    }
    for i in 0..m {
        for j in 0..n {
            ck[(i, j)] = f[(i, j)];
        }
    }
    for i in 0..m {
        for j in 0..np {
            dk[(i, j)] = dk_val[(i, j)];
        }
    }
    *nk = n;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10id_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 1);
        let c = DMatrix::zeros(1, 0);
        let d = DMatrix::zeros(1, 1);
        let mut nk = 0;
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 2];
        assert_eq!(sb10id(0, 1, 1, &a, &b, &c, &d, 1.0, &mut nk, &mut ak, &mut bk, &mut ck, &mut dk, &mut rcond), 0);
    }
}
