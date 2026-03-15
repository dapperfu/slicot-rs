//! SB10VD — State feedback and output injection matrices for H2 optimal controller (continuous-time).
//! Assumes D12 = [0;I] and D21 = [0 I] (from SB10UD). Solves X- and Y-Riccati equations and forms F, H.

use nalgebra::DMatrix;

use crate::sb02::sb02md::{sb02md, Dico, Uplo};

/// Computes state feedback F and output injection H for H2 optimal controller.
/// X and Y are solutions of the two CARE; XYCOND(1)=rcond(X), XYCOND(2)=rcond(Y).
///
/// # Returns
/// 0 success; 1 X-Riccati failed; 2 Y-Riccati failed; < 0 invalid argument.
pub fn sb10vd(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    f: &mut DMatrix<f64>,
    h: &mut DMatrix<f64>,
    x: &mut DMatrix<f64>,
    y: &mut DMatrix<f64>,
    xycond: &mut [f64],
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -7;
    }
    if b.nrows() != n || b.ncols() != m {
        return -9;
    }
    if c.nrows() != np || c.ncols() != n {
        return -11;
    }
    let m1 = m - ncon;
    let np1 = np - nmeas;
    let m2 = ncon;
    let np2 = nmeas;
    if f.nrows() != m2 || f.ncols() != n {
        return -13;
    }
    if h.nrows() != n || h.ncols() != np2 {
        return -15;
    }
    if x.nrows() != n || x.ncols() != n {
        return -17;
    }
    if y.nrows() != n || y.ncols() != n {
        return -19;
    }
    if xycond.len() < 2 {
        return -21;
    }
    xycond[0] = 0.0;
    xycond[1] = 0.0;

    if n == 0 {
        return 0;
    }

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
    let mut c2 = DMatrix::zeros(np2, n);
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

    // X-Riccati: A'*X + X*A - X*B2*B2'*X + C1'*C1 = 0
    let gx = &b2 * b2.transpose();
    let qx = c1.transpose() * &c1;
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
        xycond[0] = rcond_x;
        return 1;
    }
    xycond[0] = rcond_x;
    x.copy_from(&qx_mut);

    // Y-Riccati: A*Y + Y*A' - Y*C2'*C2*Y + B1*B1' = 0. Dual: (A')'*Z + Z*A' - Z*C2'*C2*Z + B1*B1' = 0 with Z=Y.
    let at = a.transpose();
    let gy = &c2.transpose() * &c2;
    let qy = &b1 * b1.transpose();
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
        xycond[1] = rcond_y;
        return 2;
    }
    xycond[1] = rcond_y;
    y.copy_from(&qy_mut);

    // F = -B2'*X (M2×N), i.e. F(i,j) = - sum_k B2(k,i)*X(k,j) = -(B2'*X)(i,j)
    let b2t_x = b2.transpose() * &*x;
    for i in 0..m2 {
        for j in 0..n {
            f[(i, j)] = -b2t_x[(i, j)];
        }
    }
    // H = -Y*C2' (N×NP2)
    let yc2 = &*y * c2.transpose();
    for i in 0..n {
        for j in 0..np2 {
            h[(i, j)] = -yc2[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10vd_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let mut f = DMatrix::zeros(1, 0);
        let mut h = DMatrix::zeros(0, 1);
        let mut x = DMatrix::zeros(0, 0);
        let mut y = DMatrix::zeros(0, 0);
        let mut xycond = [0.0; 2];
        assert_eq!(
            sb10vd(0, 2, 2, 1, 1, &a, &b, &c, &mut f, &mut h, &mut x, &mut y, &mut xycond),
            0
        );
    }

    #[test]
    fn test_sb10vd_1x1() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 2, &[0.0, 1.0]);
        let c = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let mut f = DMatrix::zeros(1, 1);
        let mut h = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        let mut y = DMatrix::zeros(1, 1);
        let mut xycond = [0.0; 2];
        let info = sb10vd(1, 2, 2, 1, 1, &a, &b, &c, &mut f, &mut h, &mut x, &mut y, &mut xycond);
        assert!(info == 0, "info = {}", info);
    }
}
