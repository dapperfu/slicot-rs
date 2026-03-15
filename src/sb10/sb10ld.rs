//! SB10LD — Closed-loop system matrices from open-loop P and controller K.
//! Computes G = LFT(P,K): state [x; xk], inputs = exogenous (M-NCON), outputs = (NP-NMEAS).

use nalgebra::DMatrix;

/// Forms closed-loop system (AC, BC, CC, DC) from plant P and controller K.
/// Assumes well-posed: Inp2 - D22*DK and Im2 - DK*D22 nonsingular.
///
/// # Returns
/// 0 success; 1 Inp2 - D22*DK singular; 2 Im2 - DK*D22 singular; < 0 invalid argument.
pub fn sb10ld(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    ak: &DMatrix<f64>,
    bk: &DMatrix<f64>,
    ck: &DMatrix<f64>,
    dk: &DMatrix<f64>,
    ac: &mut DMatrix<f64>,
    bc: &mut DMatrix<f64>,
    cc: &mut DMatrix<f64>,
    dc: &mut DMatrix<f64>,
) -> i32 {
    let m1 = m - ncon;
    let np1 = np - nmeas;
    let m2 = ncon;
    let np2 = nmeas;
    if ac.nrows() != 2 * n || ac.ncols() != 2 * n {
        return -21;
    }
    if bc.nrows() != 2 * n || bc.ncols() != m1 {
        return -23;
    }
    if cc.nrows() != np1 || cc.ncols() != 2 * n {
        return -25;
    }
    if dc.nrows() != np1 || dc.ncols() != m1 {
        return -27;
    }

    let mut d22 = DMatrix::zeros(np2, m2);
    for i in 0..np2 {
        for j in 0..m2 {
            d22[(i, j)] = d[(np1 + i, m1 + j)];
        }
    }
    let inp2 = DMatrix::identity(np2, np2);
    let im2 = DMatrix::identity(m2, m2);
    let inp2_minus_d22_dk = &inp2 - &d22 * dk;
    let im2_minus_dk_d22 = &im2 - dk * &d22;
    let inv_inp2 = match inp2_minus_d22_dk.try_inverse() {
        Some(inv) => inv,
        None => return 1,
    };
    let inv_im2 = match im2_minus_dk_d22.try_inverse() {
        Some(inv) => inv,
        None => return 2,
    };

    let mut b1 = DMatrix::zeros(n, m1);
    let mut b2 = DMatrix::zeros(n, m2);
    for i in 0..n {
        for j in 0..m1 {
            b1[(i, j)] = b[(i, j)];
        }
        for j in 0..m2 {
            b2[(i, j)] = b[(i, m1 + j)];
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
    let mut d11 = DMatrix::zeros(np1, m1);
    let mut d12 = DMatrix::zeros(np1, m2);
    let mut d21 = DMatrix::zeros(np2, m1);
    for i in 0..np1 {
        for j in 0..m1 {
            d11[(i, j)] = d[(i, j)];
        }
        for j in 0..m2 {
            d12[(i, j)] = d[(i, m1 + j)];
        }
    }
    for i in 0..np2 {
        for j in 0..m1 {
            d21[(i, j)] = d[(np1 + i, j)];
        }
    }

    // LFT: x_dot = A*x + B1*w + B2*u, y = C1*x + D11*w + D12*u, y_meas = C2*x + D21*w + D22*u
    // u = CK*xk + DK*y_meas, xk_dot = AK*xk + BK*y_meas
    // y_meas = C2*x + D21*w + D22*u => u = CK*xk + DK*C2*x + DK*D21*w + DK*D22*u => (I-DK*D22)*u = CK*xk + DK*C2*x + DK*D21*w => u = inv(I-DK*D22)*(CK*xk + DK*C2*x + DK*D21*w)
    let dk_c2 = dk * &c2;
    let dk_d21 = dk * &d21;
    let u_coeff_xk = &inv_im2 * ck;
    let u_coeff_x = &inv_im2 * &dk_c2;
    let u_coeff_w = &inv_im2 * &dk_d21;

    let b2_u_coeff_xk = &b2 * &u_coeff_xk;
    let b2_u_coeff_x = &b2 * &u_coeff_x;
    let b2_u_coeff_w = &b2 * &u_coeff_w;

    // x_dot = A*x + B1*w + B2*u = A*x + B1*w + B2*u_coeff_xk*xk + B2*u_coeff_x*x + B2*u_coeff_w*w = (A + B2*u_coeff_x)*x + B2*u_coeff_xk*xk + (B1 + B2*u_coeff_w)*w
    // xk_dot = AK*xk + BK*y_meas = AK*xk + BK*(C2*x + D21*w + D22*u) = AK*xk + BK*C2*x + BK*D21*w + BK*D22*(u_coeff_xk*xk + u_coeff_x*x + u_coeff_w*w) = (AK + BK*D22*u_coeff_xk)*xk + (BK*C2 + BK*D22*u_coeff_x)*x + (BK*D21 + BK*D22*u_coeff_w)*w
    let bk_d22 = bk * &d22;
    let ac_11 = &*a + &b2_u_coeff_x;
    let ac_12 = &b2_u_coeff_xk;
    let ac_21 = &*bk * &c2 + &bk_d22 * &u_coeff_x;
    let ac_22 = &*ak + &bk_d22 * &u_coeff_xk;

    for i in 0..n {
        for j in 0..n {
            ac[(i, j)] = ac_11[(i, j)];
        }
        for j in 0..n {
            ac[(i, n + j)] = ac_12[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..n {
            ac[(n + i, j)] = ac_21[(i, j)];
        }
        for j in 0..n {
            ac[(n + i, n + j)] = ac_22[(i, j)];
        }
    }

    let bc_1 = &b1 + &b2_u_coeff_w;
    let bc_2 = &*bk * &d21 + &bk_d22 * &u_coeff_w;
    for i in 0..n {
        for j in 0..m1 {
            bc[(i, j)] = bc_1[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..m1 {
            bc[(n + i, j)] = bc_2[(i, j)];
        }
    }

    let cc_1 = &c1 + &d12 * &u_coeff_x;
    let cc_2 = &d12 * &u_coeff_xk;
    for i in 0..np1 {
        for j in 0..n {
            cc[(i, j)] = cc_1[(i, j)];
        }
        for j in 0..n {
            cc[(i, n + j)] = cc_2[(i, j)];
        }
    }

    let dc_mat = &d11 + &d12 * &u_coeff_w;
    for i in 0..np1 {
        for j in 0..m1 {
            dc[(i, j)] = dc_mat[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10ld_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let d = DMatrix::zeros(2, 2);
        let ak = DMatrix::zeros(0, 0);
        let bk = DMatrix::zeros(0, 1);
        let ck = DMatrix::zeros(1, 0);
        let dk = DMatrix::zeros(1, 1);
        let mut ac = DMatrix::zeros(0, 0);
        let mut bc = DMatrix::zeros(0, 1);
        let mut cc = DMatrix::zeros(1, 0);
        let mut dc = DMatrix::zeros(1, 1);
        assert_eq!(
            sb10ld(0, 2, 2, 1, 1, &a, &b, &c, &d, &ak, &bk, &ck, &dk, &mut ac, &mut bc, &mut cc, &mut dc),
            0
        );
    }
}
