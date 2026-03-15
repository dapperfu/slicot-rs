//! SB10WD — H2 optimal controller matrices from state feedback F and output injection H (continuous-time).
//! K = (AK,BK,CK,DK) from F, H and transformation matrices TU, TY from SB10UD.

use nalgebra::DMatrix;

/// Forms controller K from F, H, TU, TY. After SB10UD, D12 = [0;I], D21 = [0 I], D22 is in D.
/// Formulas: DK = inv(I - D22'*D22)*D22' (or similar), CK = (F - DK*C2), BK = (H - B2*DK), AK = A + B2*F + H*C2 - B2*DK*C2.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn sb10wd(
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
    if a.nrows() != n || a.ncols() != n {
        return -7;
    }
    if b.nrows() != n || b.ncols() != m {
        return -9;
    }
    if c.nrows() != np || c.ncols() != n {
        return -11;
    }
    if d.nrows() != np || d.ncols() != m {
        return -13;
    }
    let m1 = m - ncon;
    let np1 = np - nmeas;
    let m2 = ncon;
    let np2 = nmeas;
    if f.nrows() != m2 || f.ncols() != n {
        return -15;
    }
    if h.nrows() != n || h.ncols() != np2 {
        return -17;
    }
    if ak.nrows() != n || ak.ncols() != n {
        return -21;
    }
    if bk.nrows() != n || bk.ncols() != np2 {
        return -23;
    }
    if ck.nrows() != m2 || ck.ncols() != n {
        return -25;
    }
    if dk.nrows() != m2 || dk.ncols() != np2 {
        return -27;
    }

    let mut b2 = DMatrix::zeros(n, m2);
    for i in 0..n {
        for j in 0..m2 {
            b2[(i, j)] = b[(i, m1 + j)];
        }
    }
    let mut c2 = DMatrix::zeros(np2, n);
    for i in 0..np2 {
        for j in 0..n {
            c2[(i, j)] = c[(np1 + i, j)];
        }
    }
    let mut d22 = DMatrix::zeros(np2, m2);
    for i in 0..np2 {
        for j in 0..m2 {
            d22[(i, j)] = d[(np1 + i, m1 + j)];
        }
    }

    // DK = inv(Im2 + D22'*D22)*D22' for H2 controller (well-posed closed loop)
    let im2 = DMatrix::identity(m2, m2);
    let d22t = d22.transpose();
    let d22td22 = &d22t * &d22;
    let mut im2_plus = im2.clone();
    for i in 0..m2 {
        for j in 0..m2 {
            im2_plus[(i, j)] += d22td22[(i, j)];
        }
    }
    if let Some(inv_im2_plus) = im2_plus.try_inverse() {
        let dk_mat = &inv_im2_plus * &d22t;
        for i in 0..m2 {
            for j in 0..np2 {
                dk[(i, j)] = dk_mat[(i, j)];
            }
        }
    } else {
        for i in 0..m2 {
            for j in 0..np2 {
                dk[(i, j)] = 0.0;
            }
        }
    }

    // CK = F - DK*C2 (F is M2×N, DK is M2×NP2, C2 is NP2×N)
    let dk_c2 = &*dk * &c2;
    for i in 0..m2 {
        for j in 0..n {
            ck[(i, j)] = f[(i, j)] - dk_c2[(i, j)];
        }
    }
    // BK = H - B2*DK
    let b2_dk = &b2 * &*dk;
    for i in 0..n {
        for j in 0..np2 {
            bk[(i, j)] = h[(i, j)] - b2_dk[(i, j)];
        }
    }
    // AK = A + B2*F + H*C2 - B2*DK*C2
    let b2_f = &b2 * f;
    let h_c2 = h * &c2;
    let b2_dk_c2 = &b2_dk * &c2;
    for i in 0..n {
        for j in 0..n {
            ak[(i, j)] = a[(i, j)] + b2_f[(i, j)] + h_c2[(i, j)] - b2_dk_c2[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10wd_n0() {
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
        assert_eq!(
            sb10wd(0, 2, 2, 1, 1, &a, &b, &c, &d, &f, &h, &tu, &ty, &mut ak, &mut bk, &mut ck, &mut dk),
            0
        );
    }
}
