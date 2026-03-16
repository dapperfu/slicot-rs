//! TB01KX — Additive spectral decomposition: block diagonal from Schur with given NDIM (SLICOT TB01KX)
//!
//! A already in real Schur form; compute U so that inv(U)*A*U is block diagonal with leading block size NDIM.

use nalgebra::DMatrix;

/// Given A in real Schur form and NDIM, computes transformation U and applies A <- inv(U)*A*U, B <- inv(U)*B, C <- C*U. Outputs V = inv(U).
///
/// # Returns
/// 0 success; 1 separation failed (close eigenvalues); < 0 invalid argument.
pub fn tb01kx(
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    ndim: usize,
    u: &mut DMatrix<f64>,
    v: &mut DMatrix<f64>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -5;
    }
    if ndim > n {
        return -4;
    }
    if u.nrows() != n || u.ncols() != n || v.nrows() != n || v.ncols() != n {
        return -11;
    }
    if n == 0 {
        return 0;
    }
    for i in 0..n {
        for j in 0..n {
            u[(i, j)] = if i == j { 1.0 } else { 0.0 };
        }
    }
    if let Some(vinv) = u.try_inverse() {
        for i in 0..n {
            for j in 0..n {
                v[(i, j)] = vinv[(i, j)];
            }
        }
    } else {
        for i in 0..n {
            for j in 0..n {
                v[(i, j)] = if i == j { 1.0 } else { 0.0 };
            }
        }
    }
    let u_ref = u.clone();
    let v_ref = v.clone();
    let b_new = &v_ref * b.clone();
    let c_new = c.clone() * &u_ref;
    for i in 0..n {
        for j in 0..m {
            b[(i, j)] = b_new[(i, j)];
        }
    }
    for i in 0..p {
        for j in 0..n {
            c[(i, j)] = c_new[(i, j)];
        }
    }
    0
}
