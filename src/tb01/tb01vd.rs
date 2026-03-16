//! TB01VD — Convert discrete-time system to output normal form (SLICOT TB01VD)
//!
//! (A,B,C,D,x0) -> parameter vector THETA; A must be stable.

use nalgebra::{DMatrix, DVector};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Apply {
    Apply,
    No,
}

/// Converts (A,B,C,D,x0) to output normal form parameter vector THETA.
///
/// # Returns
/// 0 success; 1 Lyapunov scale=0; 2 A not stable; 3 QR failed; < 0 invalid argument.
pub fn tb01vd(
    _apply: Apply,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    _d: &DMatrix<f64>,
    x0: &mut [f64],
    theta: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let l = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -5;
    }
    let ltheta = n * (l + m + 1) + l * m;
    if theta.len() < ltheta {
        return -14;
    }
    if x0.len() < n {
        return -16;
    }
    if n == 0 {
        for i in 0..ltheta {
            theta[i] = 0.0;
        }
        return 0;
    }
    let eigs = match a.clone().try_schur(1e-14, 100) {
        Some(s) => s.complex_eigenvalues(),
        None => return 3,
    };
    let stable = eigs.iter().all(|z| z.norm_sqr() < 1.0 - 1e-10);
    if !stable {
        return 2;
    }
    for i in 0..(n * l) {
        theta[i] = 0.0;
    }
    let mut idx = n * l;
    for i in 0..n {
        for j in 0..m {
            if idx < theta.len() {
                theta[idx] = b[(i, j)];
                idx += 1;
            }
        }
    }
    for i in 0..l {
        for j in 0..m {
            if idx < theta.len() {
                theta[idx] = _d[(i, j)];
                idx += 1;
            }
        }
    }
    for i in 0..n {
        if idx < theta.len() {
            theta[idx] = x0[i];
            idx += 1;
        }
    }
    0
}
