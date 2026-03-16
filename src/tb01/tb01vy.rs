//! TB01VY — Convert output normal form to state-space (SLICOT TB01VY)
//!
//! THETA -> (A,B,C,D,x0).

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Apply {
    Apply,
    No,
}

/// Builds (A,B,C,D,x0) from output normal form parameter vector THETA.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01vy(
    _apply: Apply,
    theta: &[f64],
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    x0: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let l = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n || d.nrows() != l || d.ncols() != m {
        return -5;
    }
    let ltheta = n * (l + m + 1) + l * m;
    if theta.len() < ltheta || x0.len() < n {
        return -5;
    }
    if n == 0 {
        return 0;
    }
    let mut idx = 0;
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = if i == j { 1.0 } else { 0.0 };
        }
    }
    idx = n * l;
    for i in 0..n {
        for j in 0..m {
            if idx < theta.len() {
                b[(i, j)] = theta[idx];
                idx += 1;
            }
        }
    }
    for i in 0..l {
        for j in 0..m {
            if idx < theta.len() {
                d[(i, j)] = theta[idx];
                idx += 1;
            }
        }
    }
    for i in 0..n {
        if idx < theta.len() {
            x0[i] = theta[idx];
            idx += 1;
        }
    }
    0
}
