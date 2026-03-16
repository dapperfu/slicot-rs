//! TB01XZ — Dual transformation P*A'*P, P*C', B'*P (complex case) (SLICOT TB01XZ)

use nalgebra::DMatrix;
use num_complex::Complex64;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobD {
    Present,
    Zero,
}

/// Applies dual transformation to complex (A,B,C,D): A <- P*A'*P, B <- P*C', C <- B'*P, optionally D <- D'.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01xz(
    jobd: JobD,
    a: &mut DMatrix<Complex64>,
    b: &mut DMatrix<Complex64>,
    c: &mut DMatrix<Complex64>,
    d: Option<&mut DMatrix<Complex64>>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -6;
    }
    if jobd == JobD::Present && d.is_none() {
        return -14;
    }
    if let Some(ref dd) = d {
        if jobd == JobD::Present && (dd.nrows() != p || dd.ncols() != m) {
            return -14;
        }
    }
    if n == 0 {
        if jobd == JobD::Present {
            if let Some(ref mut dd) = d {
                let dc = dd.clone();
                for j in 0..m {
                    for i in 0..p {
                        dd[(j, i)] = dc[(i, j)];
                    }
                }
            }
        }
        return 0;
    }
    let ac = a.clone();
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = ac[(n - 1 - j, n - 1 - i)];
        }
    }
    let cc = c.clone();
    let b_orig = b.clone();
    for i in 0..n {
        for k in 0..p.min(b.ncols()) {
            b[(i, k)] = cc[(k, n - 1 - i)];
        }
    }
    for k in 0..m.min(c.nrows()) {
        for j in 0..n {
            c[(k, j)] = b_orig[(n - 1 - j, k)];
        }
    }
    if jobd == JobD::Present {
        if let Some(ref mut dd) = d {
            let dc = dd.clone();
            for j in 0..m {
                for i in 0..p {
                    dd[(j, i)] = dc[(i, j)];
                }
            }
        }
    }
    0
}
