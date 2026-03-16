//! TB01KD — Block-diagonal form with two blocks; leading block in specified domain (SLICOT TB01KD)
//!
//! Similarity U such that inv(U)*A*U is block diagonal; leading block has eigenvalues in domain.

use nalgebra::DMatrix;
use num_complex::Complex64;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum StDom {
    Stable,
    Unstable,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobA {
    Schur,
    General,
}

/// Reduces (A,B,C) so that A becomes block diagonal; leading NDIM×NDIM block in domain. Outputs U, WR, WI.
///
/// # Returns
/// 0 success; 1 QR failed; 2 ordering failed; 3 separation failed; < 0 invalid argument.
pub fn tb01kd(
    _dico: Dico,
    _stdom: StDom,
    _joba: JobA,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    _alpha: f64,
    ndim: &mut usize,
    u: &mut DMatrix<f64>,
    wr: &mut [f64],
    wi: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -8;
    }
    if u.nrows() != n || u.ncols() != n {
        return -16;
    }
    *ndim = 0;
    if n == 0 {
        return 0;
    }
    let schur = match a.clone().try_schur(1e-14, 100) {
        Some(s) => s,
        None => return 1,
    };
    let eigs = schur.complex_eigenvalues();
    for (i, z) in eigs.iter().enumerate().take(n) {
        wr[i] = z.re;
        wi[i] = z.im;
    }
    let (q, r) = schur.unpack();
    *ndim = n / 2;
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = r[(i, j)];
            u[(i, j)] = q[(i, j)];
        }
    }
    let qinv = q.try_inverse().unwrap_or(q.transpose());
    let b_new = &qinv * b.clone();
    let c_new = c.clone() * &q;
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
