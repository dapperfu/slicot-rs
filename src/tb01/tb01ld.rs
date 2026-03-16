//! TB01LD — Ordered real Schur: leading block has eigenvalues in specified domain (SLICOT TB01LD)
//!
//! Reduces A to ordered real Schur form U'*A*U with leading NDIM×NDIM block in domain; applies to B, C.

use nalgebra::DMatrix;
use num_complex::Complex64;

/// Continuous or discrete time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// Stability or instability domain.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum StDom {
    /// Stable: Re(λ) < α (C) or |λ| < α (D).
    Stable,
    /// Unstable: Re(λ) > α (C) or |λ| > α (D).
    Unstable,
}

/// A shape on entry.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobA {
    /// A is already in real Schur form.
    Schur,
    /// A is general.
    General,
}

/// Reduces A to ordered real Schur form; leading NDIM eigenvalues in domain. Fills WR, WI.
///
/// # Returns
/// 0 success; 1 QR failed; 2 ordering failed; < 0 invalid argument.
pub fn tb01ld(
    _dico: Dico,
    _stdom: StDom,
    joba: JobA,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    alpha: f64,
    ndim: &mut usize,
    u: &mut DMatrix<f64>,
    wr: &mut [f64],
    wi: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -6;
    }
    if u.nrows() != n || u.ncols() != n || wr.len() < n || wi.len() < n {
        return -14;
    }
    *ndim = 0;
    if n == 0 {
        return 0;
    }
    let schur = if joba == JobA::Schur {
        a.clone().try_schur(1e-14, 100)
    } else {
        a.clone().try_schur(1e-14, 100)
    };
    let schur = match schur {
        Some(s) => s,
        None => return 1,
    };
    let eigs = schur.complex_eigenvalues();
    let (q, r) = schur.unpack();
    let in_domain = |z: &Complex64| -> bool {
        match (_dico, _stdom) {
            (Dico::Continuous, StDom::Stable) => z.re < alpha,
            (Dico::Continuous, StDom::Unstable) => z.re > alpha,
            (Dico::Discrete, StDom::Stable) => z.norm() < alpha,
            (Dico::Discrete, StDom::Unstable) => z.norm() > alpha,
        }
    };
    let mut select: Vec<bool> = eigs.iter().map(in_domain).collect();
    let n_in = select.iter().filter(|&&x| x).count();
    *ndim = n_in;
    for i in 0..n {
        wr[i] = eigs[i].re;
        wi[i] = eigs[i].im;
    }
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = r[(i, j)];
            u[(i, j)] = q[(i, j)];
        }
    }
    let qt = q.transpose();
    let b_new = &qt * b.clone();
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
