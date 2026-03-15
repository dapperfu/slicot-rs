//! SG02AD — Solution of continuous- or discrete-time algebraic Riccati equations for descriptor systems (SLICOT SG02).
//!
//! Solves for X the continuous-time equation
//!   Q + A'XE + E'XA - (L+E'XB)R^{-1}(L+E'XB)' = 0
//! or discrete-time equation
//!   E'XE = A'XA - (L+A'XB)(R+B'XB)^{-1}(L+A'XB)' + Q
//! using the method of deflating subspaces (generalized Schur).

use nalgebra::DMatrix;

/// DICO: continuous or discrete.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// JOBB: B and R given, or G = B*R^{-1}*B' given.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobb {
    B,
    G,
}

/// UPLO: upper or lower triangle of symmetric matrices.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// SORT: stable or unstable eigenvalues first.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Sort {
    StableFirst,
    UnstableFirst,
}

fn sym_full(n: usize, a: &DMatrix<f64>, uplo: Uplo) -> DMatrix<f64> {
    let mut f = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            f[(i, j)] = match uplo {
                Uplo::Upper => if j >= i { a[(i, j)] } else { a[(j, i)] },
                Uplo::Lower => if i >= j { a[(i, j)] } else { a[(j, i)] },
            };
        }
    }
    f
}

/// Solves the descriptor Riccati equation. Supports JOBB='G' path (G given) and JOBB='B' with L=0.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1 singular pencil; 2 QZ failed; 3 reorder failed; 4 roundoff changed eigenvalues; 5 dimension not N; 6 spectrum too close to boundary; 7 singular when solving for X.
pub fn sg02ad(
    dico: Dico,
    jobb: Jobb,
    _fact: &str,
    uplo: Uplo,
    _jobl: char,
    _scal: char,
    sort: Sort,
    _acc: char,
    n: usize,
    m: usize,
    _p: usize,
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &DMatrix<f64>,
    q: &DMatrix<f64>,
    r: Option<&DMatrix<f64>>,
    _l: Option<&DMatrix<f64>>,
    rcondu: &mut f64,
    x: &mut DMatrix<f64>,
    alfar: &mut [f64],
    alfai: &mut [f64],
    beta: &mut [f64],
    _s: &mut DMatrix<f64>,
    _t: &mut DMatrix<f64>,
    u: &mut DMatrix<f64>,
    _tol: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
    _bwork: &mut [bool],
    _iwarn: &mut i32,
    info: &mut i32,
) -> i32 {
    *info = 0;
    *rcondu = 1.0;
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || e.nrows() != n || e.ncols() != n {
        *info = -10;
        return -10;
    }
    if q.nrows() != n || q.ncols() != n {
        *info = -17;
        return -17;
    }
    if x.nrows() != n || x.ncols() != n {
        *info = -22;
        return -22;
    }

    let q_full = sym_full(n, q, uplo);
    let e_full = e.clone();
    let a_full = a.clone();

    let (m_pencil, n_pencil) = if jobb == Jobb::G {
        let g_full = sym_full(n, b, uplo);
        let m_mat = DMatrix::from_fn(2 * n, 2 * n, |i, j| {
            if i < n && j < n {
                a_full[(i, j)]
            } else if i < n && j >= n {
                0.0
            } else if i >= n && j < n {
                q_full[(i - n, j)]
            } else {
                -e_full[(j - n, i - n)]
            }
        });
        let mut n_mat = DMatrix::from_fn(2 * n, 2 * n, |i, j| {
            if i < n && j < n {
                e_full[(i, j)]
            } else if i < n && j >= n {
                0.0
            } else if i >= n && j < n {
                0.0
            } else {
                -a_full[(j - n, i - n)]
            }
        });
        if dico == Dico::Continuous {
            let gxe = &g_full * &e_full;
            for i in 0..n {
                for j in 0..n {
                    n_mat[(n + i, n + j)] -= gxe[(i, j)];
                }
            }
        }
        (m_mat, n_mat)
    } else {
        if m == 0 {
            *info = 1;
            return 1;
        }
        let r_mat = match r {
            Some(rr) => rr,
            None => {
                *info = -20;
                return -20;
            }
        };
        let r_full = sym_full(m, r_mat, uplo);
        let r_inv = match r_full.try_inverse() {
            Some(inv) => inv,
            None => {
                *info = 1;
                return 1;
            }
        };
        let g_full = b * &r_inv * b.transpose();
        let m_mat = DMatrix::from_fn(2 * n, 2 * n, |i, j| {
            if i < n && j < n {
                a_full[(i, j)]
            } else if i < n && j >= n {
                0.0
            } else if i >= n && j < n {
                q_full[(i - n, j)]
            } else {
                -e_full[(j - n, i - n)]
            }
        });
        let mut n_mat = DMatrix::from_fn(2 * n, 2 * n, |i, j| {
            if i < n && j < n {
                e_full[(i, j)]
            } else if i < n && j >= n {
                0.0
            } else if i >= n && j < n {
                0.0
            } else {
                -a_full[(j - n, i - n)]
            }
        });
        if dico == Dico::Continuous {
            let gxe = &g_full * &e_full;
            for i in 0..n {
                for j in 0..n {
                    n_mat[(n + i, n + j)] -= gxe[(i, j)];
                }
            }
        }
        (m_mat, n_mat)
    };

    let n_inv = match n_pencil.try_inverse() {
        Some(inv) => inv,
        None => {
            *info = 1;
            return 1;
        }
    };
    let standard = n_inv * &m_pencil;
    let schur = match standard.try_schur(1e-14, 200) {
        Some(s) => s,
        None => {
            *info = 2;
            return 2;
        }
    };
    let eigs = schur.complex_eigenvalues();
    let (q_schur, t_schur) = schur.unpack();
    let mut order: Vec<usize> = (0..2 * n).collect();
    let stable = |i: usize| {
        let z = eigs[i];
        if dico == Dico::Continuous {
            z.re < 0.0 || (z.re.abs() < 1e-14 && z.im == 0.0)
        } else {
            z.norm_sqr() < 1.0 - 1e-10
        }
    };
    order.sort_by(|&i, &j| {
        let si = stable(i);
        let sj = stable(j);
        match (sort == Sort::StableFirst, si, sj) {
            (true, true, false) => std::cmp::Ordering::Less,
            (true, false, true) => std::cmp::Ordering::Greater,
            (false, true, false) => std::cmp::Ordering::Greater,
            (false, false, true) => std::cmp::Ordering::Less,
            _ => i.cmp(&j),
        }
    });
    let mut u_ordered = DMatrix::zeros(2 * n, 2 * n);
    for (j, &idx) in order.iter().enumerate() {
        for i in 0..2 * n {
            u_ordered[(i, j)] = q_schur[(i, idx)];
        }
    }
    let u1 = u_ordered.view((0, 0), (n, n));
    let u2 = u_ordered.view((n, 0), (n, n));
    let u1_mat = DMatrix::from_fn(n, n, |i, j| u1[(i, j)]);
    let u2_mat = DMatrix::from_fn(n, n, |i, j| u2[(i, j)]);
    let u1_inv = match u1_mat.try_inverse() {
        Some(inv) => inv,
        None => {
            *info = 7;
            return 7;
        }
    };
    let sol = &u2_mat * &u1_inv;
    for i in 0..n {
        for j in 0..n {
            x[(i, j)] = sol[(i, j)];
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let s = (x[(i, j)] + x[(j, i)]) / 2.0;
            x[(i, j)] = s;
            x[(j, i)] = s;
        }
    }
    for (i, &idx) in order.iter().take(2 * n).enumerate() {
        if i < alfar.len() {
            alfar[i] = eigs[idx].re;
        }
        if i < alfai.len() {
            alfai[i] = eigs[idx].im;
        }
        if i < beta.len() {
            beta[i] = 1.0;
        }
    }
    for i in 0..2 * n {
        for j in 0..2 * n {
            u[(i, j)] = u_ordered[(i, j)];
        }
    }
    let _ = t_schur;
    0
}
