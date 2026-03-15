//! SG02CW — Residual of continuous- or discrete-time (generalized) algebraic Riccati equations (SLICOT SG02).
//!
//! Computes R = residual matrix and/or C = closed-loop matrix from
//! R = op(A)'*X*op(E) + op(E)'*X*op(A) ± op(E)'*X*G*X*op(E) + Q (continuous)
//! or discrete-time variants.

use nalgebra::DMatrix;

/// DICO: continuous or discrete.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// JOB: compute R only, C only, both, or norms too.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Job {
    Both,
    ROnly,
    COnly,
    Norms,
    RAndNorms,
}

/// JOBE: E general or identity.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobE {
    General,
    Identity,
}

/// FLAG: plus or minus in quadratic term.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Flag {
    Plus,
    Minus,
}

/// JOBG: how quadratic term is given (G, D, F, or H,K).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobG {
    G,
    D,
    F,
    H,
}

/// UPLO: upper or lower triangle of symmetric matrices.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// TRANS: op(W) = W or W'.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Trans {
    NoTrans,
    Trans,
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

/// Computes residual R and/or closed-loop matrix C.
///
/// # Returns
/// 0 on success; < 0 invalid argument index.
pub fn sg02cw(
    dico: Dico,
    job: Job,
    jobe: JobE,
    flag: Flag,
    jobg: JobG,
    uplo: Uplo,
    trans: Trans,
    n: usize,
    m: usize,
    a: &DMatrix<f64>,
    e: Option<&DMatrix<f64>>,
    g: &DMatrix<f64>,
    x: &DMatrix<f64>,
    f: Option<&DMatrix<f64>>,
    k: Option<&DMatrix<f64>>,
    _xe: Option<&DMatrix<f64>>,
    r: &mut DMatrix<f64>,
    c: Option<&mut DMatrix<f64>>,
    norms: Option<&mut [f64]>,
    _dwork: &mut [f64],
    _ldwork: i32,
    info: &mut i32,
) -> i32 {
    *info = 0;
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || x.nrows() != n || x.ncols() != n {
        *info = -8;
        return -8;
    }
    let trans_a = trans == Trans::Trans;
    let sign = if flag == Flag::Plus { 1.0 } else { -1.0 };

    let x_full = sym_full(n, x, uplo);
    let q_full = if job != Job::COnly {
        sym_full(n, r, uplo)
    } else {
        DMatrix::zeros(n, n)
    };

    let (need_r, need_c) = match job {
        Job::ROnly | Job::RAndNorms => (true, false),
        Job::COnly => (false, true),
        Job::Both | Job::Norms => (true, true),
    };

    let e_mat = if jobe == JobE::General {
        e.map(|em| sym_full(n, em, uplo))
    } else {
        None
    };

    let (op_a, op_e) = if trans_a {
        (a.transpose(), e_mat.as_ref().map(|em| em.transpose()))
    } else {
        (a.clone(), e_mat.as_ref().cloned())
    };

    let op_a_x = &op_a * &x_full;
    let op_a_x_op_e = if let Some(ref em) = op_e {
        &op_a_x * em
    } else {
        op_a_x.clone()
    };
    let op_e_x = if let Some(ref em) = op_e {
        em * &x_full
    } else {
        x_full.clone()
    };

    if need_r {
        let mut res = op_a_x_op_e.clone() + op_a_x_op_e.transpose();
        if dico == Dico::Discrete {
            let op_e_x_op_e = if let Some(ref em) = op_e {
                em * &op_e_x
            } else {
                x_full.clone()
            };
            res -= op_e_x_op_e;
        }

        let quad_term = match jobg {
            JobG::G => {
                let g_full = sym_full(n, g, uplo);
                if dico == Dico::Continuous {
                    if op_e.is_some() {
                        op_e_x.transpose() * &g_full * &op_e_x
                    } else {
                        &x_full * &g_full * &x_full
                    }
                } else {
                    op_a_x.transpose() * &g_full * &op_a_x
                }
            }
            JobG::D => {
                let g_eff = g * g.transpose();
                if dico == Dico::Continuous {
                    if op_e.is_some() {
                        op_e_x.transpose() * &g_eff * &op_e_x
                    } else {
                        &x_full * &g_eff * &x_full
                    }
                } else {
                    op_a_x.transpose() * &g_eff * &op_a_x
                }
            }
            JobG::F => {
                let ff = f.unwrap_or(g);
                ff * ff.transpose()
            }
            JobG::H => {
                let hh = f.unwrap_or(g);
                let kk = match k {
                    Some(kmat) => kmat,
                    None => {
                        *info = -18;
                        return -18;
                    }
                };
                hh * kk
            }
        };
        res += sign * quad_term;
        res += &q_full;

        for i in 0..n {
            for j in 0..n {
                r[(i, j)] = res[(i, j)];
            }
        }
        if uplo == Uplo::Upper {
            for i in 0..n {
                for j in 0..i {
                    r[(i, j)] = res[(j, i)];
                }
            }
        } else {
            for i in 0..n {
                for j in (i + 1)..n {
                    r[(i, j)] = res[(j, i)];
                }
            }
        }
    }

    if need_c {
        if let Some(cmat) = c {
            let gx_term = match jobg {
                JobG::G => {
                    let g_full = sym_full(n, g, uplo);
                    &g_full * &x_full
                }
                JobG::D => (g * g.transpose()) * &x_full,
                JobG::F => {
                    let ff = f.unwrap_or(g);
                    ff * (ff.transpose() * &x_full)
                }
                JobG::H => {
                    let kk = match k {
                        Some(kmat) => kmat,
                        None => {
                            *info = -18;
                            return -18;
                        }
                    };
                    g * kk
                }
            };
            let clos = &op_a + sign * &gx_term;
            cmat.copy_from(&clos);
        }
    }

    if let Some(ns) = norms {
        if job == Job::Norms || job == Job::RAndNorms {
            let mut n1 = 0.0_f64;
            for i in 0..n {
                for j in 0..n {
                    n1 += op_a_x_op_e[(i, j)].powi(2);
                }
            }
            ns[0] = n1.sqrt();
            if ns.len() > 1 {
                ns[1] = 0.0;
            }
        }
    }
    0
}
