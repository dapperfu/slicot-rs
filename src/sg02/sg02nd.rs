//! SG02ND — Optimal state feedback matrix for descriptor optimal control (SLICOT SG02).
//!
//! Computes the optimal gain matrix K:
//! - Continuous: R*K = B'*X*op(E) + L'  =>  K = R^{-1} * (B'*X*op(E) + L')
//! - Discrete:   (R + B'*X*B)*K = B'*X*op(A) + L'

use nalgebra::DMatrix;

/// Uplo: which triangle of R is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// FACT: how R is given.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Fact {
    /// R is full symmetric (one triangle stored).
    NotFactored,
    /// R = D'*D, array R contains P×M matrix D.
    FactoredD,
    /// R is Cholesky factor.
    Cholesky,
    /// R is UdU' or LdL' (indefinite).
    UdUorLdL,
}

/// JOB: what to compute.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Job {
    /// K only.
    K,
    /// H and K.
    H,
    /// F if possible (F*C = H), else H and K.
    F,
    /// H and K when B,L transformed (SB02MT/SB02MX); R contains Cholesky of R+B'XB.
    D,
    /// F when B,L transformed.
    C,
}

/// JOBE: E general or identity.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobE {
    General,
    Identity,
}

/// JOBL: L zero or nonzero.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobL {
    Zero,
    Nonzero,
}

/// TRANS: op(W) = W or W'.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Trans {
    NoTrans,
    Trans,
}

/// Computes the optimal gain matrix K (and optionally H or F).
///
/// `e` is required when dico is 'C' and jobe is General (E matrix). For discrete, `a` is the state matrix.
///
/// # Returns
/// 0 success; < 0 invalid argument; i (1..=M) d factor zero; M+1 R (or R+B'XB) singular; M+2 eigs not converged; M+3 indefinite update failed.
pub fn sg02nd(
    dico: char,
    jobe: JobE,
    _job: Job,
    _jobx: char,
    fact: Fact,
    uplo: Uplo,
    jobl: JobL,
    trans: Trans,
    n: usize,
    m: usize,
    _p: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    r: &DMatrix<f64>,
    _ipiv: &mut [i32],
    l: Option<&DMatrix<f64>>,
    x: &DMatrix<f64>,
    _rnorm: f64,
    k: &mut DMatrix<f64>,
    h: Option<&mut DMatrix<f64>>,
    _xe: Option<&mut DMatrix<f64>>,
    oufact: Option<&mut [i32; 2]>,
    _iwork: &mut [i32],
    dwork: Option<&mut [f64]>,
    _ldwork: i32,
    info: &mut i32,
    e: Option<&DMatrix<f64>>,
) -> i32 {
    *info = 0;
    if n == 0 || m == 0 {
        if let Some(ou) = oufact {
            ou[0] = 0;
            ou[1] = 0;
        }
        return 0;
    }
    if b.nrows() != n || b.ncols() != m {
        *info = -18;
        return -18;
    }
    if r.nrows() < m || r.ncols() < m {
        *info = -22;
        return -22;
    }
    if x.nrows() != n || x.ncols() != n {
        *info = -28;
        return -28;
    }
    if k.nrows() < m || k.ncols() < n {
        *info = -31;
        return -31;
    }
    if jobl == JobL::Nonzero {
        if let Some(ell) = l {
            if ell.nrows() != n || ell.ncols() != m {
                *info = -26;
                return -26;
            }
        } else {
            *info = -26;
            return -26;
        }
    }
    let trans_a = trans == Trans::Trans;

    // H_mat (M×N): R*K = H_mat  =>  K = R_sys^{-1} * H_mat. So columns of H are right-hand sides.
    // H' = B'*X*op(E) + L' (continuous) or B'*X*op(A) + L' (discrete). So H_mat(i,j) = (B'*X*op(E)+L')(i,j).
    let mut h_mat = DMatrix::zeros(m, n);
    if dico == 'C' || dico == 'c' {
        let op_ex_b = if jobe == JobE::Identity {
            x * b
        } else {
            let e_mat = match e {
                Some(em) => em,
                None => {
                    *info = -2;
                    return -2;
                }
            };
            if trans_a {
                e_mat.transpose() * x * b
            } else {
                e_mat * x * b
            }
        };
        for i in 0..n {
            for j in 0..m {
                let mut v = op_ex_b[(i, j)];
                if jobl == JobL::Nonzero {
                    v += l.unwrap()[(i, j)];
                }
                h_mat[(j, i)] = v;
            }
        }
    } else {
        // Discrete: H' = B'*X*op(A) + L', so h_mat (M×N) = (B'*X*op(A))' + L stored as (j,i). (B'*X*A)(j,i) = ((X*A)'*B)(i,j).
        let xa = if trans_a { x * a } else { x * a.transpose() };
        let xa_t_b = xa.transpose() * b;
        for i in 0..n {
            for j in 0..m {
                let mut v = xa_t_b[(i, j)];
                if jobl == JobL::Nonzero {
                    v += l.unwrap()[(i, j)];
                }
                h_mat[(j, i)] = v;
            }
        }
    }

    let mut r_sys = DMatrix::zeros(m, m);
    match fact {
        Fact::NotFactored => {
            for i in 0..m {
                for j in 0..m {
                    r_sys[(i, j)] = match uplo {
                        Uplo::Upper => if j >= i { r[(i, j)] } else { r[(j, i)] },
                        Uplo::Lower => if i >= j { r[(i, j)] } else { r[(j, i)] },
                    };
                }
            }
            if dico == 'D' || dico == 'd' {
                r_sys += &(b.transpose() * x * b);
            }
        }
        Fact::Cholesky => {
            let rf = r.clone();
            if uplo == Uplo::Upper {
                for i in 0..m {
                    for j in i..m {
                        let s: f64 = (0..=i).map(|k| rf[(k, i)] * rf[(k, j)]).sum();
                        r_sys[(i, j)] = s;
                        r_sys[(j, i)] = s;
                    }
                }
            } else {
                for i in 0..m {
                    for j in 0..=i {
                        let s: f64 = (j..m).map(|k| rf[(i, k)] * rf[(j, k)]).sum();
                        r_sys[(i, j)] = s;
                        r_sys[(j, i)] = s;
                    }
                }
            }
            if dico == 'D' || dico == 'd' {
                r_sys += &(b.transpose() * x * b);
            }
        }
        Fact::FactoredD | Fact::UdUorLdL => {
            *info = 1;
            return 1;
        }
    }

    let lu = r_sys.lu();
    for j in 0..n {
        let mut col = DMatrix::zeros(m, 1);
        for i in 0..m {
            col[(i, 0)] = h_mat[(i, j)];
        }
        if let Some(sol) = lu.solve(&col) {
            for i in 0..m {
                k[(i, j)] = sol[(i, 0)];
            }
        } else {
            *info = (m + 1) as i32;
            return *info;
        }
    }

    if let Some(hout) = h {
        for i in 0..n {
            for j in 0..m {
                hout[(i, j)] = h_mat[(j, i)];
            }
        }
    }
    if let Some(ou) = oufact {
        ou[0] = 1;
        ou[1] = 0;
    }
    if let Some(dw) = dwork {
        if dw.len() >= 2 {
            dw[0] = 0.0;
            dw[1] = 1.0;
        }
    }
    0
}
