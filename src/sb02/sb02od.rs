//! SB02OD — Solution of continuous/discrete-time algebraic Riccati equation (general) (SLICOT).
//!
//! Driver that uses SB02MT (convert to standard form) and SB02MD (Newton iteration for CARE).

use nalgebra::DMatrix;

use super::sb02md::{sb02md, Dico as MdDico, Uplo as MdUplo};
use super::sb02mt::{Fact as MtFact, JobG, JobL as MtJobL, Uplo as MtUplo};

/// Solves the CARE (or DARE placeholder) via SB02MT + SB02MD.
///
/// JOBB='B': B and R given; we form G = B*R^{-1}*B' and optionally A_bar, Q_bar, then solve.
/// JOBB='G': B is G (n×n), Q is Q; solve directly.
/// FACT: 'N' = Q and R given; 'C' = Q = C'*C (q is P×N C); 'D' = R = D'*D (r is P×M D); 'B' = both.
/// JOBL: 'Z' = L zero, 'N' = L nonzero (l required).
///
/// # Returns
/// 0 success; &lt; 0 invalid argument; 1–6 from SB02MD.
pub fn sb02od(
    dico: char,
    jobb: char,
    fact: &str,
    uplo: char,
    jobl: char,
    _sort: char,
    n: usize,
    m: usize,
    p: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    q: &DMatrix<f64>,
    r: &DMatrix<f64>,
    l: Option<&DMatrix<f64>>,
    x: &mut DMatrix<f64>,
    _tol: f64,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -10;
    }
    if x.nrows() != n || x.ncols() != n {
        return -18;
    }
    let jobb_b = jobb == 'B' || jobb == 'b';
    let uplo_u = uplo == 'U' || uplo == 'u';
    let mt_uplo = if uplo_u { MtUplo::Upper } else { MtUplo::Lower };
    let md_uplo = if uplo_u { MdUplo::Upper } else { MdUplo::Lower };
    let dico_c = dico == 'C' || dico == 'c';
    let md_dico = if dico_c { MdDico::Continuous } else { MdDico::Discrete };

    if jobb_b {
        if b.nrows() != n || b.ncols() != m {
            return -12;
        }
        let (q_full, r_full) = match fact {
            "N" | "n" => {
                if q.nrows() != n || q.ncols() != n {
                    return -15;
                }
                if r.nrows() != m || r.ncols() != m {
                    return -16;
                }
                let mut qm = DMatrix::zeros(n, n);
                for i in 0..n {
                    for j in 0..n {
                        qm[(i, j)] = if uplo_u {
                            if j >= i { q[(i, j)] } else { q[(j, i)] }
                        } else {
                            if i >= j { q[(i, j)] } else { q[(j, i)] }
                        };
                    }
                }
                let mut rm = DMatrix::zeros(m, m);
                for i in 0..m {
                    for j in 0..m {
                        rm[(i, j)] = if uplo_u {
                            if j >= i { r[(i, j)] } else { r[(j, i)] }
                        } else {
                            if i >= j { r[(i, j)] } else { r[(j, i)] }
                        };
                    }
                }
                (qm, rm)
            }
            "C" | "c" => {
                if q.nrows() != p || q.ncols() != n {
                    return -15;
                }
                let c = q;
                let qm = c.transpose() * c;
                if r.nrows() != m || r.ncols() != m {
                    return -16;
                }
                let mut rm = DMatrix::zeros(m, m);
                for i in 0..m {
                    for j in 0..m {
                        rm[(i, j)] = if uplo_u {
                            if j >= i { r[(i, j)] } else { r[(j, i)] }
                        } else {
                            if i >= j { r[(i, j)] } else { r[(j, i)] }
                        };
                    }
                }
                (qm, rm)
            }
            "D" | "d" => {
                if q.nrows() != n || q.ncols() != n {
                    return -15;
                }
                let mut qm = DMatrix::zeros(n, n);
                for i in 0..n {
                    for j in 0..n {
                        qm[(i, j)] = if uplo_u {
                            if j >= i { q[(i, j)] } else { q[(j, i)] }
                        } else {
                            if i >= j { q[(i, j)] } else { q[(j, i)] }
                        };
                    }
                }
                if r.nrows() != p || r.ncols() != m {
                    return -16;
                }
                let d = r;
                let rm = d.transpose() * d;
                (qm, rm)
            }
            "B" | "b" | "Both" => {
                if q.nrows() != p || q.ncols() != n {
                    return -15;
                }
                if r.nrows() != p || r.ncols() != m {
                    return -16;
                }
                let c = q;
                let d = r;
                let qm = c.transpose() * c;
                let rm = d.transpose() * d;
                (qm, rm)
            }
            _ => return -4,
        };

        let jobl_n = jobl == 'N' || jobl == 'n';
        let mt_jobl = if jobl_n { MtJobL::Nonzero } else { MtJobL::Zero };
        let mut a_work = a.clone();
        let mut b_work = b.clone();
        let mut q_work = q_full;
        let mut r_work = r_full;
        let mut l_work = if jobl_n {
            match l {
                Some(ell) => {
                    if ell.nrows() != n || ell.ncols() != m {
                        return -17;
                    }
                    ell.clone()
                }
                None => return -17,
            }
        } else {
            DMatrix::zeros(n, m)
        };
        let mut oufact = 0i32;
        let mut g = DMatrix::zeros(n, n);
        let info_mt = super::sb02mt::sb02mt(
            JobG::Compute,
            mt_jobl,
            MtFact::NotFactored,
            mt_uplo,
            n,
            m,
            &mut a_work,
            &mut b_work,
            &mut q_work,
            &mut r_work,
            &mut l_work,
            &mut oufact,
            &mut g,
        );
        if info_mt != 0 {
            return info_mt;
        }
        let mut rcond = 0.0;
        let mut wr = vec![0.0; 2 * n];
        let mut wi = vec![0.0; 2 * n];
        let mut s = DMatrix::zeros(2 * n, 2 * n);
        let mut u = DMatrix::zeros(2 * n, 2 * n);
        let info_md = sb02md(
            md_dico,
            'D',
            md_uplo,
            'N',
            'S',
            n,
            &mut a_work,
            &g,
            &mut q_work,
            &mut rcond,
            &mut wr,
            &mut wi,
            &mut s,
            &mut u,
        );
        if info_md != 0 {
            return info_md;
        }
        x.copy_from(&q_work);
        return 0;
    }

    // JOBB = 'G': B is G (n×n), Q is Q (n×n)
    if b.nrows() != n || b.ncols() != n {
        return -12;
    }
    if q.nrows() != n || q.ncols() != n {
        return -15;
    }
    let mut g_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            g_full[(i, j)] = if uplo_u {
                if j >= i { b[(i, j)] } else { b[(j, i)] }
            } else {
                if i >= j { b[(i, j)] } else { b[(j, i)] }
            };
        }
    }
    let mut q_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            q_full[(i, j)] = if uplo_u {
                if j >= i { q[(i, j)] } else { q[(j, i)] }
            } else {
                if i >= j { q[(i, j)] } else { q[(j, i)] }
            };
        }
    }
    let mut a_work = a.clone();
    let mut q_work = q_full;
    let mut rcond = 0.0;
    let mut wr = vec![0.0; 2 * n];
    let mut wi = vec![0.0; 2 * n];
    let mut s = DMatrix::zeros(2 * n, 2 * n);
    let mut u = DMatrix::zeros(2 * n, 2 * n);
    let info_md = sb02md(
        md_dico,
        'D',
        md_uplo,
        'N',
        'S',
        n,
        &mut a_work,
        &g_full,
        &mut q_work,
        &mut rcond,
        &mut wr,
        &mut wi,
        &mut s,
        &mut u,
    );
    if info_md != 0 {
        return info_md;
    }
    x.copy_from(&q_work);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02od_jobb_b_fact_n_jobl_z() {
        // Stable A so Newton iteration converges
        let a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let b = DMatrix::from_row_slice(2, 1, &[0.0, 1.0]);
        let q = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let r = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut x = DMatrix::zeros(2, 2);
        assert_eq!(
            sb02od('C', 'B', "N", 'U', 'Z', 'S', 2, 1, 1, &a, &b, &q, &r, None, &mut x, 0.0),
            0
        );
        assert!(x[(0, 0)] > 0.0 && x[(1, 1)] > 0.0);
    }

    #[test]
    fn test_sb02od_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 0);
        let q = DMatrix::zeros(0, 0);
        let r = DMatrix::zeros(0, 0);
        let mut x = DMatrix::zeros(0, 0);
        assert_eq!(
            sb02od('C', 'B', "N", 'U', 'Z', 'S', 0, 0, 0, &a, &b, &q, &r, None, &mut x, 0.0),
            0
        );
    }
}
