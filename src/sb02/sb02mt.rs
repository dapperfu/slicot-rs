//! SB02MT — Compute G = B*R^{-1}*B', optionally A_bar = A - B*R^{-1}*L', Q_bar = Q - L*R^{-1}*L' (SLICOT).

use nalgebra::DMatrix;

/// Which triangle of symmetric matrices is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// Whether to compute G.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobG {
    Compute,
    No,
}

/// Whether L is zero.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobL {
    Zero,
    Nonzero,
}

/// How R is given.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Fact {
    NotFactored,
    Cholesky,
    UdUorLdL,
}

/// Computes G = B*R^{-1}*B', and optionally A_bar = A - B*R^{-1}*L', Q_bar = Q - L*R^{-1}*L'.
/// R is M×M (upper or lower triangle per uplo). On exit OUFACT indicates Cholesky (1) or UdU/LdL (2) used.
///
/// # Returns
/// 0 on success; &lt; 0 invalid argument; i (1..=M) if d factor zero; M+1 if R numerically singular.
pub fn sb02mt(
    jobg: JobG,
    jobl: JobL,
    fact: Fact,
    uplo: Uplo,
    n: usize,
    m: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    q: &mut DMatrix<f64>,
    r: &mut DMatrix<f64>,
    l: &mut DMatrix<f64>,
    oufact: &mut i32,
    g: &mut DMatrix<f64>,
) -> i32 {
    if b.nrows() != n || b.ncols() != m {
        return -10;
    }
    if r.nrows() != m || r.ncols() != m {
        return -14;
    }
    if jobg == JobG::Compute && (g.nrows() != n || g.ncols() != n) {
        return -18;
    }
    if jobl == JobL::Nonzero && (a.nrows() != n || a.ncols() != n) {
        return -6;
    }
    if jobl == JobL::Nonzero && (q.nrows() != n || q.ncols() != n) {
        return -12;
    }
    if jobl == JobL::Nonzero && (l.nrows() != n || l.ncols() != m) {
        return -16;
    }
    if m == 0 {
        *oufact = 0;
        if jobg == JobG::Compute {
            g.fill(0.0);
        }
        return 0;
    }

    let mut r_full = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            r_full[(i, j)] = match uplo {
                Uplo::Upper => if j >= i { r[(i, j)] } else { r[(j, i)] },
                Uplo::Lower => if i >= j { r[(i, j)] } else { r[(j, i)] },
            };
        }
    }

    let ch = match r_full.cholesky() {
        Some(c) => c,
        None => {
            *oufact = 2;
            return (m + 1) as i32;
        }
    };
    *oufact = 1;

    // Solve R*X = B' => X = R^{-1}*B' (X is M×N)
    let bt = b.transpose();
    let x = ch.solve(&bt);
    // G = B*X (N×N) = B*R^{-1}*B'
    if jobg == JobG::Compute {
        let g_mat = &*b * &x;
        for i in 0..n {
            for j in 0..n {
                g[(i, j)] = g_mat[(i, j)];
            }
        }
        if uplo == Uplo::Upper {
            for i in 0..n {
                for j in 0..i {
                    g[(i, j)] = g[(j, i)];
                }
            }
        } else {
            for i in 0..n {
                for j in (i + 1)..n {
                    g[(i, j)] = g[(j, i)];
                }
            }
        }
    }

    if jobl == JobL::Nonzero {
        // A_bar = A - B*R^{-1}*L' = A - B*X_l where R*X_l = L' => X_l = R^{-1}*L'
        let lt = l.transpose();
        let x_l = ch.solve(&lt);
        let bl = &*b * &x_l;
        for i in 0..n {
            for j in 0..n {
                a[(i, j)] -= bl[(i, j)];
            }
        }
        // Q_bar = Q - L*R^{-1}*L' = Q - L*X_l
        let lxl = &*l * &x_l;
        for i in 0..n {
            for j in 0..n {
                q[(i, j)] -= lxl[(i, j)];
            }
        }
        if uplo == Uplo::Upper {
            for i in 0..n {
                for j in 0..i {
                    q[(i, j)] = q[(j, i)];
                }
            }
        } else {
            for i in 0..n {
                for j in (i + 1)..n {
                    q[(i, j)] = q[(j, i)];
                }
            }
        }
        // B := B*chol(R)^{-1} (solve L'*Y = B' => Y = (L')^{-1}*B' => B_new = Y')
        let lt = ch.l().transpose();
        let lu_lt = lt.lu();
        let bt = b.transpose();
        if let Some(b_new_t) = lu_lt.solve(&bt) {
            b.copy_from(&b_new_t.transpose());
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02mt_m0() {
        let mut a = DMatrix::zeros(2, 2);
        let mut b = DMatrix::zeros(2, 0);
        let mut q = DMatrix::zeros(2, 2);
        let mut r = DMatrix::zeros(0, 0);
        let mut l = DMatrix::zeros(2, 0);
        let mut oufact = -1;
        let mut g = DMatrix::zeros(2, 2);
        assert_eq!(
            sb02mt(JobG::Compute, JobL::Zero, Fact::NotFactored, Uplo::Upper, 2, 0, &mut a, &mut b, &mut q, &mut r, &mut l, &mut oufact, &mut g),
            0
        );
        assert_eq!(oufact, 0);
    }

    #[test]
    fn test_sb02mt_g_basic() {
        let mut a = DMatrix::zeros(1, 1);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut q = DMatrix::zeros(1, 1);
        let mut r = DMatrix::from_row_slice(1, 1, &[4.0]); // R = 4 => R^{-1} = 0.25
        let mut l = DMatrix::zeros(1, 1);
        let mut oufact = -1;
        let mut g = DMatrix::zeros(1, 1);
        assert_eq!(
            sb02mt(JobG::Compute, JobL::Zero, Fact::NotFactored, Uplo::Upper, 1, 1, &mut a, &mut b, &mut q, &mut r, &mut l, &mut oufact, &mut g),
            0
        );
        assert_eq!(oufact, 1);
        // G = B*R^{-1}*B' = 1 * 0.25 * 1 = 0.25
        assert!((g[(0, 0)] - 0.25).abs() < 1e-10);
    }
}
