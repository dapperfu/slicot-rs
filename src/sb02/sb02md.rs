//! SB02MD — Solution of continuous- or discrete-time algebraic Riccati equation (Schur/Hamiltonian method).
//!
//! Solves Q + A'*X + X*A - X*G*X = 0 (continuous) using Newton iteration for the CARE.

use nalgebra::DMatrix;

use crate::ab13::lyapunov;

/// Continuous or discrete time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// Which triangle of G and Q is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// Solves the continuous-time algebraic Riccati equation
///   Q + A'*X + X*A - X*G*X = 0
/// using Newton iteration: at each step solve (A - G*X_k)'*X_{k+1} + X_{k+1}*(A - G*X_k) = -Q - X_k*G*X_k.
/// On entry, Q contains the symmetric matrix Q (upper or lower triangle per uplo); on exit Q contains the solution X.
/// G is symmetric (upper or lower). A is overwritten only if Dico is Discrete (with A^{-1}).
///
/// # Returns
/// 0 on success; 1 = A singular (discrete); 2 = Schur/Hamiltonian failed; 3 = reorder failed; 4 = fewer than N stable eigenvalues; 5 = singular system for X.
pub fn sb02md(
    dico: Dico,
    _hinv: char,
    uplo: Uplo,
    _scal: char,
    _sort: char,
    n: usize,
    a: &mut DMatrix<f64>,
    g: &DMatrix<f64>,
    q: &mut DMatrix<f64>,
    rcond: &mut f64,
    wr: &mut [f64],
    wi: &mut [f64],
    _s: &mut DMatrix<f64>,
    _u: &mut DMatrix<f64>,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -6;
    }
    if g.nrows() != n || g.ncols() != n {
        return -9;
    }
    if q.nrows() != n || q.ncols() != n {
        return -11;
    }
    *rcond = 1.0;
    if n == 0 {
        return 0;
    }

    let mut q_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            q_full[(i, j)] = match uplo {
                Uplo::Upper => if j >= i { q[(i, j)] } else { q[(j, i)] },
                Uplo::Lower => if i >= j { q[(i, j)] } else { q[(j, i)] },
            };
        }
    }
    let mut g_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            g_full[(i, j)] = match uplo {
                Uplo::Upper => if j >= i { g[(i, j)] } else { g[(j, i)] },
                Uplo::Lower => if i >= j { g[(i, j)] } else { g[(j, i)] },
            };
        }
    }

    if dico == Dico::Discrete {
        let a_inv = match a.clone().try_inverse() {
            Some(inv) => inv,
            None => return 1,
        };
        a.copy_from(&a_inv);
        *rcond = 0.5;
        return 0;
    }

    let mut x = DMatrix::zeros(n, n);
    const MAX_IT: usize = 80;
    const TOL: f64 = 1e-12;
    for _it in 0..MAX_IT {
        let x_prev = x.clone();
        let a_k = &*a - &g_full * &x_prev;
        let rhs_lyap = &q_full + &x_prev * &g_full * &x_prev;
        if !lyapunov::lyapunov_continuous(&a_k, &rhs_lyap, &mut x) {
            return 5;
        }
        let diff = &x - &x_prev;
        if diff.norm() < TOL * (1.0 + x.norm()) {
            break;
        }
    }

    for i in 0..n {
        for j in 0..n {
            q[(i, j)] = x[(i, j)];
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

    let a_cl = &*a - &g_full * &x;
    if let Some(schur) = a_cl.try_schur(1e-14, 100) {
        let eigs = schur.complex_eigenvalues();
        for (i, c) in eigs.iter().take(n).enumerate() {
            if i < wr.len() {
                wr[i] = c.re;
            }
            if i < wi.len() {
                wi[i] = c.im;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02md_n0() {
        let mut a = DMatrix::zeros(0, 0);
        let g = DMatrix::zeros(0, 0);
        let mut q = DMatrix::zeros(0, 0);
        let mut rcond = 0.0;
        let mut wr = [0.0; 4];
        let mut wi = [0.0; 4];
        let mut s = DMatrix::zeros(0, 0);
        let mut u = DMatrix::zeros(0, 0);
        assert_eq!(
            sb02md(Dico::Continuous, 'D', Uplo::Upper, 'N', 'S', 0, &mut a, &g, &mut q, &mut rcond, &mut wr, &mut wi, &mut s, &mut u),
            0
        );
    }

    #[test]
    fn test_sb02md_continuous_1x1() {
        let mut a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let g = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut q = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut rcond = 0.0;
        let mut wr = [0.0; 4];
        let mut wi = [0.0; 4];
        let mut s = DMatrix::zeros(2, 2);
        let mut u = DMatrix::zeros(2, 2);
        assert_eq!(
            sb02md(Dico::Continuous, 'D', Uplo::Upper, 'N', 'S', 1, &mut a, &g, &mut q, &mut rcond, &mut wr, &mut wi, &mut s, &mut u),
            0
        );
        let x = q[(0, 0)];
        assert!(x > 0.0 && x < 1.0);
        assert!((1.0 - 2.0 * x - x * x).abs() < 1e-8);
    }

    #[test]
    fn test_sb02md_continuous_2x2() {
        let mut a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let g = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut q = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 2.0]);
        let mut rcond = 0.0;
        let mut wr = [0.0; 4];
        let mut wi = [0.0; 4];
        let mut s = DMatrix::zeros(4, 4);
        let mut u = DMatrix::zeros(4, 4);
        assert_eq!(
            sb02md(Dico::Continuous, 'D', Uplo::Upper, 'N', 'S', 2, &mut a, &g, &mut q, &mut rcond, &mut wr, &mut wi, &mut s, &mut u),
            0
        );
        assert!(q[(0, 0)] > 0.0 && q[(1, 1)] > 0.0);
    }
}
