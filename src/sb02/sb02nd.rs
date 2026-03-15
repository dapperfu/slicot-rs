//! SB02ND — Optimal state feedback matrix for LQR (SLICOT).
//!
//! Continuous: F = R^{-1}*(B'*X + L'). Discrete: F = (R+B'*X*B)^{-1}*(B'*X*A + L').

use nalgebra::DMatrix;

/// Continuous or discrete time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// How R is given (not used in this implementation; R is always full).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Fact {
    NotFactored,
    Cholesky,
    UdUorLdL,
}

/// Which triangle of R is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// Whether L is zero.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobL {
    Zero,
    Nonzero,
}

/// Computes optimal feedback F. R is M×M (symmetric, upper or lower). X is the CARE/DARE solution.
///
/// # Returns
/// 0 on success; &lt; 0 invalid argument; M+1 if R (or R+B'*X*B) singular.
pub fn sb02nd(
    dico: Dico,
    _fact: Fact,
    uplo: Uplo,
    jobl: JobL,
    n: usize,
    m: usize,
    a: &DMatrix<f64>,
    b: &mut DMatrix<f64>,
    r: &mut DMatrix<f64>,
    l: Option<&DMatrix<f64>>,
    x: &DMatrix<f64>,
    f: &mut DMatrix<f64>,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -8;
    }
    if b.nrows() != n || b.ncols() != m {
        return -10;
    }
    if r.nrows() != m || r.ncols() != m {
        return -12;
    }
    if x.nrows() != n || x.ncols() != n {
        return -16;
    }
    if f.nrows() != m || f.ncols() != n {
        return -18;
    }
    if n == 0 || m == 0 {
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

    let l_mat = match jobl {
        JobL::Zero => None,
        JobL::Nonzero => {
            let ell = match l {
                Some(e) => e,
                None => return -14,
            };
            if ell.nrows() != n || ell.ncols() != m {
                return -14;
            }
            Some(ell)
        }
    };

    let (coef, rhs_final) = match dico {
        Dico::Continuous => {
            let bt_x = b.transpose() * x;
            let rhs = match &l_mat {
                None => bt_x,
                Some(ell) => bt_x + ell.transpose(),
            };
            (r_full.clone(), rhs)
        }
        Dico::Discrete => {
            let btxb = b.transpose() * &*x * b.clone();
            let r_plus_btxb = &r_full + &btxb;
            let mut rhs_d = b.transpose() * x * a;
            if let Some(ell) = &l_mat {
                rhs_d += ell.transpose();
            }
            (r_plus_btxb, rhs_d)
        }
    };
    let ch = match coef.cholesky() {
        Some(c) => c,
        None => return (m + 1) as i32,
    };
    let z = ch.solve(&rhs_final);
    f.copy_from(&z.transpose());
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02nd_continuous() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut r = DMatrix::from_row_slice(1, 1, &[1.0]);
        let x = DMatrix::from_row_slice(1, 1, &[0.5]);
        let mut f = DMatrix::zeros(1, 1);
        assert_eq!(
            sb02nd(Dico::Continuous, Fact::NotFactored, Uplo::Upper, JobL::Zero, 1, 1, &a, &mut b, &mut r, None, &x, &mut f),
            0
        );
        assert!((f[(0, 0)] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_sb02nd_n0() {
        let a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut r = DMatrix::zeros(0, 0);
        let x = DMatrix::zeros(0, 0);
        let mut f = DMatrix::zeros(0, 0);
        assert_eq!(sb02nd(Dico::Continuous, Fact::NotFactored, Uplo::Upper, JobL::Zero, 0, 0, &a, &mut b, &mut r, None, &x, &mut f), 0);
    }
}
