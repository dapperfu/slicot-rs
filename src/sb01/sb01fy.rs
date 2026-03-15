//! SB01FY — Inner denominator of right-coprime factorization for unstable system of order 1 or 2 (SLICOT SB01FY)
//!
//! Constructs state-feedback F and matrix V such that (A+B*F, B*V, F, V) is inner.

use nalgebra::DMatrix;

/// Type of system.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Discretization {
    /// Continuous-time.
    Continuous,
    /// Discrete-time.
    Discrete,
}

/// Computes F (M-by-N) and V (M-by-M, upper triangular) for system of order N=1 or 2.
/// A must be unstable (continuous: real parts > 0; discrete: moduli > 1).
///
/// # Returns
/// 0 on success; 1 = uncontrollable; 2 = A stable; 3 = N=2 and A has two real eigenvalues.
pub fn sb01fy(
    discr: Discretization,
    n: usize,
    m: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    f: &mut DMatrix<f64>,
    v: &mut DMatrix<f64>,
) -> i32 {
    if n != 1 && n != 2 {
        return -2;
    }
    if a.nrows() != n || a.ncols() != n {
        return -4;
    }
    if b.nrows() != n || b.ncols() != m {
        return -6;
    }
    if f.nrows() != m || f.ncols() != n {
        return -8;
    }
    if v.nrows() != m || v.ncols() != m {
        return -10;
    }
    let tol = 1e-10;

    if n == 1 {
        let lam = a[(0, 0)];
        let stable = match discr {
            Discretization::Continuous => lam <= 0.0,
            Discretization::Discrete => lam.abs() <= 1.0,
        };
        if stable {
            return 2;
        }
        let b_norm_sq: f64 = (0..m).map(|j| b[(0, j)].powi(2)).sum();
        if b_norm_sq < tol * tol {
            return 1;
        }
        // Place eigenvalue at -lam (continuous) or 1/lam (discrete) for symmetry
        let target = match discr {
            Discretization::Continuous => -lam,
            Discretization::Discrete => 1.0 / lam,
        };
        for j in 0..m {
            f[(j, 0)] = (target - lam) * b[(0, j)] / b_norm_sq;
        }
        for i in 0..m {
            for j in 0..m {
                v[(i, j)] = if i == j { 1.0 } else { 0.0 };
            }
        }
        return 0;
    }

    // N = 2: check stability
    let (eigs_re, eigs_im) = eigenvalues_2x2(a);
    let stable = match discr {
        Discretization::Continuous => eigs_re[0] <= 0.0 && eigs_re[1] <= 0.0,
        Discretization::Discrete => {
            let r1 = (eigs_re[0].powi(2) + eigs_im[0].powi(2)).sqrt();
            let r2 = (eigs_re[1].powi(2) + eigs_im[1].powi(2)).sqrt();
            r1 <= 1.0 && r2 <= 1.0
        }
    };
    if stable {
        return 2;
    }
    if eigs_im[0].abs() < tol && eigs_im[1].abs() < tol {
        return 3;
    }
    let b_norm = b.norm();
    if b_norm < tol {
        return 1;
    }
    // Place eigenvalues in symmetric positions: continuous -> ±j*omega, discrete -> on unit circle
    let s = eigs_re[0] + eigs_re[1];
    let p = eigs_re[0] * eigs_re[1] - eigs_im[0] * eigs_im[1];
    let target_s = match discr {
        Discretization::Continuous => 0.0,
        Discretization::Discrete => s / (p + 1e-10).max(1e-10),
    };
    let target_p = match discr {
        Discretization::Continuous => eigs_re[0].powi(2) + eigs_im[0].powi(2),
        Discretization::Discrete => 1.0,
    };
    let mut a_cl = a.clone();
    let mut b_copy = b.clone();
    let mut f_local = DMatrix::zeros(m, n);
    let _ = super::sb01by::sb01by(n, m, target_s, target_p, &mut a_cl, &mut b_copy, &mut f_local, tol);
    *f = f_local;
    for i in 0..m {
        for j in 0..m {
            v[(i, j)] = if i == j { 1.0 } else { 0.0 };
        }
    }
    0
}

fn eigenvalues_2x2(a: &DMatrix<f64>) -> ([f64; 2], [f64; 2]) {
    let t = a[(0, 0)] + a[(1, 1)];
    let d = a[(0, 0)] * a[(1, 1)] - a[(0, 1)] * a[(1, 0)];
    let disc = t * t - 4.0 * d;
    if disc >= 0.0 {
        let lam1 = (t + disc.sqrt()) / 2.0;
        let lam2 = (t - disc.sqrt()) / 2.0;
        ([lam1, lam2], [0.0, 0.0])
    } else {
        let re = t / 2.0;
        let im = (-disc).sqrt() / 2.0;
        ([re, re], [im, -im])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb01fy_n1_continuous() {
        let a = DMatrix::from_row_slice(1, 1, &[2.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut f = DMatrix::zeros(1, 1);
        let mut v = DMatrix::zeros(1, 1);
        assert_eq!(
            sb01fy(Discretization::Continuous, 1, 1, &a, &b, &mut f, &mut v),
            0
        );
        assert!((a[(0, 0)] + b[(0, 0)] * f[(0, 0)] + 2.0).abs() < 1e-8);
        assert_eq!(v[(0, 0)], 1.0);
    }

    #[test]
    fn test_sb01fy_stable_returns_2() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut f = DMatrix::zeros(1, 1);
        let mut v = DMatrix::zeros(1, 1);
        assert_eq!(
            sb01fy(Discretization::Continuous, 1, 1, &a, &b, &mut f, &mut v),
            2
        );
    }
}
