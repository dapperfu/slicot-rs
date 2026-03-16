//! ZGEGV — Generalized complex eigenvalue problem (SLICOT/LAPACK auxiliary).
//!
//! Solves A*x = lambda*B*x. Returns eigenvalues as (alpha, beta) with lambda = alpha/beta.
//! When B is nonsingular, eigenvalues are those of B^{-1}*A (Schur decomposition).

use nalgebra::linalg::LU;
use nalgebra::DMatrix;
use num_complex::Complex64;

/// Computes generalized eigenvalues of the pencil (A, B). Eigenvalues are returned as
/// alpha[i]/beta[i]. When B is nonsingular, these are the eigenvalues of B^{-1}*A.
///
/// # Arguments
/// - `a`, `b`: N×N matrices (only leading N×N is used).
/// - `alpha`: length N, filled with numerator of eigenvalue.
/// - `beta`: length N, filled with denominator (0 => infinite eigenvalue).
///
/// # Returns
/// 0 = success; 1 = B singular or numerical failure; < 0 = invalid argument.
pub fn zgegv(
    n: usize,
    a: &DMatrix<Complex64>,
    b: &DMatrix<Complex64>,
    alpha: &mut [Complex64],
    beta: &mut [Complex64],
) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n {
        return -2;
    }
    if b.nrows() < n || b.ncols() < n {
        return -3;
    }
    if alpha.len() < n || beta.len() < n {
        return -4;
    }
    let a_sub = a.view((0, 0), (n, n));
    let b_sub = b.view((0, 0), (n, n));
    let b_mat = b_sub.into_owned();
    let lu = LU::new(b_mat.clone());
    let b_inv = match lu.try_inverse() {
        Some(inv) => inv,
        None => return 1,
    };
    let a_owned = a_sub.into_owned();
    let c = &b_inv * &a_owned;
    let schur = match c.try_schur(1e-14, 200) {
        Some(s) => s,
        None => return 1,
    };
    let (_, t) = schur.unpack();
    for i in 0..n {
        alpha[i] = t[(i, i)];
        beta[i] = Complex64::new(1.0, 0.0);
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zgegv_n_zero() {
        let a = DMatrix::from_element(1, 1, Complex64::new(1.0, 0.0));
        let b = DMatrix::from_element(1, 1, Complex64::new(1.0, 0.0));
        let mut alpha = [Complex64::new(0.0, 0.0)];
        let mut beta = [Complex64::new(0.0, 0.0)];
        assert_eq!(zgegv(0, &a, &b, &mut alpha, &mut beta), 0);
    }

    #[test]
    fn test_zgegv_1x1() {
        let a = DMatrix::from_element(1, 1, Complex64::new(2.0, 0.0));
        let b = DMatrix::from_element(1, 1, Complex64::new(1.0, 0.0));
        let mut alpha = [Complex64::new(0.0, 0.0)];
        let mut beta = [Complex64::new(0.0, 0.0)];
        assert_eq!(zgegv(1, &a, &b, &mut alpha, &mut beta), 0);
        assert!((alpha[0].re - 2.0).abs() < 1e-10);
        assert!(beta[0].re.abs() - 1.0 < 1e-10);
    }

    #[test]
    fn test_zgegv_singular_b() {
        let a = DMatrix::from_element(2, 2, Complex64::new(1.0, 0.0));
        let b = DMatrix::from_element(2, 2, Complex64::new(0.0, 0.0));
        let mut alpha = [Complex64::new(0.0, 0.0); 2];
        let mut beta = [Complex64::new(0.0, 0.0); 2];
        assert_eq!(zgegv(2, &a, &b, &mut alpha, &mut beta), 1);
    }
}
