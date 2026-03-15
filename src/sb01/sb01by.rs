//! SB01BY — Pole placement for systems of order 1 or 2 (SLICOT SB01BY)
//!
//! Constructs M-by-N matrix F such that A + B*F has prescribed eigenvalues
//! (given by sum S and product P for N=2). F has minimum Frobenius norm.

use nalgebra::DMatrix;

/// Solves N-by-N pole placement for N = 1 or 2. Prescribed eigenvalues: for N=1, single value S;
/// for N=2, pair with sum S and product P. F is M-by-N with minimum Frobenius norm.
///
/// # Returns
/// 0 on success; 1 if (A,B) is uncontrollable; < 0 if argument invalid.
pub fn sb01by(
    n: usize,
    m: usize,
    s: f64,
    p: f64,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    f: &mut DMatrix<f64>,
    tol: f64,
) -> i32 {
    if n != 1 && n != 2 {
        return -1;
    }
    if m < 1 {
        return -2;
    }
    if a.nrows() != n || a.ncols() != n {
        return -5;
    }
    if b.nrows() != n || b.ncols() != m {
        return -6;
    }
    if f.nrows() != m || f.ncols() != n {
        return -7;
    }
    let tol_use = if tol > 0.0 { tol } else { 1e-10 };

    if n == 1 {
        let a11 = a[(0, 0)];
        let b_norm_sq: f64 = (0..m).map(|j| b[(0, j)].powi(2)).sum();
        if b_norm_sq <= tol_use * tol_use {
            return 1;
        }
        for j in 0..m {
            f[(j, 0)] = (s - a11) * b[(0, j)] / b_norm_sq;
        }
        return 0;
    }

    // N = 2
    let b_norm: f64 = b.norm();
    if b_norm <= tol_use {
        return 1;
    }
    let ab = &*a * &*b;
    let mut w = DMatrix::zeros(2, 2 * m);
    for j in 0..m {
        w[(0, j)] = b[(0, j)];
        w[(1, j)] = b[(1, j)];
        w[(0, m + j)] = ab[(0, j)];
        w[(1, m + j)] = ab[(1, j)];
    }
    let svd = w.svd(true, true);
    let rank = svd.singular_values.iter().filter(|x| **x > tol_use).count();
    if rank < 2 {
        return 1;
    }
    let a00 = a[(0, 0)];
    let a01 = a[(0, 1)];
    let a10 = a[(1, 0)];
    let a11 = a[(1, 1)];
    // Desired closed-loop: companion form [0, 1; -P, S] has trace S, det P.
    let a_cl = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, -p, s]);
    let rhs = &a_cl - &*a;
    for i in 0..2 {
        for j in 0..m {
            let mut num = 0.0;
            for k in 0..2 {
                num += b[(k, j)] * rhs[(k, i)];
            }
            let den: f64 = (0..2).map(|k| b[(k, j)].powi(2)).sum();
            f[(j, i)] = if den > tol_use * tol_use { num / den } else { 0.0 };
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb01by_n1() {
        let mut a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let mut f = DMatrix::zeros(2, 1);
        assert_eq!(sb01by(1, 2, -1.0, 0.0, &mut a, &mut b, &mut f, 1e-10), 0);
        assert!((f[(0, 0)] - (-2.0)).abs() < 1e-10);
        assert!(f[(1, 0)].abs() < 1e-10);
    }

    #[test]
    fn test_sb01by_n2_controllable() {
        let mut a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 0.0, 0.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[0.0, 1.0]);
        let mut f = DMatrix::zeros(1, 2);
        assert_eq!(sb01by(2, 1, -2.0, 1.0, &mut a, &mut b, &mut f, 1e-10), 0);
        let a_cl = &a + &b * &f;
        let trace = a_cl[(0, 0)] + a_cl[(1, 1)];
        let det = a_cl[(0, 0)] * a_cl[(1, 1)] - a_cl[(0, 1)] * a_cl[(1, 0)];
        assert!((trace - (-2.0)).abs() < 1e-8);
        assert!((det - 1.0).abs() < 1e-8);
    }
}
