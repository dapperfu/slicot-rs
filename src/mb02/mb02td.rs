//! MB02TD — Estimate RCOND for upper Hessenberg H factored by MB02SD (SLICOT).
//!
//! RCOND = 1 / (||H||_1 * ||inv(H)||_1).

use nalgebra::DMatrix;
use std::f64;

use crate::mb02::mb02rd::mb02rd;

/// Computes RCOND. H is the N×N factored matrix from MB02SD; ipiv from MB02SD.
/// Returns 0 and sets rcond, or >0 if singular (rcond = 0), or <0 if invalid input.
pub fn mb02td(n: usize, h: &DMatrix<f64>, ipiv: &[i32], rcond: &mut f64) -> i32 {
    if n == 0 {
        *rcond = 1.0;
        return 0;
    }
    if h.nrows() != n || h.ncols() != n || ipiv.len() < n {
        return -1;
    }

    // 1-norm of H = P'*L*U is same as 1-norm of L*U
    let norm_h = one_norm_lu(n, h);
    if norm_h == 0.0 {
        *rcond = 0.0;
        return 0;
    }

    // X = inv(H) by solving H*X = I
    let mut ident = DMatrix::identity(n, n);
    let info = mb02rd(n, h, ipiv, &mut ident);
    if info != 0 {
        *rcond = 0.0;
        return info;
    }
    let norm_inv = one_norm_matrix(&ident);
    *rcond = if norm_inv == 0.0 {
        0.0
    } else {
        1.0 / (norm_h * norm_inv).min(f64::MAX)
    };
    0
}

/// 1-norm of L*U stored in h (L lower bidiagonal, U upper).
fn one_norm_lu(n: usize, h: &DMatrix<f64>) -> f64 {
    let mut max_col = 0.0f64;
    for j in 0..n {
        let mut col: Vec<f64> = vec![0.0; n];
        for k in 0..=j {
            col[k] = h[(k, j)];
        }
        for i in 1..n {
            col[i] = h[(i, i - 1)] * col[i - 1] + col[i];
        }
        let sum_abs: f64 = col.iter().map(|x| x.abs()).sum();
        max_col = max_col.max(sum_abs);
    }
    max_col
}

fn one_norm_matrix(a: &DMatrix<f64>) -> f64 {
    (0..a.ncols())
        .map(|j| (0..a.nrows()).map(|i| a[(i, j)].abs()).sum::<f64>())
        .fold(0.0f64, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mb02::mb02sd::mb02sd_matrix;

    #[test]
    fn test_mb02td_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let ipiv: Vec<i32> = vec![];
        let mut rcond = 0.0;
        assert_eq!(mb02td(0, &a, &ipiv, &mut rcond), 0);
        assert_eq!(rcond, 1.0);
    }

    #[test]
    fn test_mb02td_2x2() {
        let mut h = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 1.0, 3.0]);
        let mut ipiv = vec![0i32; 2];
        assert_eq!(mb02sd_matrix(&mut h, &mut ipiv), 0);
        let mut rcond = 0.0;
        assert_eq!(mb02td(2, &h, &ipiv, &mut rcond), 0);
        assert!(rcond > 0.0 && rcond <= 1.0);
    }
}