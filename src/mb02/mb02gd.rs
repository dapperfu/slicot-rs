//! MB02GD — Cholesky factorization for banded/symmetric positive definite (SLICOT).
//!
//! Dense fallback: A = R'*R, overwrites b with upper Cholesky R.

use nalgebra::DMatrix;

/// Dense fallback: computes upper Cholesky R of symmetric positive definite A (A = R'*R).
/// a is N×N, b is N×N; on success b contains R. Returns 0, or 1 if A is not positive definite.
pub fn mb02gd(n: usize, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || b.nrows() != n || b.ncols() != n {
        return -1;
    }
    let ch = match a.clone().cholesky() {
        Some(c) => c,
        None => return 1,
    };
    let l = ch.l();
    for i in 0..n {
        for j in 0..n {
            b[(i, j)] = if j >= i { l[(j, i)] } else { 0.0 };
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb02gd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb02gd(0, &a, &mut b), 0);
    }

    #[test]
    fn test_mb02gd_2x2() {
        let a = DMatrix::from_row_slice(2, 2, &[4.0, 2.0, 2.0, 3.0]);
        let mut b = DMatrix::zeros(2, 2);
        assert_eq!(mb02gd(2, &a, &mut b), 0);
        let rtr = b.transpose() * &b;
        assert!((rtr[(0, 0)] - 4.0).abs() < 1e-10);
        assert!((rtr[(1, 1)] - 3.0).abs() < 1e-10);
    }
}
