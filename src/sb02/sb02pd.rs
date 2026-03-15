//! SB02PD — Symmetric part of matrix: Q := (A + A')/2 (SLICOT support).

use nalgebra::DMatrix;

/// Computes the symmetric part: Q(i,j) = (A(i,j) + A(j,i))/2. A and Q are n×n.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02pd(n: usize, a: &DMatrix<f64>, q: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n || q.nrows() < n || q.ncols() < n {
        return -3;
    }
    for i in 0..n {
        for j in 0..n {
            q[(i, j)] = (a[(i, j)] + a[(j, i)]) * 0.5;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02pd() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 4.0, 2.0, 3.0]);
        let mut q = DMatrix::zeros(2, 2);
        assert_eq!(sb02pd(2, &a, &mut q), 0);
        assert!((q[(0, 1)] - 3.0).abs() < 1e-10);
        assert!((q[(1, 0)] - 3.0).abs() < 1e-10);
    }
}
