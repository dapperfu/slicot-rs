//! SG02CV — Copy matrix block for descriptor Riccati (SLICOT SG02).
//!
//! Copies the leading N×N block of matrix A into B. Used as an auxiliary
//! in descriptor Riccati and related routines.

use nalgebra::DMatrix;

/// Copies the leading N×N part of A into the leading N×N part of B.
/// LDA/LDB are ignored when using DMatrix (full matrix storage).
///
/// # Returns
/// 0 on success; -1 if dimensions are invalid (n out of range or A/B too small).
pub fn sg02cv(
    n: usize,
    a: &DMatrix<f64>,
    _lda: usize,
    b: &mut DMatrix<f64>,
    _ldb: usize,
) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n || b.nrows() < n || b.ncols() < n {
        return -1;
    }
    for i in 0..n {
        for j in 0..n {
            b[(i, j)] = a[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_sg02cv_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(sg02cv(0, &a, 0, &mut b, 0), 0);
    }

    #[test]
    fn test_sg02cv_copy() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let mut b = DMatrix::<f64>::zeros(2, 2);
        assert_eq!(sg02cv(2, &a, 2, &mut b, 2), 0);
        assert_eq!(b[(0, 0)], 1.0);
        assert_eq!(b[(1, 0)], 3.0);
        assert_eq!(b[(0, 1)], 2.0);
        assert_eq!(b[(1, 1)], 4.0);
    }

    #[test]
    fn test_sg02cv_invalid() {
        let a = DMatrix::<f64>::zeros(2, 2);
        let mut b = DMatrix::<f64>::zeros(1, 1);
        assert_eq!(sg02cv(2, &a, 2, &mut b, 1), -1);
    }
}
