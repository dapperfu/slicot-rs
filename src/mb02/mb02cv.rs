//! MB02CV — Copy matrix block (SLICOT MB02). Dense fallback.

use nalgebra::DMatrix;

/// Copies A(0:n,0:n) to B(0:n,0:n). LDA/LDB ignored when using DMatrix. Returns 0, -1 if invalid.
pub fn mb02cv(n: usize, a: &DMatrix<f64>, _lda: usize, b: &mut DMatrix<f64>, _ldb: usize) -> i32 {
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
    fn test_mb02cv_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb02cv(0, &a, 0, &mut b, 0), 0);
    }
}
