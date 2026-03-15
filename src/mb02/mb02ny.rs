//! MB02NY — 1-norm of matrix (SLICOT MB02). Dense fallback.

use nalgebra::DMatrix;

/// Sets *b to the 1-norm (max column sum of abs) of A. Returns 0, -1 if invalid.
pub fn mb02ny(n: usize, a: &DMatrix<f64>, b: &mut f64) -> i32 {
    if n == 0 {
        *b = 0.0;
        return 0;
    }
    if a.nrows() < n || a.ncols() < n {
        return -1;
    }
    *b = (0..n)
        .map(|j| (0..n).map(|i| a[(i, j)].abs()).sum::<f64>())
        .fold(0.0f64, f64::max);
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb02ny_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = 0.0;
        assert_eq!(mb02ny(0, &a, &mut b), 0);
    }
}
