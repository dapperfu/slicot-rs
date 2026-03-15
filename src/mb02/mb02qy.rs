//! MB02QY — Frobenius norm of matrix (SLICOT MB02). Dense fallback.

use nalgebra::DMatrix;

/// Sets *b to the Frobenius norm of A(0:n,0:n). Returns 0, -1 if invalid.
pub fn mb02qy(n: usize, a: &DMatrix<f64>, b: &mut f64) -> i32 {
    if n == 0 {
        *b = 0.0;
        return 0;
    }
    if a.nrows() < n || a.ncols() < n {
        return -1;
    }
    *b = (0..n)
        .flat_map(|i| (0..n).map(move |j| a[(i, j)] * a[(i, j)]))
        .sum::<f64>()
        .sqrt();
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb02qy_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = 0.0;
        assert_eq!(mb02qy(0, &a, &mut b), 0);
    }
}
