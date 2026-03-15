//! SB02RU — Extract diagonal of N×N matrix (SLICOT support).

use nalgebra::DMatrix;

/// Copies the diagonal of A into x: x(i) = A(i,i). x must have length at least n.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02ru(n: usize, a: &DMatrix<f64>, x: &mut [f64]) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n || x.len() < n {
        return -3;
    }
    for i in 0..n {
        x[i] = a[(i, i)];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ru() {
        let a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut x = [0.0];
        assert_eq!(sb02ru(1, &a, &mut x), 0);
        assert!((x[0] - 1.0).abs() < 1e-10);
    }
}
