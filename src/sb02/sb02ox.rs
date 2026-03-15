//! SB02OX — Extract diagonal of N×N matrix (SLICOT support).

use nalgebra::DMatrix;

/// Copies the diagonal of A into b: b(i) = A(i,i). b must have length at least n.
///
/// # Returns
/// 0 on success; &lt; 0 if invalid.
pub fn sb02ox(n: usize, a: &DMatrix<f64>, b: &mut [f64]) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n || b.len() < n {
        return -3;
    }
    for i in 0..n {
        b[i] = a[(i, i)];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02ox() {
        let a = DMatrix::from_row_slice(1, 1, &[7.0]);
        let mut b = [0.0];
        assert_eq!(sb02ox(1, &a, &mut b), 0);
        assert!((b[0] - 7.0).abs() < 1e-10);
    }
}
