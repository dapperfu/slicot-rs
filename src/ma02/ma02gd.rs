//! MA02GD — Column interchanges (SLICOT MA02GD)
//
// Swaps columns K and IPIV(K) for K = K1..K2. Used after DGETRF for X*A = B.

use nalgebra::DMatrix;

/// Applies column interchanges: for k in k1..=k2, if ipiv[k] != k, swap columns k and ipiv[k].
/// Indices 0-based. ipiv.len() must be > k2.
pub fn ma02gd(a: &mut DMatrix<f64>, k1: usize, k2: usize, ipiv: &[i32]) -> i32 {
    let n = a.nrows();
    let m = a.ncols();
    if n == 0 {
        return 0;
    }
    for j in k1..=k2.min(m.saturating_sub(1)) {
        let jp = ipiv.get(j).copied().unwrap_or(j as i32) as usize;
        if j < jp && jp < m {
            for i in 0..n {
                a.swap((i, j), (i, jp));
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02gd_swap_cols() {
        let mut a = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let ipiv = [0, 2, 1]; // swap col 1 and 2
        assert_eq!(ma02gd(&mut a, 1, 2, &ipiv), 0);
        assert_eq!(a[(0, 1)], 3.0);
        assert_eq!(a[(0, 2)], 2.0);
    }
}
