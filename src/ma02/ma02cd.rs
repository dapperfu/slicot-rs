//! MA02CD — Pertranspose of central band of a square matrix (SLICOT MA02CD)
//
// Reverses the order of elements on each antidiagonal within the band
// (KL subdiagonals, main diagonal, KU superdiagonals). Equivalent to P*B'*P.

use nalgebra::DMatrix;

/// Pertransposes the central band (KL subdiagonals, main diagonal, KU superdiagonals) of A in place.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn ma02cd(a: &mut DMatrix<f64>, kl: usize, ku: usize) -> i32 {
    let n = a.nrows();
    if a.ncols() != n {
        return -4;
    }
    if n <= 1 {
        return 0;
    }

    // Pertranspose the KL subdiagonals (subdiagonal i has indices (i+j, j) for j=0..n-1-i)
    for i in 1..=kl.min(n.saturating_sub(2)) {
        let i1 = (n - i) / 2;
        if i1 > 0 {
            for j in 0..i1 {
                let r1 = i + j;
                let c1 = j;
                let r2 = n - 1 - i1 + j;
                let c2 = n - 1 - i1 - i + j;
                if r2 < n && c2 < n {
                    a.swap((r1, c1), (r2, c2));
                }
            }
        }
    }
    // Pertranspose the KU superdiagonals (superdiagonal i has indices (j, i+j) for j=0..n-1-i)
    for i in 1..=ku.min(n.saturating_sub(2)) {
        let i1 = (n - i) / 2;
        if i1 > 0 {
            for j in 0..i1 {
                let r1 = j;
                let c1 = i + j;
                let r2 = n - 1 - i1 - i + j;
                let c2 = n - 1 - i1 + j;
                if r2 < n && c2 < n {
                    a.swap((r1, c1), (r2, c2));
                }
            }
        }
    }
    // Pertranspose the main diagonal (first half with second half)
    let i1 = n / 2;
    if i1 > 0 {
        for j in 0..i1 {
            let r1 = j;
            let c1 = j;
            let r2 = n - i1 - j;
            let c2 = n - i1 - j;
            a.swap((r1, c1), (r2, c2));
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02cd_diag_only() {
        // 3x3: reverse main diagonal (0,0),(1,1),(2,2) -> swap (0,0)-(2,2), (1,1) stays
        let mut a = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]);
        assert_eq!(ma02cd(&mut a, 0, 0), 0);
        assert_eq!(a[(0, 0)], 3.0);
        assert_eq!(a[(1, 1)], 2.0);
        assert_eq!(a[(2, 2)], 1.0);
    }
}
