//! MA02BZ — Reverse the order of rows and/or columns (SLICOT MA02BZ)
//
// Complex version in SLICOT; this is the real-matrix implementation (same logic as MA02BD).

use nalgebra::DMatrix;

/// Which side to reverse.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02BzSide {
    /// Reverse rows (pre-multiply by P).
    Left,
    /// Reverse columns (post-multiply by P).
    Right,
    /// Reverse both rows and columns.
    Both,
}

/// Reverses the order of rows and/or columns of A in place (P has ones on the anti-diagonal).
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn ma02bz(side: Ma02BzSide, a: &mut DMatrix<f64>) -> i32 {
    let m = a.nrows();
    let n = a.ncols();
    if m == 0 || n == 0 {
        return 0;
    }
    if side == Ma02BzSide::Left || side == Ma02BzSide::Both {
        let m2 = m / 2;
        for j in 0..n {
            for i in 0..m2 {
                let k = m - 1 - i;
                a.swap((i, j), (k, j));
            }
        }
    }
    if side == Ma02BzSide::Right || side == Ma02BzSide::Both {
        let n2 = n / 2;
        for i in 0..m {
            for j in 0..n2 {
                let k = n - 1 - j;
                a.swap((i, j), (i, k));
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02bz_both() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        assert_eq!(ma02bz(Ma02BzSide::Both, &mut a), 0);
        assert_eq!(a[(0, 0)], 4.0);
        assert_eq!(a[(1, 1)], 1.0);
    }
}
