//! MA02SD — Compute a scale factor for a real matrix (SLICOT-style)
//
// No SLICOT doc found; provides scale = max(1, max_ij |A(i,j)|) so that A/scale has max norm <= 1.
// Empty or zero matrix => scale = 1.0.

use nalgebra::DMatrix;

/// Computes scale = max(1.0, max_ij |A(i,j)|). Returns 0 on success.
pub fn ma02sd(a: &DMatrix<f64>, scale: &mut f64) -> i32 {
    let m = a.nrows();
    let n = a.ncols();
    if m == 0 || n == 0 {
        *scale = 1.0;
        return 0;
    }
    let mut max_abs = 0.0_f64;
    for i in 0..m {
        for j in 0..n {
            max_abs = max_abs.max(a[(i, j)].abs());
        }
    }
    *scale = 1.0_f64.max(max_abs);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02sd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut scale = 0.0;
        assert_eq!(ma02sd(&a, &mut scale), 0);
        assert_eq!(scale, 1.0);
    }

    #[test]
    fn test_ma02sd_zero_matrix() {
        let a = DMatrix::<f64>::zeros(2, 2);
        let mut scale = 0.0;
        assert_eq!(ma02sd(&a, &mut scale), 0);
        assert_eq!(scale, 1.0);
    }

    #[test]
    fn test_ma02sd_nonzero() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, -4.0, 2.0, 3.0]);
        let mut scale = 0.0;
        assert_eq!(ma02sd(&a, &mut scale), 0);
        assert!((scale - 4.0).abs() < 1e-10);
    }
}
