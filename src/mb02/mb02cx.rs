//! MB02CX — Solve op(A)*X = B (SLICOT MB02). Dense fallback.

use nalgebra::DMatrix;

use crate::mb02::common;

/// Solves A*X = B (trans b"N") or A'*X = B (trans b"T"), overwrites B with X. Returns 0, 1 if singular, -1 if invalid.
pub fn mb02cx(trans: &[u8], n: usize, k: usize, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    if k == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || b.nrows() != n || b.ncols() != k {
        return -1;
    }
    if trans.first() == Some(&b'T') || trans.first() == Some(&b't') {
        common::solve_atx_b(n, a, b)
    } else {
        common::solve_ax_b(n, a, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb02cx_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb02cx(b"N", 0, 0, &a, &mut b), 0);
    }
}
