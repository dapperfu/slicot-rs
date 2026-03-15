//! Shared helpers for MB02 dense fallbacks.

use nalgebra::DMatrix;

/// Solves A*X = B, overwrites B with X. A is n×n, B is n×nrhs.
/// Returns 0 on success, 1 if singular, -1 if invalid dimensions.
pub fn solve_ax_b(n: usize, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || b.nrows() != n {
        return -1;
    }
    let lu = a.clone().lu();
    let bc = b.clone();
    let x = match lu.solve(&bc) {
        Some(s) => s,
        None => return 1,
    };
    b.copy_from(&x);
    0
}

/// Solves A'*X = B, overwrites B with X.
pub fn solve_atx_b(n: usize, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n || b.nrows() != n {
        return -1;
    }
    let lu = a.clone().lu();
    let bt = b.transpose();
    let y = match lu.solve(&bt) {
        Some(s) => s,
        None => return 1,
    };
    b.copy_from(&y.transpose());
    0
}
