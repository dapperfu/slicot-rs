//! AB05PD — Addition of two system state-space models (SLICOT AB05PD).
//!
//! Computes G = G1 + alpha*G2: A = diag(A1,A2), B = [B1;B2],
//! C = [C1, alpha*C2], D = D1 + alpha*D2. This implementation uses alpha = 1.

use nalgebra::DMatrix;

/// Addition of two systems: G = G1 + G2 (alpha=1).
/// Returns 0 on success; < 0 invalid argument index.
pub fn ab05pd(
    n1: usize,
    m1: usize,
    p1: usize,
    n2: usize,
    a1: &DMatrix<f64>,
    b1: &DMatrix<f64>,
    c1: &DMatrix<f64>,
    d1: &DMatrix<f64>,
    a2: &DMatrix<f64>,
    b2: &DMatrix<f64>,
    c2: &DMatrix<f64>,
    d2: &DMatrix<f64>,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
) -> i32 {
    let n = n1 + n2;
    if n1 > 0 && (a1.nrows() != n1 || a1.ncols() != n1) {
        return -5;
    }
    if n1 > 0 && (b1.nrows() != n1 || b1.ncols() != m1) {
        return -6;
    }
    if c1.nrows() != p1 || (n1 > 0 && c1.ncols() != n1) {
        return -7;
    }
    if d1.nrows() != p1 || d1.ncols() != m1 {
        return -8;
    }
    if n2 > 0 && (a2.nrows() != n2 || a2.ncols() != n2) {
        return -9;
    }
    if n2 > 0 && (b2.nrows() != n2 || b2.ncols() != m1) {
        return -10;
    }
    if c2.nrows() != p1 || (n2 > 0 && c2.ncols() != n2) {
        return -11;
    }
    if d2.nrows() != p1 || d2.ncols() != m1 {
        return -12;
    }
    if a.nrows() != n || a.ncols() != n {
        return -13;
    }
    if b.nrows() != n || b.ncols() != m1 {
        return -14;
    }
    if c.nrows() != p1 || c.ncols() != n {
        return -15;
    }
    if d.nrows() != p1 || d.ncols() != m1 {
        return -16;
    }
    if n == 0 || (m1 == 0 && p1 == 0) {
        return 0;
    }
    let alpha = 1.0;
    if n1 > 0 {
        a.view_mut((0, 0), (n1, n1)).copy_from(a1);
    }
    if n2 > 0 {
        a.view_mut((n1, n1), (n2, n2)).copy_from(a2);
        if n1 > 0 {
            a.view_mut((0, n1), (n1, n2)).fill(0.0);
            a.view_mut((n1, 0), (n2, n1)).fill(0.0);
        }
    }
    if n1 > 0 && m1 > 0 {
        b.rows_mut(0, n1).copy_from(b1);
    }
    if n2 > 0 && m1 > 0 {
        b.rows_mut(n1, n2).copy_from(b2);
    }
    if p1 > 0 && n1 > 0 {
        c.columns_mut(0, n1).copy_from(c1);
    }
    if p1 > 0 && n2 > 0 {
        c.columns_mut(n1, n2).copy_from(&(alpha * c2));
    }
    if p1 > 0 && m1 > 0 {
        *d = d1 + alpha * d2;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab05pd_trivial() {
        let a1 = DMatrix::zeros(0, 0);
        let b1 = DMatrix::zeros(0, 0);
        let c1 = DMatrix::zeros(0, 0);
        let d1 = DMatrix::zeros(0, 0);
        let a2 = DMatrix::zeros(0, 0);
        let b2 = DMatrix::zeros(0, 0);
        let c2 = DMatrix::zeros(0, 0);
        let d2 = DMatrix::zeros(0, 0);
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut c = DMatrix::zeros(0, 0);
        let mut d = DMatrix::zeros(0, 0);
        assert_eq!(
            ab05pd(
                0, 0, 0, 0,
                &a1, &b1, &c1, &d1,
                &a2, &b2, &c2, &d2,
                &mut a, &mut b, &mut c, &mut d
            ),
            0
        );
    }

    #[test]
    fn test_ab05pd_simple() {
        let n1 = 1;
        let n2 = 1;
        let m1 = 1;
        let p1 = 1;
        let a1 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let b1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d1 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let a2 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let b2 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c2 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d2 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut a = DMatrix::zeros(2, 2);
        let mut b = DMatrix::zeros(2, 1);
        let mut c = DMatrix::zeros(1, 2);
        let mut d = DMatrix::zeros(1, 1);
        let info = ab05pd(
            n1, m1, p1, n2,
            &a1, &b1, &c1, &d1,
            &a2, &b2, &c2, &d2,
            &mut a, &mut b, &mut c, &mut d,
        );
        assert_eq!(info, 0);
        assert!((a[(0, 0)] - 0.0).abs() < 1e-10);
        assert!((a[(1, 1)] - 0.0).abs() < 1e-10);
        assert!((d[(0, 0)] - 0.0).abs() < 1e-10);
    }
}
