//! AB05QD — Append two systems in state-space form (SLICOT AB05QD).
//!
//! Builds G = diag(G1, G2): A = diag(A1,A2), B = [B1 0; 0 B2],
//! C = [C1 0; 0 C2], D = [D1 0; 0 D2]. Output dimensions: n = n1+n2, m = m1+m2, p = p1+p2.

use nalgebra::DMatrix;

/// Append two systems. Returns 0 on success; < 0 invalid argument index.
pub fn ab05qd(
    n1: usize,
    m1: usize,
    p1: usize,
    n2: usize,
    m2: usize,
    p2: usize,
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
    let m = m1 + m2;
    let p = p1 + p2;
    if n1 > 0 && (a1.nrows() != n1 || a1.ncols() != n1) {
        return -8;
    }
    if n1 > 0 && (b1.nrows() != n1 || b1.ncols() != m1) {
        return -9;
    }
    if p1 > 0 && (c1.nrows() != p1 || (n1 > 0 && c1.ncols() != n1)) {
        return -10;
    }
    if p1 > 0 && d1.ncols() != m1 {
        return -11;
    }
    if d1.nrows() != p1 {
        return -11;
    }
    if n2 > 0 && (a2.nrows() != n2 || a2.ncols() != n2) {
        return -13;
    }
    if n2 > 0 && (b2.nrows() != n2 || b2.ncols() != m2) {
        return -14;
    }
    if p2 > 0 && (c2.nrows() != p2 || (n2 > 0 && c2.ncols() != n2)) {
        return -15;
    }
    if d2.nrows() != p2 || d2.ncols() != m2 {
        return -16;
    }
    if a.nrows() != n || a.ncols() != n {
        return -18;
    }
    if b.nrows() != n || b.ncols() != m {
        return -20;
    }
    if c.nrows() != p || c.ncols() != n {
        return -22;
    }
    if d.nrows() != p || d.ncols() != m {
        return -24;
    }
    if n == 0 && (m == 0 || p == 0) {
        return 0;
    }
    a.fill(0.0);
    b.fill(0.0);
    c.fill(0.0);
    d.fill(0.0);
    if n1 > 0 {
        a.view_mut((0, 0), (n1, n1)).copy_from(a1);
        if m1 > 0 {
            b.view_mut((0, 0), (n1, m1)).copy_from(b1);
        }
        if p1 > 0 {
            c.view_mut((0, 0), (p1, n1)).copy_from(c1);
            if m1 > 0 {
                d.view_mut((0, 0), (p1, m1)).copy_from(d1);
            }
        }
    }
    if n2 > 0 {
        a.view_mut((n1, n1), (n2, n2)).copy_from(a2);
        if m2 > 0 {
            b.view_mut((n1, m1), (n2, m2)).copy_from(b2);
        }
        if p2 > 0 {
            c.view_mut((p1, n1), (p2, n2)).copy_from(c2);
            if m2 > 0 {
                d.view_mut((p1, m1), (p2, m2)).copy_from(d2);
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab05qd_trivial() {
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
            ab05qd(
                0, 0, 0, 0, 0, 0,
                &a1, &b1, &c1, &d1,
                &a2, &b2, &c2, &d2,
                &mut a, &mut b, &mut c, &mut d
            ),
            0
        );
    }

    #[test]
    fn test_ab05qd_append() {
        let n1 = 1;
        let m1 = 1;
        let p1 = 1;
        let n2 = 1;
        let m2 = 1;
        let p2 = 1;
        let a1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let b1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d1 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let a2 = DMatrix::from_row_slice(1, 1, &[2.0]);
        let b2 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c2 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d2 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut a = DMatrix::zeros(2, 2);
        let mut b = DMatrix::zeros(2, 2);
        let mut c = DMatrix::zeros(2, 2);
        let mut d = DMatrix::zeros(2, 2);
        let info = ab05qd(
            n1, m1, p1, n2, m2, p2,
            &a1, &b1, &c1, &d1,
            &a2, &b2, &c2, &d2,
            &mut a, &mut b, &mut c, &mut d,
        );
        assert_eq!(info, 0);
        assert!((a[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((a[(1, 1)] - 2.0).abs() < 1e-10);
    }
}
