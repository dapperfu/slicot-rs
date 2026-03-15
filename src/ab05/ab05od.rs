//! AB05OD — Parallel interconnection (rowwise concatenation) of two systems (SLICOT AB05OD).
//!
//! Forms (A,B,C,D) with A = diag(A1,A2), B = [B1 0; 0 B2], C = [C1; alpha*C2], D = [D1 alpha*D2].

use nalgebra::DMatrix;

/// Parallel connection. Returns 0 on success; <0 invalid argument.
pub fn ab05od(
    n1: usize,
    m1: usize,
    p1: usize,
    n2: usize,
    m2: usize,
    alpha: f64,
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
    if a1.nrows() != n1 || a1.ncols() != n1 {
        return -8;
    }
    if b1.nrows() != n1 || b1.ncols() != m1 {
        return -10;
    }
    if c1.nrows() != p1 || c1.ncols() != n1 {
        return -12;
    }
    if d1.nrows() != p1 || d1.ncols() != m1 {
        return -14;
    }
    if a2.nrows() != n2 || a2.ncols() != n2 {
        return -16;
    }
    if b2.nrows() != n2 || b2.ncols() != m2 {
        return -18;
    }
    if c2.nrows() != p1 || c2.ncols() != n2 {
        return -20;
    }
    if d2.nrows() != p1 || d2.ncols() != m2 {
        return -22;
    }
    let n = n1 + n2;
    let m = m1 + m2;
    if a.nrows() != n || a.ncols() != n {
        return -25;
    }
    if b.nrows() != n || b.ncols() != m {
        return -27;
    }
    if c.nrows() != p1 || c.ncols() != n {
        return -29;
    }
    if d.nrows() != p1 || d.ncols() != m {
        return -31;
    }
    if n == 0 && (m == 0 || p1 == 0) {
        return 0;
    }
    a.fill(0.0);
    b.fill(0.0);
    if n1 > 0 {
        a.view_mut((0, 0), (n1, n1)).copy_from(a1);
        if m1 > 0 {
            b.view_mut((0, 0), (n1, m1)).copy_from(b1);
        }
    }
    if n2 > 0 {
        a.view_mut((n1, n1), (n2, n2)).copy_from(a2);
        if m2 > 0 {
            b.view_mut((n1, m1), (n2, m2)).copy_from(b2);
        }
    }
    if n1 > 0 {
        c.view_mut((0, 0), (p1, n1)).copy_from(c1);
    }
    if n2 > 0 {
        let alpha_c2 = alpha * c2;
        c.view_mut((0, n1), (p1, n2)).copy_from(&alpha_c2);
    }
    if m1 > 0 {
        d.view_mut((0, 0), (p1, m1)).copy_from(d1);
    }
    if m2 > 0 {
        let alpha_d2 = alpha * d2;
        d.view_mut((0, m1), (p1, m2)).copy_from(&alpha_d2);
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab05od_simple() {
        let a1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let b1 = DMatrix::from_row_slice(1, 1, &[0.5]);
        let c1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d1 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let a2 = DMatrix::from_row_slice(1, 1, &[2.0]);
        let b2 = DMatrix::from_row_slice(1, 1, &[0.5]);
        let c2 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d2 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut a = DMatrix::zeros(2, 2);
        let mut b = DMatrix::zeros(2, 2);
        let mut c = DMatrix::zeros(1, 2);
        let mut d = DMatrix::zeros(1, 2);
        let info = ab05od(
            1, 1, 1, 1, 1, 1.0,
            &a1, &b1, &c1, &d1,
            &a2, &b2, &c2, &d2,
            &mut a, &mut b, &mut c, &mut d,
        );
        assert_eq!(info, 0);
        assert_eq!(a[(0, 0)], 1.0);
        assert_eq!(a[(1, 1)], 2.0);
    }
}
