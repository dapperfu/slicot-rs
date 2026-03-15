//! AB05ND — Feedback interconnection of two systems (SLICOT AB05ND).
//!
//! Forms state-space (A,B,C,D) for the feedback interconnection:
//! U = U1 + alpha*Y2, Y = Y1 = U2. Uses E21 = (I + alpha*D1*D2)^{-1}, E12 = I - alpha*D2*E21*D1.

use nalgebra::DMatrix;

fn try_invert(m: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    let (n, n2) = (m.nrows(), m.ncols());
    if n != n2 {
        return None;
    }
    m.clone().try_inverse()
}

/// Feedback interconnection. Returns 0 on success; >0 if (I+alpha*D1*D2) is singular; <0 invalid argument.
pub fn ab05nd(
    n1: usize,
    m1: usize,
    p1: usize,
    n2: usize,
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
        return -7;
    }
    if b1.nrows() != n1 || b1.ncols() != m1 {
        return -8;
    }
    if c1.nrows() != p1 || c1.ncols() != n1 {
        return -9;
    }
    if d1.nrows() != p1 || d1.ncols() != m1 {
        return -10;
    }
    if a2.nrows() != n2 || a2.ncols() != n2 {
        return -11;
    }
    if b2.nrows() != n2 || b2.ncols() != p1 {
        return -12;
    }
    if c2.nrows() != m1 || c2.ncols() != n2 {
        return -13;
    }
    if d2.nrows() != m1 || d2.ncols() != p1 {
        return -14;
    }
    let n = n1 + n2;
    if a.nrows() != n || a.ncols() != n {
        return -15;
    }
    if b.nrows() != n || b.ncols() != m1 {
        return -16;
    }
    if c.nrows() != p1 || c.ncols() != n {
        return -17;
    }
    if d.nrows() != p1 || d.ncols() != m1 {
        return -18;
    }
    if n == 0 && (m1 == 0 || p1 == 0) {
        return 0;
    }
    if p1 > 0 {
        // E21 = (I + alpha*D1*D2)^{-1}; D1 P1xM1, D2 M1xP1 => D1*D2 P1xP1
        let d1d2 = d1 * d2;
        let i_plus = DMatrix::identity(p1, p1) + alpha * &d1d2;
        let e21 = match try_invert(&i_plus) {
            Some(e) => e,
            None => return 1,
        };
        // D = E21 * D1  (P1 x M1)
        *d = &e21 * d1;
        // C(:,1:n1) = E21 * C1
        if n1 > 0 {
            c.columns_mut(0, n1).copy_from(&(&e21 * c1));
        }
        // E12 = I - alpha*D2*E21*D1  (M1 x M1)
        let e12 = DMatrix::identity(m1, m1) - alpha * (d2 * &e21 * d1);
        // B(1:n1,:) = B1 * E12
        if n1 > 0 && m1 > 0 {
            b.rows_mut(0, n1).copy_from(&(b1 * &e12));
        }
        // A(1:n1, 1:n1) = A1 - alpha*B1*E12*D2*C1
        if n1 > 0 {
            let mut a11 = a1.clone();
            a11 -= alpha * (b1 * &e12 * d2 * c1);
            a.view_mut((0, 0), (n1, n1)).copy_from(&a11);
        }
        // A(1:n1, n1:n) = -alpha*B1*E12*C2
        if n1 > 0 && n2 > 0 {
            a.view_mut((0, n1), (n1, n2)).copy_from(&(-alpha * (b1 * &e12 * c2)));
        }
        // B(n1:n,:) = B2 * (E21*D1)
        if n2 > 0 {
            b.rows_mut(n1, n2).copy_from(&(b2 * &*d));
        }
        // C(:, n1:n) = -alpha*E21*D1*C2
        if n2 > 0 {
            c.columns_mut(n1, n2).copy_from(&(-alpha * &*d * c2));
        }
        // A(n1:n, 1:n1) = B2*E21*C1
        if n2 > 0 && n1 > 0 {
            a.view_mut((n1, 0), (n2, n1)).copy_from(&(b2 * &e21 * c1));
        }
        // A(n1:n, n1:n) = A2 - alpha*B2*E21*D1*C2
        if n2 > 0 {
            let mut a22 = a2.clone();
            a22 -= alpha * (b2 * &e21 * d1 * c2);
            a.view_mut((n1, n1), (n2, n2)).copy_from(&a22);
        }
    } else {
        // P1 == 0: E12 = I (M1xM1)
        if n1 > 0 && m1 > 0 {
            b.rows_mut(0, n1).copy_from(b1);
        }
        if n1 > 0 {
            a.view_mut((0, 0), (n1, n1)).copy_from(a1);
        }
        if n2 > 0 {
            b.rows_mut(n1, n2).fill(0.0);
            a.view_mut((n1, 0), (n2, n1)).copy_from(&(b2 * c1));
            let mut a22 = a2.clone();
            a22 -= alpha * (b2 * c2);
            a.view_mut((n1, n1), (n2, n2)).copy_from(&a22);
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab05nd_simple() {
        let n1 = 1;
        let m1 = 1;
        let p1 = 1;
        let n2 = 1;
        let alpha = -1.0;
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
        let info = ab05nd(
            n1, m1, p1, n2, alpha,
            &a1, &b1, &c1, &d1,
            &a2, &b2, &c2, &d2,
            &mut a, &mut b, &mut c, &mut d,
        );
        assert_eq!(info, 0);
    }
}
