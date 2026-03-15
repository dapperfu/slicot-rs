//! AB05MD — Cascade interconnection of two state-space systems (SLICOT AB05MD)
//!
//! Forms (A,B,C,D) for the cascaded system from (A1,B1,C1,D1) and (A2,B2,C2,D2).

use nalgebra::DMatrix;

/// Block structure of the combined state matrix A.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    /// Lower block diagonal: A = [A1 0; B2*C1 A2].
    Lower,
    /// Upper block diagonal: A = [A2 B2*C1; 0 A1].
    Upper,
}

/// Forms the state-space model (A,B,C,D) for the cascade of two systems.
///
/// System 1: x1' = A1*x1 + B1*u, v = C1*x1 + D1*u.
/// System 2: x2' = A2*x2 + B2*v, y = C2*x2 + D2*v.
/// Result: x' = A*x + B*u, y = C*x + D*u.
///
/// # Arguments
/// * `uplo` - Lower or Upper block form of A.
/// * `a1,b1,c1,d1` - First system (N1×N1, N1×M1, P1×N1, P1×M1); P1 = inputs to second system.
/// * `a2,b2,c2,d2` - Second system (N2×N2, N2×P1, P2×N2, P2×P1).
/// * `a,b,c,d` - Output matrices (must be sized (N1+N2)×(N1+N2), (N1+N2)×M1, P2×(N1+N2), P2×M1).
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn ab05md(
    uplo: Uplo,
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
    let n1 = a1.nrows();
    let m1 = b1.ncols();
    let p1 = c1.nrows();
    let n2 = a2.nrows();
    let p2 = c2.nrows();
    if a1.ncols() != n1 || b1.nrows() != n1 || c1.ncols() != n1 || d1.nrows() != p1 || d1.ncols() != m1 {
        return -9;
    }
    if a2.ncols() != n2 || b2.nrows() != n2 || b2.ncols() != p1 || c2.ncols() != n2 || d2.nrows() != p2 || d2.ncols() != p1 {
        return -17;
    }
    let n = n1 + n2;
    if a.nrows() != n || a.ncols() != n || b.nrows() != n || b.ncols() != m1 || c.nrows() != p2 || c.ncols() != n || d.nrows() != p2 || d.ncols() != m1 {
        return -26;
    }
    if n == 0 && m1 == 0 && p2 == 0 {
        return 0;
    }

    match uplo {
        Uplo::Lower => {
            a.view_mut((0, 0), (n1, n1)).copy_from(a1);
            if n2 > 0 {
                a.view_mut((n1, n1), (n2, n2)).copy_from(a2);
            }
            a.view_mut((0, n1), (n1, n2)).fill(0.0);
            if n1 > 0 && n2 > 0 && p1 > 0 {
                a.view_mut((n1, 0), (n2, n1)).copy_from(&(b2 * c1));
            }
            if n1 > 0 && m1 > 0 {
                b.view_mut((0, 0), (n1, m1)).copy_from(b1);
            }
            if n2 > 0 && m1 > 0 && p1 > 0 {
                b.view_mut((n1, 0), (n2, m1)).copy_from(&(b2 * d1));
            }
            if n1 > 0 && p1 > 0 && p2 > 0 {
                c.view_mut((0, 0), (p2, n1)).copy_from(&(d2 * c1));
            }
            if n2 > 0 && p2 > 0 {
                c.view_mut((0, n1), (p2, n2)).copy_from(c2);
            }
            if p2 > 0 && m1 > 0 && p1 > 0 {
                d.copy_from(&(d2 * d1));
            }
        }
        Uplo::Upper => {
            if n2 > 0 {
                a.view_mut((0, 0), (n2, n2)).copy_from(a2);
            }
            if n1 > 0 {
                a.view_mut((n2, n2), (n1, n1)).copy_from(a1);
            }
            a.view_mut((n2, 0), (n1, n2)).fill(0.0);
            if n1 > 0 && n2 > 0 && p1 > 0 {
                a.view_mut((0, n2), (n2, n1)).copy_from(&(b2 * c1));
            }
            if n2 > 0 && m1 > 0 && p1 > 0 {
                b.view_mut((0, 0), (n2, m1)).copy_from(&(b2 * d1));
            }
            if n1 > 0 && m1 > 0 {
                b.view_mut((n2, 0), (n1, m1)).copy_from(b1);
            }
            if n2 > 0 && p2 > 0 {
                c.view_mut((0, 0), (p2, n2)).copy_from(c2);
            }
            if n1 > 0 && p1 > 0 && p2 > 0 {
                c.view_mut((0, n2), (p2, n1)).copy_from(&(d2 * c1));
            }
            if p2 > 0 && m1 > 0 && p1 > 0 {
                d.copy_from(&(d2 * d1));
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab05md_lower() {
        let a1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let b1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c1 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d1 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let a2 = DMatrix::from_row_slice(1, 1, &[2.0]);
        let b2 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c2 = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d2 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut a = DMatrix::zeros(2, 2);
        let mut b = DMatrix::zeros(2, 1);
        let mut c = DMatrix::zeros(1, 2);
        let mut d = DMatrix::zeros(1, 1);
        let info = ab05md(Uplo::Lower, &a1, &b1, &c1, &d1, &a2, &b2, &c2, &d2, &mut a, &mut b, &mut c, &mut d);
        assert_eq!(info, 0);
        assert!((a[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((a[(1, 1)] - 2.0).abs() < 1e-10);
        assert!((a[(1, 0)] - 1.0).abs() < 1e-10); // B2*C1
    }
}
