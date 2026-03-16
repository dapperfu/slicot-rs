//! AB07ND — Compute the inverse (Ai, Bi, Ci, Di) of a given system (A, B, C, D) (SLICOT AB07ND).
//!
//! Ai = A - B*D^{-1}*C, Bi = -B*D^{-1}, Ci = D^{-1}*C, Di = D^{-1}.

use nalgebra::DMatrix;

/// Compute inverse system. A, B, C, D are overwritten by Ai, Bi, Ci, Di.
///
/// # Returns
/// * 0: success
/// * < 0: invalid argument (-i)
/// * i (1..=M): D is singular, (i,i) zero; RCOND set to 0
/// * M+1: D numerically singular (RCOND < eps)
pub fn ab07nd(
    n: usize,
    m: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    rcond: &mut f64,
) -> i32 {
    if n > 0 && (a.nrows() != n || a.ncols() != n) {
        return -4;
    }
    if n > 0 && (b.nrows() != n || b.ncols() != m) {
        return -6;
    }
    if c.nrows() != m || (n > 0 && c.ncols() != n) {
        return -8;
    }
    if d.nrows() != m || d.ncols() != m {
        return -10;
    }

    if m == 0 {
        *rcond = 1.0;
        return 0;
    }

    // Factorize D (LU) and invert
    let d_norm = d.norm();
    let d_copy = d.clone();
    let lu = d_copy.lu();
    if let Some(d_inv) = lu.try_inverse() {
        // Reciprocal condition number (1-norm): RCOND = 1/(||D||_1 * ||D^{-1}||_1)
        let dinv_norm = d_inv.norm();
        *rcond = if d_norm > 0.0 && dinv_norm > 0.0 {
            1.0 / (d_norm * dinv_norm)
        } else {
            0.0
        };
        let eps = f64::EPSILON;
        if *rcond < eps {
            return (m + 1) as i32;
        }

        *d = d_inv.clone();

        if n > 0 {
            // Bi = -B * D^{-1}
            let bi = -(&*b) * &d_inv;
            b.copy_from(&bi);
            // Ai = A + Bi * C = A - B*D^{-1}*C
            let bc = &*b * &*c;
            *a += &bc;
            // Ci = D^{-1} * C
            let ci = &d_inv * &*c;
            c.copy_from(&ci);
        }
        0
    } else {
        *rcond = 0.0;
        (m + 1) as i32
    }
}

/// Convenience wrapper for benchmarking: (n, m) -> INFO. Uses zero matrices.
#[inline]
pub fn ab07nd_nm(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    if m == 0 {
        return 0;
    }
    let mut a = DMatrix::zeros(n, n);
    let mut b = DMatrix::zeros(n, m);
    let mut c = DMatrix::zeros(m, n);
    let mut d = DMatrix::identity(m, m); // non-singular
    let mut rcond = 0.0;
    ab07nd(n, m, &mut a, &mut b, &mut c, &mut d, &mut rcond)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab07nd_trivial() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut c = DMatrix::zeros(0, 0);
        let mut d = DMatrix::zeros(0, 0);
        let mut rcond = 0.0;
        assert_eq!(ab07nd(0, 0, &mut a, &mut b, &mut c, &mut d, &mut rcond), 0);
        assert_eq!(rcond, 1.0);
    }

    #[test]
    fn test_ab07nd_m0() {
        let mut a = DMatrix::zeros(2, 2);
        let mut b = DMatrix::zeros(2, 0);
        let mut c = DMatrix::zeros(0, 2);
        let mut d = DMatrix::zeros(0, 0);
        let mut rcond = 0.0;
        assert_eq!(ab07nd(2, 0, &mut a, &mut b, &mut c, &mut d, &mut rcond), 0);
        assert_eq!(rcond, 1.0);
    }

    #[test]
    fn test_ab07nd_inverse_1x1() {
        // D=1 -> Di=1, Bi=-B, Ai=A-B*C, Ci=C
        let mut a = DMatrix::from_row_slice(1, 1, &[2.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut d = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut rcond = 0.0;
        assert_eq!(ab07nd(1, 1, &mut a, &mut b, &mut c, &mut d, &mut rcond), 0);
        assert!((d[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((b[(0, 0)] - (-1.0)).abs() < 1e-10);
        assert!((c[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((a[(0, 0)] - (2.0 - 1.0 * 1.0)).abs() < 1e-10);
    }
}
