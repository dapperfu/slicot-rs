//! ZLATZM — Apply complex Householder matrix (SLICOT/LAPACK auxiliary).
//!
//! P = I - tau*u*u^H, u = [1; v]. Overwrites C with P*C (Left) or C*P (Right).

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// Apply from left (P*C) or right (C*P).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ZlatzmSide {
    /// P * C (u has length m, C is m×n).
    Left,
    /// C * P (u has length n, C is m×n).
    Right,
}

/// Applies the Householder reflector P = I - tau*u*u^H to C (complex).
///
/// - **Left:** `c` is m×n, `v` has length m-1 (u = [1; v]), overwrites `c` with P*c.
/// - **Right:** `c` is m×n, `v` has length n-1 (u = [1; v]), overwrites `c` with c*P.
///
/// `incv` is the stride between elements of `v` (must be > 0). No-op if tau == 0 or min(m,n) == 0.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn zlatzm(
    side: ZlatzmSide,
    v: &[Complex64],
    incv: usize,
    tau: Complex64,
    c: &mut DMatrix<Complex64>,
) -> i32 {
    let m = c.nrows();
    let n = c.ncols();
    if incv == 0 {
        return -5;
    }
    if m == 0 || n == 0 || tau == Complex64::new(0.0, 0.0) {
        return 0;
    }
    match side {
        ZlatzmSide::Left => {
            if m > 1 && v.len() < (m - 1) * incv {
                return -4;
            }
            // w := u^H * C  (conjugate of u)
            let mut work = DVector::zeros(n);
            for j in 0..n {
                work[j] = c[(0, j)];
                for i in 1..m {
                    work[j] += v[(i - 1) * incv].conj() * c[(i, j)];
                }
            }
            for j in 0..n {
                c[(0, j)] -= tau * work[j];
            }
            for j in 0..n {
                for i in 1..m {
                    c[(i, j)] -= tau * v[(i - 1) * incv] * work[j];
                }
            }
        }
        ZlatzmSide::Right => {
            if n > 1 && v.len() < (n - 1) * incv {
                return -4;
            }
            // w := C * u
            let mut work = DVector::zeros(m);
            for i in 0..m {
                work[i] = c[(i, 0)];
                for j in 1..n {
                    work[i] += c[(i, j)] * v[(j - 1) * incv];
                }
            }
            for i in 0..m {
                c[(i, 0)] -= tau * work[i];
            }
            for i in 0..m {
                for j in 1..n {
                    c[(i, j)] -= tau * work[i] * v[(j - 1) * incv].conj();
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zlatzm_left_tau_zero_no_op() {
        let mut c = DMatrix::from_fn(2, 2, |i, j| Complex64::new((i * 2 + j) as f64 + 1.0, 0.0));
        let v = [Complex64::new(0.5, 0.0)];
        assert_eq!(
            zlatzm(ZlatzmSide::Left, &v, 1, Complex64::new(0.0, 0.0), &mut c),
            0
        );
        assert_eq!(c[(0, 0)].re, 1.0);
        assert_eq!(c[(1, 1)].re, 4.0);
    }

    #[test]
    fn test_zlatzm_right_tau_zero_no_op() {
        let mut c = DMatrix::from_fn(2, 2, |i, j| Complex64::new((i * 2 + j) as f64 + 1.0, 0.0));
        let v = [Complex64::new(0.5, 0.0)];
        assert_eq!(
            zlatzm(ZlatzmSide::Right, &v, 1, Complex64::new(0.0, 0.0), &mut c),
            0
        );
    }

    #[test]
    fn test_zlatzm_left_apply() {
        let mut c = DMatrix::from_fn(2, 2, |i, j| Complex64::new(1.0, 0.0));
        let v = [Complex64::new(1.0, 0.0)];
        let tau = Complex64::new(0.5, 0.0);
        assert_eq!(zlatzm(ZlatzmSide::Left, &v, 1, tau, &mut c), 0);
        // P = I - 0.5*u*u^H, u = [1,1]. u^H*u = 2, so diagonal change
        assert!(c[(0, 0)].re != 1.0 || c[(1, 1)].re != 1.0);
    }
}
