//! DLATZM — Apply Householder matrix from DTZRQF (SLICOT/LAPACK auxiliary)
//
// P = I - tau*u*u', u = [1; v]. Overwrites C with P*C (Left) or C*P (Right).

use nalgebra::{DMatrix, DVector};

/// Apply from left (P*C) or right (C*P).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum DlatzmSide {
    /// P * C (u has length m, C is m×n).
    Left,
    /// C * P (u has length n, C is m×n).
    Right,
}

/// Applies the Householder reflector P = I - tau*u*u' to C.
///
/// - **Left:** `c` is m×n, `v` has length m-1 (u = [1; v]), overwrites `c` with P*c.
/// - **Right:** `c` is m×n, `v` has length n-1 (u = [1; v]), overwrites `c` with c*P.
///
/// `incv` is the stride between elements of `v` (must be > 0). No-op if tau == 0 or min(m,n) == 0.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn dlatzm(
    side: DlatzmSide,
    v: &[f64],
    incv: usize,
    tau: f64,
    c: &mut DMatrix<f64>,
) -> i32 {
    let m = c.nrows();
    let n = c.ncols();
    if incv == 0 {
        return -5;
    }
    if m == 0 || n == 0 || tau == 0.0 {
        return 0;
    }
    match side {
        DlatzmSide::Left => {
            // V dimension in Fortran: 1 + (M-1)*abs(INCV); when m > 1 need v long enough
            if m > 1 && v.len() < 1 + (m - 1) * incv {
                return -4;
            }
            // w := (C1 + v'*C2)'  (C1 = first row, C2 = rows 1..m-1)
            let mut work = DVector::zeros(n);
            for j in 0..n {
                work[j] = c[(0, j)];
            }
            for j in 0..n {
                for i in 1..m {
                    work[j] += v[(i - 1) * incv] * c[(i, j)];
                }
            }
            // [C1; C2] := [C1; C2] - tau * u * w'
            for j in 0..n {
                c[(0, j)] -= tau * work[j];
            }
            for j in 0..n {
                for i in 1..m {
                    c[(i, j)] -= tau * v[(i - 1) * incv] * work[j];
                }
            }
        }
        DlatzmSide::Right => {
            if n > 1 && v.len() < 1 + (n - 1) * incv {
                return -4;
            }
            // w := C1 + C2*v  (C1 = first column, C2 = columns 1..n-1)
            let mut work = DVector::zeros(m);
            for i in 0..m {
                work[i] = c[(i, 0)];
            }
            for i in 0..m {
                for j in 1..n {
                    work[i] += c[(i, j)] * v[(j - 1) * incv];
                }
            }
            // [C1, C2] := [C1, C2] - tau * w * [1, v']
            for i in 0..m {
                c[(i, 0)] -= tau * work[i];
            }
            for i in 0..m {
                for j in 1..n {
                    c[(i, j)] -= tau * work[i] * v[(j - 1) * incv];
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
    fn test_dlatzm_left_tau_zero_no_op() {
        let mut c = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let v = [0.5];
        assert_eq!(dlatzm(DlatzmSide::Left, &v, 1, 0.0, &mut c), 0);
        assert_eq!(c[(0, 0)], 1.0);
        assert_eq!(c[(1, 1)], 4.0);
    }

    #[test]
    fn test_dlatzm_right_tau_zero_no_op() {
        let mut c = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let v = [0.5];
        assert_eq!(dlatzm(DlatzmSide::Right, &v, 1, 0.0, &mut c), 0);
        assert_eq!(c[(0, 0)], 1.0);
    }
}
