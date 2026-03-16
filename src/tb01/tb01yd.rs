//! TB01YD — Special similarity transformation P*A*P, P*B, C*P (SLICOT TB01YD)
//!
//! Applies transformation with P = anti-diagonal identity (1 on secondary diagonal).

use nalgebra::DMatrix;

/// Applies P*A*P, P*B, C*P in place. P has 1 on the secondary diagonal (i+j = n-1).
///
/// # Returns
/// 0 on success; &lt; 0 if the i-th argument had an illegal value.
pub fn tb01yd(
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n {
        return -5;
    }
    if b.nrows() != n {
        return -7;
    }
    if c.ncols() != n {
        return -9;
    }
    if n == 0 {
        return 0;
    }
    // A <- P*A*P: (P*A*P)[i,j] = A[n-1-j, n-1-i]
    let ac = a.clone();
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = ac[(n - 1 - j, n - 1 - i)];
        }
    }
    // B <- P*B: (P*B)[i,k] = B[n-1-i, k]
    let bc = b.clone();
    for i in 0..n {
        for k in 0..m {
            b[(i, k)] = bc[(n - 1 - i, k)];
        }
    }
    // C <- C*P: (C*P)[k,j] = C[k, n-1-j]
    let cc = c.clone();
    for k in 0..p {
        for j in 0..n {
            c[(k, j)] = cc[(k, n - 1 - j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb01yd_identity() {
        let n = 2;
        let m = 1;
        let p = 1;
        let mut a = DMatrix::identity(n, n);
        let mut b = DMatrix::from_fn(n, m, |i, _| if i == 0 { 1.0 } else { 0.0 });
        let mut c = DMatrix::from_fn(p, n, |_, j| if j == 0 { 1.0 } else { 0.0 });
        let info = tb01yd(&mut a, &mut b, &mut c);
        assert_eq!(info, 0);
        // P*I*P = I (identity, since P*P = I for anti-diagonal P)
        assert!((a[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((a[(0, 1)] - 0.0).abs() < 1e-10);
        assert!((a[(1, 0)] - 0.0).abs() < 1e-10);
        assert!((a[(1, 1)] - 1.0).abs() < 1e-10);
    }
}
