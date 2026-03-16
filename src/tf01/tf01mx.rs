//! TF01MX — Output sequence from system matrix S (SLICOT TF01MX)
//!
//! Computes y(1),...,y(NY) from (x(k+1); y(k)) = S * (x(k); u(k)) with S = [A B; C D].

use nalgebra::DMatrix;
use nalgebra::DVector;

/// U is NY×M (row k = u(k)'); Y is NY×P (row k = y(k)').
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tf01mx(
    n: usize,
    m: usize,
    p: usize,
    ny: usize,
    s: &DMatrix<f64>,
    u: &DMatrix<f64>,
    x: &mut DVector<f64>,
    y: &mut DMatrix<f64>,
) -> i32 {
    if s.nrows() != n + p || s.ncols() != n + m {
        return -6;
    }
    if u.nrows() != ny || u.ncols() != m {
        return -9;
    }
    if x.len() != n {
        return -10;
    }
    if y.nrows() != ny || y.ncols() != p {
        return -12;
    }
    if n == 0 && (ny == 0 || p == 0) {
        return 0;
    }

    let mut xcur = x.clone();
    for k in 0..ny {
        let uk = u.row(k).transpose();
        let mut xu = DVector::zeros(n + m);
        for i in 0..n {
            xu[i] = xcur[i];
        }
        for i in 0..m {
            xu[n + i] = uk[i];
        }
        let out = s * &xu;
        for i in 0..p {
            y[(k, i)] = out[n + i];
        }
        for i in 0..n {
            xcur[i] = out[i];
        }
    }
    x.copy_from(&xcur);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tf01mx_simple() {
        let n = 2;
        let m = 1;
        let p = 1;
        let ny = 3;
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let c = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut s = DMatrix::zeros(n + p, n + m);
        s.view_mut((0, 0), (n, n)).copy_from(&a);
        s.view_mut((0, n), (n, m)).copy_from(&b);
        s.view_mut((n, 0), (p, n)).copy_from(&c);
        s.view_mut((n, n), (p, m)).copy_from(&d);
        let u = DMatrix::from_row_slice(3, 1, &[1.0, 0.0, 0.0]);
        let mut x = DVector::from_row_slice(&[0.0, 0.0]);
        let mut y = DMatrix::zeros(ny, p);
        assert_eq!(tf01mx(n, m, p, ny, &s, &u, &mut x, &mut y), 0);
        assert!((y[(0, 0)] - 0.0).abs() < 1e-10);
        assert!((y[(1, 0)] - 1.0).abs() < 1e-10);
        assert!((y[(2, 0)] - 1.0).abs() < 1e-10);
    }
}
