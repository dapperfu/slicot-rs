//! TF01MY — Output sequence (A,B,C,D) with row-wise U and Y (SLICOT TF01MY)
//!
//! Same recurrence as TF01MD but U is NY×M (row k = u(k)'), Y is NY×P (row k = y(k)').

use nalgebra::{DMatrix, DVector};

/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tf01my(
    n: usize,
    m: usize,
    p: usize,
    ny: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    u: &DMatrix<f64>,
    x: &mut DVector<f64>,
    y: &mut DMatrix<f64>,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -5;
    }
    if b.nrows() != n || b.ncols() != m || c.nrows() != p || c.ncols() != n || d.nrows() != p || d.ncols() != m {
        return -6;
    }
    if u.nrows() != ny || u.ncols() != m {
        return -9;
    }
    if x.len() != n {
        return -11;
    }
    if y.nrows() != ny || y.ncols() != p {
        return -12;
    }
    if n == 0 || ny == 0 {
        return 0;
    }

    let mut xcur = x.clone();
    for k in 0..ny {
        let uk = u.row(k).transpose();
        let yk = c * &xcur + d * &uk;
        for i in 0..p {
            y[(k, i)] = yk[i];
        }
        xcur = a * &xcur + b * &uk;
    }
    x.copy_from(&xcur);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tf01my_simple() {
        let n = 2;
        let m = 1;
        let p = 1;
        let ny = 3;
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let c = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let u = DMatrix::from_row_slice(3, 1, &[1.0, 0.0, 0.0]);
        let mut x = DVector::from_row_slice(&[0.0, 0.0]);
        let mut y = DMatrix::zeros(ny, p);
        assert_eq!(tf01my(n, m, p, ny, &a, &b, &c, &d, &u, &mut x, &mut y), 0);
        assert!((y[(0, 0)] - 0.0).abs() < 1e-10);
        assert!((y[(1, 0)] - 1.0).abs() < 1e-10);
    }
}
