//! TF01ND — Output sequence with A upper or lower Hessenberg (SLICOT TF01ND)
//!
//! Same recurrence as TF01MD; A is Hessenberg for efficiency.

use nalgebra::{DMatrix, DVector};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// U is M×NY (column k = u(k)); Y is P×NY (column k = y(k)).
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tf01nd(
    uplo: Uplo,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    u: &DMatrix<f64>,
    x: &mut DVector<f64>,
    y: &mut DMatrix<f64>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    let ny = u.ncols();
    if a.ncols() != n {
        return -5;
    }
    if b.nrows() != n || c.ncols() != n || d.nrows() != p || d.ncols() != m {
        return -6;
    }
    if u.nrows() != m {
        return -9;
    }
    if x.len() != n {
        return -11;
    }
    if y.nrows() != p || y.ncols() != ny {
        return -12;
    }
    if n == 0 || ny == 0 {
        return 0;
    }

    let mut xcur = x.clone();
    for k in 0..ny {
        let uk = u.column(k);
        let yk = c * &xcur + d * &uk;
        y.column_mut(k).copy_from(&yk);
        xcur = a * &xcur + b * &uk;
    }
    x.copy_from(&xcur);
    let _ = uplo;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tf01nd_smoke() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 0.0, 0.5]);
        let b = DMatrix::from_row_slice(2, 1, &[0.0, 1.0]);
        let c = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let u = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let mut x = DVector::from_row_slice(&[0.0, 0.0]);
        let mut y = DMatrix::zeros(1, 2);
        assert_eq!(tf01nd(Uplo::Upper, &a, &b, &c, &d, &u, &mut x, &mut y), 0);
        assert!(y[(0, 0)].is_finite());
    }
}
