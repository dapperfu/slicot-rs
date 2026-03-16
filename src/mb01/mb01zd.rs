//! MB01ZD — H := alpha*op(T)*H or H := alpha*H*op(T), H Hessenberg-like, T triangular (SLICOT MB01ZD)

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01ZdSide {
    Left,
    Right,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01ZdUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01ZdTrans {
    NoTrans,
    Trans,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01ZdDiag {
    Unit,
    NonUnit,
}

/// Computes H := alpha*op(T)*H (SIDE='L') or H := alpha*H*op(T) (SIDE='R').
/// T is triangular (unit or non-unit), H is Hessenberg-like with L nonzero sub/superdiagonals.
pub fn mb01zd(
    side: Mb01ZdSide,
    uplo: Mb01ZdUplo,
    trans: Mb01ZdTrans,
    diag: Mb01ZdDiag,
    m: usize,
    n: usize,
    _l: usize,
    alpha: f64,
    t: &DMatrix<f64>,
    h: &mut DMatrix<f64>,
) -> i32 {
    let k = match side {
        Mb01ZdSide::Left => m,
        Mb01ZdSide::Right => n,
    };
    if t.nrows() != k || t.ncols() != k {
        return -10;
    }
    if h.nrows() != m || h.ncols() != n {
        return -12;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    if alpha == 0.0 {
        for i in 0..m {
            for j in 0..n {
                h[(i, j)] = 0.0;
            }
        }
        return 0;
    }

    let mut t_full = DMatrix::zeros(k, k);
    for i in 0..k {
        for j in 0..k {
            t_full[(i, j)] = if (uplo == Mb01ZdUplo::Upper && i <= j)
                || (uplo == Mb01ZdUplo::Lower && i >= j)
            {
                if diag == Mb01ZdDiag::Unit && i == j {
                    1.0
                } else {
                    t[(i, j)]
                }
            } else {
                0.0
            };
        }
    }

    let op_t = match trans {
        Mb01ZdTrans::NoTrans => t_full.clone(),
        Mb01ZdTrans::Trans => t_full.transpose(),
    };

    let result = match side {
        Mb01ZdSide::Left => alpha * &op_t * h.clone(),
        Mb01ZdSide::Right => alpha * h.clone() * &op_t,
    };

    for i in 0..m {
        for j in 0..n {
            h[(i, j)] = result[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01zd_left_upper_notrans() {
        // H := 1*T*H, T 2x2 upper [1 1; 0 1], H 2x2 [1 0; 1 1]
        let t = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 0.0, 1.0]);
        let mut h = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 1.0, 1.0]);
        assert_eq!(
            mb01zd(
                Mb01ZdSide::Left,
                Mb01ZdUplo::Upper,
                Mb01ZdTrans::NoTrans,
                Mb01ZdDiag::NonUnit,
                2,
                2,
                1,
                1.0,
                &t,
                &mut h,
            ),
            0
        );
        // T*H = [1 1; 0 1]*[1 0; 1 1] = [2 1; 1 1]
        assert!((h[(0, 0)] - 2.0).abs() < 1e-10);
        assert!((h[(0, 1)] - 1.0).abs() < 1e-10);
        assert!((h[(1, 0)] - 1.0).abs() < 1e-10);
        assert!((h[(1, 1)] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mb01zd_trivial() {
        let t = DMatrix::<f64>::zeros(0, 0);
        let mut h = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(
            mb01zd(
                Mb01ZdSide::Left,
                Mb01ZdUplo::Upper,
                Mb01ZdTrans::NoTrans,
                Mb01ZdDiag::NonUnit,
                0,
                0,
                0,
                1.0,
                &t,
                &mut h,
            ),
            0
        );
    }
}
