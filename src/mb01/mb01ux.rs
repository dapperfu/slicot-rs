//! MB01UX — A := alpha*op(T)*A or A := alpha*A*op(T) (SLICOT MB01UX)
// T quasi-triangular (we implement upper triangular).

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UxSide {
    Left,  // A := alpha*op(T)*A, T is M×M
    Right, // A := alpha*A*op(T), T is N×N
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UxUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UxTrans {
    NoTrans,
    Trans,
}

/// Overwrites A with alpha*op(T)*A (Left) or alpha*A*op(T) (Right).
pub fn mb01ux(
    side: Mb01UxSide,
    uplo: Mb01UxUplo,
    trans: Mb01UxTrans,
    m: usize,
    n: usize,
    alpha: f64,
    t: &[f64],
    ldt: usize,
    a: &mut [f64],
    lda: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(side, Mb01UxSide::Left | Mb01UxSide::Right) {
        return -1;
    }
    if !matches!(uplo, Mb01UxUplo::Upper | Mb01UxUplo::Lower) {
        return -2;
    }
    if !matches!(trans, Mb01UxTrans::NoTrans | Mb01UxTrans::Trans) {
        return -3;
    }
    let k = if side == Mb01UxSide::Left { m } else { n };
    if ldt < k.max(1) || lda < m.max(1) {
        return -8;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    if alpha == 0.0 {
        for j in 0..n {
            for i in 0..m {
                a[i + j * lda] = 0.0;
            }
        }
        return 0;
    }

    let t_mat = DMatrix::from_fn(k, k, |i, j| {
        if (uplo == Mb01UxUplo::Upper && i <= j) || (uplo == Mb01UxUplo::Lower && i >= j) {
            t[i + j * ldt]
        } else {
            0.0
        }
    });
    let a_mat = DMatrix::from_fn(m, n, |i, j| a[i + j * lda]);
    let op_t = match trans {
        Mb01UxTrans::NoTrans => t_mat.clone(),
        Mb01UxTrans::Trans => t_mat.transpose(),
    };
    let result = match side {
        Mb01UxSide::Left => alpha * &op_t * &a_mat,
        Mb01UxSide::Right => alpha * &a_mat * &op_t,
    };
    for j in 0..n {
        for i in 0..m {
            a[i + j * lda] = result[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01ux_left_upper() {
        let m = 2;
        let n = 2;
        let t = [1.0, 0.0, 0.0, 1.0];
        let mut a = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; m * n];
        assert_eq!(
            mb01ux(
                Mb01UxSide::Left,
                Mb01UxUplo::Upper,
                Mb01UxTrans::NoTrans,
                m,
                n,
                1.0,
                &t,
                2,
                &mut a,
                2,
                &mut dwork,
            ),
            0
        );
        assert!((a[0] - 1.0).abs() < 1e-14);
        assert!((a[3] - 1.0).abs() < 1e-14);
    }
}
