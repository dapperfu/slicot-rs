//! MB01UY — T := alpha*op(T)*A or T := alpha*A*op(T) (SLICOT MB01UY)
// Result (M×N) overwrites the leading part of T.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UySide {
    Left,  // T := alpha*op(T)*A, T input M×M, A M×N, result M×N
    Right, // T := alpha*A*op(T), T input N×N, A M×N, result M×N
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UyUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UyTrans {
    NoTrans,
    Trans,
}

/// Overwrites the leading M×N part of T with the product.
pub fn mb01uy(
    side: Mb01UySide,
    uplo: Mb01UyUplo,
    trans: Mb01UyTrans,
    m: usize,
    n: usize,
    alpha: f64,
    t: &mut [f64],
    ldt: usize,
    a: &[f64],
    lda: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(side, Mb01UySide::Left | Mb01UySide::Right) {
        return -1;
    }
    if !matches!(uplo, Mb01UyUplo::Upper | Mb01UyUplo::Lower) {
        return -2;
    }
    if !matches!(trans, Mb01UyTrans::NoTrans | Mb01UyTrans::Trans) {
        return -3;
    }
    let k = if side == Mb01UySide::Left { m } else { n };
    if ldt < k.max(m).max(n) || lda < m.max(1) {
        return -8;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    if alpha == 0.0 {
        for j in 0..n {
            for i in 0..m {
                t[i + j * ldt] = 0.0;
            }
        }
        return 0;
    }

    let t_mat = DMatrix::from_fn(k, k, |i, j| {
        if (uplo == Mb01UyUplo::Upper && i <= j) || (uplo == Mb01UyUplo::Lower && i >= j) {
            t[i + j * ldt]
        } else {
            0.0
        }
    });
    let a_mat = DMatrix::from_fn(m, n, |i, j| a[i + j * lda]);
    let op_t = match trans {
        Mb01UyTrans::NoTrans => t_mat.clone(),
        Mb01UyTrans::Trans => t_mat.transpose(),
    };
    let result = match side {
        Mb01UySide::Left => alpha * &op_t * &a_mat,
        Mb01UySide::Right => alpha * &a_mat * &op_t,
    };
    for j in 0..n {
        for i in 0..m {
            t[i + j * ldt] = result[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01uy_left_upper() {
        let m = 2;
        let n = 2;
        let mut t = [1.0, 0.0, 0.0, 1.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; m.max(n)];
        assert_eq!(
            mb01uy(
                Mb01UySide::Left,
                Mb01UyUplo::Upper,
                Mb01UyTrans::NoTrans,
                m,
                n,
                1.0,
                &mut t,
                2,
                &a,
                2,
                &mut dwork,
            ),
            0
        );
        assert!((t[0] - 1.0).abs() < 1e-14);
        assert!((t[3] - 1.0).abs() < 1e-14);
    }
}
