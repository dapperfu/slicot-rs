//! MB01UW — A := alpha*op(H)*A or A := alpha*A*op(H) (SLICOT MB01UW)
// H upper Hessenberg.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UwSide {
    Left,  // A := alpha*op(H)*A, H is M×M
    Right, // A := alpha*A*op(H), H is N×N
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UwTrans {
    NoTrans,
    Trans,
}

/// Overwrites A with alpha*op(H)*A (Left) or alpha*A*op(H) (Right).
pub fn mb01uw(
    side: Mb01UwSide,
    trans: Mb01UwTrans,
    m: usize,
    n: usize,
    alpha: f64,
    h: &[f64],
    ldh: usize,
    a: &mut [f64],
    lda: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(side, Mb01UwSide::Left | Mb01UwSide::Right) {
        return -1;
    }
    if !matches!(trans, Mb01UwTrans::NoTrans | Mb01UwTrans::Trans) {
        return -2;
    }
    let k = if side == Mb01UwSide::Left { m } else { n };
    if ldh < k.max(1) || lda < m.max(1) {
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

    let h_mat = DMatrix::from_fn(k, k, |i, j| if i <= j + 1 { h[i + j * ldh] } else { 0.0 });
    let a_mat = DMatrix::from_fn(m, n, |i, j| a[i + j * lda]);
    let op_h = match trans {
        Mb01UwTrans::NoTrans => h_mat.clone(),
        Mb01UwTrans::Trans => h_mat.transpose(),
    };
    let result = match side {
        Mb01UwSide::Left => alpha * &op_h * &a_mat,
        Mb01UwSide::Right => alpha * &a_mat * &op_h,
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
    fn test_mb01uw_left_notrans() {
        let m = 2;
        let n = 2;
        let h = [1.0, 0.0, 0.0, 1.0];
        let mut a = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; m * n];
        assert_eq!(
            mb01uw(
                Mb01UwSide::Left,
                Mb01UwTrans::NoTrans,
                m,
                n,
                1.0,
                &h,
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
