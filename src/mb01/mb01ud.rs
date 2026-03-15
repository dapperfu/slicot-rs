//! MB01UD — B = alpha*op(H)*A or B = alpha*A*op(H) (SLICOT MB01UD)
// H upper Hessenberg.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UdSide {
    Left,  // B = alpha*op(H)*A, H is M×M
    Right, // B = alpha*A*op(H), H is N×N
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UdTrans {
    NoTrans,
    Trans,
}

/// Computes B = alpha*op(H)*A (side Left) or B = alpha*A*op(H) (side Right).
pub fn mb01ud(
    side: Mb01UdSide,
    trans: Mb01UdTrans,
    m: usize,
    n: usize,
    alpha: f64,
    h: &[f64],
    ldh: usize,
    a: &[f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
) -> i32 {
    if !matches!(side, Mb01UdSide::Left | Mb01UdSide::Right) {
        return -1;
    }
    if !matches!(trans, Mb01UdTrans::NoTrans | Mb01UdTrans::Trans) {
        return -2;
    }
    let k = if side == Mb01UdSide::Left { m } else { n };
    if ldh < k.max(1) || lda < m.max(1) || ldb < m.max(1) {
        return -10;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    if alpha == 0.0 {
        for j in 0..n {
            for i in 0..m {
                b[i + j * ldb] = 0.0;
            }
        }
        return 0;
    }

    let h_mat = DMatrix::from_fn(k, k, |i, j| {
        if i <= j + 1 {
            h[i + j * ldh]
        } else {
            0.0
        }
    });
    let a_mat = DMatrix::from_fn(m, n, |i, j| a[i + j * lda]);
    let op_h = match trans {
        Mb01UdTrans::NoTrans => h_mat.clone(),
        Mb01UdTrans::Trans => h_mat.transpose(),
    };
    let result = match side {
        Mb01UdSide::Left => alpha * &op_h * &a_mat,
        Mb01UdSide::Right => alpha * &a_mat * &op_h,
    };
    for j in 0..n {
        for i in 0..m {
            b[i + j * ldb] = result[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01ud_left_notrans() {
        let m = 2;
        let n = 2;
        let h = [1.0, 0.0, 0.0, 1.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        let mut b = vec![0.0; 4];
        assert_eq!(
            mb01ud(
                Mb01UdSide::Left,
                Mb01UdTrans::NoTrans,
                m,
                n,
                1.0,
                &h,
                2,
                &a,
                2,
                &mut b,
                2,
            ),
            0
        );
        assert!((b[0] - 1.0).abs() < 1e-14);
        assert!((b[3] - 1.0).abs() < 1e-14);
    }
}
