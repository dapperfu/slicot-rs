//! MB01RB — R = alpha*R + beta*op(A)*B (side L) or R = alpha*R + beta*B*op(A) (side R) (SLICOT MB01RB)
// Same operation as MB01RX; BLAS 2 style. We use full-matrix computation.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RbSide {
    Left,
    Right,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RbUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RbTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R. Same as MB01RX.
pub fn mb01rb(
    side: Mb01RbSide,
    uplo: Mb01RbUplo,
    trans: Mb01RbTrans,
    m: usize,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
) -> i32 {
    if !matches!(side, Mb01RbSide::Left | Mb01RbSide::Right) {
        return -1;
    }
    if !matches!(uplo, Mb01RbUplo::Upper | Mb01RbUplo::Lower) {
        return -2;
    }
    if !matches!(trans, Mb01RbTrans::NoTrans | Mb01RbTrans::Trans) {
        return -3;
    }
    if ldr < m.max(1) {
        return -9;
    }
    let (ar, ac) = match (side, trans) {
        (Mb01RbSide::Left, Mb01RbTrans::NoTrans) | (Mb01RbSide::Right, Mb01RbTrans::Trans) => (m, n),
        (Mb01RbSide::Left, Mb01RbTrans::Trans) | (Mb01RbSide::Right, Mb01RbTrans::NoTrans) => (n, m),
    };
    if lda < ar.max(1) {
        return -11;
    }
    let br = if side == Mb01RbSide::Left { n } else { m };
    if ldb < br.max(1) {
        return -13;
    }
    if m == 0 {
        return 0;
    }
    if beta == 0.0 || n == 0 {
        if alpha == 0.0 {
            for j in 0..m {
                for i in 0..m {
                    if (uplo == Mb01RbUplo::Upper && i <= j) || (uplo == Mb01RbUplo::Lower && i >= j) {
                        r[i + j * ldr] = 0.0;
                    }
                }
            }
        } else if alpha != 1.0 {
            for j in 0..m {
                for i in 0..m {
                    if (uplo == Mb01RbUplo::Upper && i <= j) || (uplo == Mb01RbUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let a_mat = DMatrix::from_fn(ar, ac, |i, j| a[i + j * lda]);
    let b_mat = DMatrix::from_fn(br, if side == Mb01RbSide::Left { m } else { n }, |i, j| b[i + j * ldb]);
    let op_a = match trans {
        Mb01RbTrans::NoTrans => a_mat.clone(),
        Mb01RbTrans::Trans => a_mat.transpose(),
    };
    let update = match side {
        Mb01RbSide::Left => &op_a * &b_mat,
        Mb01RbSide::Right => &b_mat * &op_a,
    };

    let mut r_full = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            r_full[(i, j)] = if (uplo == Mb01RbUplo::Upper && i <= j) || (uplo == Mb01RbUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..m {
        for j in 0..m {
            if (uplo == Mb01RbUplo::Upper && i <= j) || (uplo == Mb01RbUplo::Lower && i >= j) {
                r[i + j * ldr] = r_full[(i, j)];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01rb_left_upper() {
        let m = 2;
        let n = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        let b = [1.0, 0.0, 0.0, 1.0];
        assert_eq!(
            mb01rb(
                Mb01RbSide::Left,
                Mb01RbUplo::Upper,
                Mb01RbTrans::NoTrans,
                m,
                n,
                1.0,
                1.0,
                &mut r,
                2,
                &a,
                2,
                &b,
                2,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
