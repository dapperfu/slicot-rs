//! MB01RX — R = alpha*R + beta*op(A)*B (side L) or R = alpha*R + beta*B*op(A) (side R) (SLICOT MB01RX)

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RxSide {
    Left,  // R = alpha*R + beta*op(A)*B
    Right, // R = alpha*R + beta*B*op(A)
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RxUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RxTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R. R is M×M; for Side Left, A M×N and B N×M (trans N) or A N×M and B N×M (trans T).
pub fn mb01rx(
    side: Mb01RxSide,
    uplo: Mb01RxUplo,
    trans: Mb01RxTrans,
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
    if !matches!(side, Mb01RxSide::Left | Mb01RxSide::Right) {
        return -1;
    }
    if !matches!(uplo, Mb01RxUplo::Upper | Mb01RxUplo::Lower) {
        return -2;
    }
    if !matches!(trans, Mb01RxTrans::NoTrans | Mb01RxTrans::Trans) {
        return -3;
    }
    if ldr < m.max(1) {
        return -9;
    }
    let (ar, ac) = match (side, trans) {
        (Mb01RxSide::Left, Mb01RxTrans::NoTrans) | (Mb01RxSide::Right, Mb01RxTrans::Trans) => (m, n),
        (Mb01RxSide::Left, Mb01RxTrans::Trans) | (Mb01RxSide::Right, Mb01RxTrans::NoTrans) => (n, m),
    };
    if lda < ar.max(1) {
        return -11;
    }
    let (br, _bc) = if side == Mb01RxSide::Left {
        (n, m)
    } else {
        (m, n)
    };
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
                    if (uplo == Mb01RxUplo::Upper && i <= j) || (uplo == Mb01RxUplo::Lower && i >= j) {
                        r[i + j * ldr] = 0.0;
                    }
                }
            }
        } else if alpha != 1.0 {
            for j in 0..m {
                for i in 0..m {
                    if (uplo == Mb01RxUplo::Upper && i <= j) || (uplo == Mb01RxUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let a_mat = DMatrix::from_fn(ar, ac, |i, j| a[i + j * lda]);
    let b_mat = DMatrix::from_fn(br, if side == Mb01RxSide::Left { m } else { n }, |i, j| b[i + j * ldb]);

    let op_a = match trans {
        Mb01RxTrans::NoTrans => a_mat.clone(),
        Mb01RxTrans::Trans => a_mat.transpose(),
    };

    let update = match side {
        Mb01RxSide::Left => &op_a * &b_mat,
        Mb01RxSide::Right => &b_mat * &op_a,
    };

    let mut r_full = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            r_full[(i, j)] = if (uplo == Mb01RxUplo::Upper && i <= j) || (uplo == Mb01RxUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..m {
        for j in 0..m {
            if (uplo == Mb01RxUplo::Upper && i <= j) || (uplo == Mb01RxUplo::Lower && i >= j) {
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
    fn test_mb01rx_left_upper_notrans() {
        let m = 2;
        let n = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        let b = [1.0, 0.0, 0.0, 1.0];
        assert_eq!(
            mb01rx(
                Mb01RxSide::Left,
                Mb01RxUplo::Upper,
                Mb01RxTrans::NoTrans,
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
        // R = I + A*B = I + I = 2*I
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
