//! MB01RY — R = alpha*R + beta*op(H)*B or R = alpha*R + beta*B*op(H) (SLICOT MB01RY)
// H is upper Hessenberg, R and result M×M.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RySide {
    Left,  // R = alpha*R + beta*op(H)*B
    Right, // R = alpha*R + beta*B*op(H)
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RyUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RyTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R. H is M×M upper Hessenberg, B is M×M.
pub fn mb01ry(
    side: Mb01RySide,
    uplo: Mb01RyUplo,
    trans: Mb01RyTrans,
    m: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    h: &[f64],
    ldh: usize,
    b: &[f64],
    ldb: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(side, Mb01RySide::Left | Mb01RySide::Right) {
        return -1;
    }
    if !matches!(uplo, Mb01RyUplo::Upper | Mb01RyUplo::Lower) {
        return -2;
    }
    if !matches!(trans, Mb01RyTrans::NoTrans | Mb01RyTrans::Trans) {
        return -3;
    }
    if ldr < m.max(1) || ldh < m.max(1) || ldb < m.max(1) {
        return -8;
    }
    if m == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha == 0.0 {
            for j in 0..m {
                for i in 0..m {
                    if (uplo == Mb01RyUplo::Upper && i <= j) || (uplo == Mb01RyUplo::Lower && i >= j) {
                        r[i + j * ldr] = 0.0;
                    }
                }
            }
        } else if alpha != 1.0 {
            for j in 0..m {
                for i in 0..m {
                    if (uplo == Mb01RyUplo::Upper && i <= j) || (uplo == Mb01RyUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let h_mat = DMatrix::from_fn(m, m, |i, j| {
        if i <= j + 1 {
            h[i + j * ldh]
        } else {
            0.0
        }
    });
    let b_mat = DMatrix::from_fn(m, m, |i, j| b[i + j * ldb]);
    let op_h = match trans {
        Mb01RyTrans::NoTrans => h_mat.clone(),
        Mb01RyTrans::Trans => h_mat.transpose(),
    };
    let update = match side {
        Mb01RySide::Left => &op_h * &b_mat,
        Mb01RySide::Right => &b_mat * &op_h,
    };

    let mut r_full = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            r_full[(i, j)] = if (uplo == Mb01RyUplo::Upper && i <= j) || (uplo == Mb01RyUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..m {
        for j in 0..m {
            if (uplo == Mb01RyUplo::Upper && i <= j) || (uplo == Mb01RyUplo::Lower && i >= j) {
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
    fn test_mb01ry_left_upper() {
        let m = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0];
        let h = [1.0, 0.0, 0.0, 1.0]; // Hessenberg (diagonal)
        let b = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; m];
        assert_eq!(
            mb01ry(
                Mb01RySide::Left,
                Mb01RyUplo::Upper,
                Mb01RyTrans::NoTrans,
                m,
                1.0,
                1.0,
                &mut r,
                2,
                &h,
                2,
                &b,
                2,
                &mut dwork,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
