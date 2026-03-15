//! MB01RT — R = alpha*R + beta*op(E)*X*op(E)' (SLICOT MB01RT)
// E upper triangular, R and X symmetric.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RtUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RtTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R. E is N×N upper triangular, X is N×N symmetric.
pub fn mb01rt(
    uplo: Mb01RtUplo,
    trans: Mb01RtTrans,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    e: &[f64],
    lde: usize,
    x: &[f64],
    ldx: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(uplo, Mb01RtUplo::Upper | Mb01RtUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01RtTrans::NoTrans | Mb01RtTrans::Trans) {
        return -2;
    }
    if ldr < n.max(1) || lde < n.max(1) || ldx < n.max(1) {
        return -8;
    }
    if n == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha == 0.0 {
            for i in 0..n {
                for j in 0..n {
                    if (uplo == Mb01RtUplo::Upper && i <= j) || (uplo == Mb01RtUplo::Lower && i >= j) {
                        r[i + j * ldr] = 0.0;
                    }
                }
            }
        } else if alpha != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01RtUplo::Upper && i <= j) || (uplo == Mb01RtUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let e_mat = DMatrix::from_fn(n, n, |i, j| if i <= j { e[i + j * lde] } else { 0.0 });
    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            x_full[(i, j)] = if (uplo == Mb01RtUplo::Upper && i <= j) || (uplo == Mb01RtUplo::Lower && i >= j) {
                x[i + j * ldx]
            } else {
                x[j + i * ldx]
            };
        }
    }
    let op_e = match trans {
        Mb01RtTrans::NoTrans => e_mat.clone(),
        Mb01RtTrans::Trans => e_mat.transpose(),
    };
    let update = &op_e * &x_full * op_e.transpose();

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01RtUplo::Upper && i <= j) || (uplo == Mb01RtUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01RtUplo::Upper && i <= j) || (uplo == Mb01RtUplo::Lower && i >= j) {
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
    fn test_mb01rt_upper_notrans() {
        let n = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0];
        let e = [1.0, 0.0, 0.0, 1.0];
        let x = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; n * n];
        assert_eq!(
            mb01rt(
                Mb01RtUplo::Upper,
                Mb01RtTrans::NoTrans,
                n,
                1.0,
                1.0,
                &mut r,
                2,
                &e,
                2,
                &x,
                2,
                &mut dwork,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
