//! MB01RH — R := alpha*R + beta*op(H)*X*op(H)' (SLICOT MB01RH)
// H upper Hessenberg, R and X symmetric.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RhUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RhTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R.
pub fn mb01rh(
    uplo: Mb01RhUplo,
    trans: Mb01RhTrans,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    h: &[f64],
    ldh: usize,
    x: &[f64],
    ldx: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(uplo, Mb01RhUplo::Upper | Mb01RhUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01RhTrans::NoTrans | Mb01RhTrans::Trans) {
        return -2;
    }
    if ldr < n.max(1) || ldh < n.max(1) || ldx < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha == 0.0 {
            for i in 0..n {
                for j in 0..n {
                    if (uplo == Mb01RhUplo::Upper && i <= j) || (uplo == Mb01RhUplo::Lower && i >= j) {
                        r[i + j * ldr] = 0.0;
                    }
                }
            }
        } else if alpha != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01RhUplo::Upper && i <= j) || (uplo == Mb01RhUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let h_mat = DMatrix::from_fn(n, n, |i, j| if i <= j + 1 { h[i + j * ldh] } else { 0.0 });
    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            x_full[(i, j)] = if (uplo == Mb01RhUplo::Upper && i <= j) || (uplo == Mb01RhUplo::Lower && i >= j) {
                x[i + j * ldx]
            } else {
                x[j + i * ldx]
            };
        }
    }
    let op_h = match trans {
        Mb01RhTrans::NoTrans => h_mat.clone(),
        Mb01RhTrans::Trans => h_mat.transpose(),
    };
    let update = &op_h * &x_full * op_h.transpose();

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01RhUplo::Upper && i <= j) || (uplo == Mb01RhUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01RhUplo::Upper && i <= j) || (uplo == Mb01RhUplo::Lower && i >= j) {
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
    fn test_mb01rh_upper_notrans() {
        let n = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0];
        let h = [1.0, 0.0, 0.0, 1.0];
        let x = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; n * n];
        assert_eq!(
            mb01rh(
                Mb01RhUplo::Upper,
                Mb01RhTrans::NoTrans,
                n,
                1.0,
                1.0,
                &mut r,
                2,
                &h,
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
