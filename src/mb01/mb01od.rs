//! MB01OD — R := alpha*R + beta*(op(H)*X*op(E)' + op(E)*X*op(H)') (SLICOT MB01OD)
// H upper Hessenberg, E upper triangular, X symmetric.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OdUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OdTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R.
pub fn mb01od(
    uplo: Mb01OdUplo,
    trans: Mb01OdTrans,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    h: &[f64],
    ldh: usize,
    x: &[f64],
    ldx: usize,
    e: &[f64],
    lde: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(uplo, Mb01OdUplo::Upper | Mb01OdUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01OdTrans::NoTrans | Mb01OdTrans::Trans) {
        return -2;
    }
    if ldr < n.max(1) || ldh < n.max(1) || ldx < n.max(1) || lde < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01OdUplo::Upper && i <= j) || (uplo == Mb01OdUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let h_mat = DMatrix::from_fn(n, n, |i, j| if i <= j + 1 { h[i + j * ldh] } else { 0.0 });
    let e_mat = DMatrix::from_fn(n, n, |i, j| if i <= j { e[i + j * lde] } else { 0.0 });
    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            x_full[(i, j)] = if (uplo == Mb01OdUplo::Upper && i <= j) || (uplo == Mb01OdUplo::Lower && i >= j) {
                x[i + j * ldx]
            } else {
                x[j + i * ldx]
            };
        }
    }
    let update = match trans {
        Mb01OdTrans::NoTrans => &h_mat * &x_full * e_mat.transpose() + &e_mat * &x_full * h_mat.transpose(),
        Mb01OdTrans::Trans => h_mat.transpose() * &x_full * &e_mat + e_mat.transpose() * &x_full * &h_mat,
    };

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01OdUplo::Upper && i <= j) || (uplo == Mb01OdUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01OdUplo::Upper && i <= j) || (uplo == Mb01OdUplo::Lower && i >= j) {
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
    fn test_mb01od_upper_notrans() {
        let n = 2;
        let mut r = vec![0.0, 0.0, 0.0, 0.0];
        let h = [1.0, 0.0, 0.0, 1.0];
        let x = [1.0, 0.0, 0.0, 1.0];
        let e = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; n * n];
        assert_eq!(
            mb01od(
                Mb01OdUplo::Upper,
                Mb01OdTrans::NoTrans,
                n,
                0.0,
                1.0,
                &mut r,
                2,
                &h,
                2,
                &x,
                2,
                &e,
                2,
                &mut dwork,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
