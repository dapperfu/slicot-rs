//! MB01OC — R := alpha*R + beta*H*X + beta*X*H' (or H'*X + X*H) (SLICOT MB01OC)
// H upper Hessenberg, R and X symmetric.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OcUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OcTrans {
    NoTrans, // R := alpha*R + beta*H*X + beta*X*H'
    Trans,   // R := alpha*R + beta*H'*X + beta*X*H
}

/// Overwrites the triangle of R.
pub fn mb01oc(
    uplo: Mb01OcUplo,
    trans: Mb01OcTrans,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    h: &[f64],
    ldh: usize,
    x: &[f64],
    ldx: usize,
) -> i32 {
    if !matches!(uplo, Mb01OcUplo::Upper | Mb01OcUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01OcTrans::NoTrans | Mb01OcTrans::Trans) {
        return -2;
    }
    if ldr < n.max(1) || ldh < n.max(1) || ldx < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01OcUplo::Upper && i <= j) || (uplo == Mb01OcUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let h_mat = DMatrix::from_fn(n, n, |i, j| {
        if i <= j + 1 {
            h[i + j * ldh]
        } else {
            0.0
        }
    });
    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            x_full[(i, j)] = if (uplo == Mb01OcUplo::Upper && i <= j) || (uplo == Mb01OcUplo::Lower && i >= j) {
                x[i + j * ldx]
            } else {
                x[j + i * ldx]
            };
        }
    }
    let update = match trans {
        Mb01OcTrans::NoTrans => &h_mat * &x_full + &x_full * h_mat.transpose(),
        Mb01OcTrans::Trans => h_mat.transpose() * &x_full + &x_full * &h_mat,
    };

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01OcUplo::Upper && i <= j) || (uplo == Mb01OcUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01OcUplo::Upper && i <= j) || (uplo == Mb01OcUplo::Lower && i >= j) {
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
    fn test_mb01oc_upper_notrans() {
        let n = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0];
        let h = [1.0, 0.0, 0.0, 1.0];
        let x = [1.0, 0.0, 0.0, 1.0];
        assert_eq!(
            mb01oc(
                Mb01OcUplo::Upper,
                Mb01OcTrans::NoTrans,
                n,
                1.0,
                1.0,
                &mut r,
                2,
                &h,
                2,
                &x,
                2,
            ),
            0
        );
        // R = I + H*I + I*H' = I + H + H' = I + 2*I = 3*I for H=I
        assert!((r[0] - 3.0).abs() < 1e-14);
        assert!((r[3] - 3.0).abs() < 1e-14);
    }
}
