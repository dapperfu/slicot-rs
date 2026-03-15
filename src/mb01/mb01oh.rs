//! MB01OH — R := alpha*R + beta*H*A' + beta*A*H' (SLICOT MB01OH)
// H and A upper Hessenberg.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OhUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OhTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R.
pub fn mb01oh(
    uplo: Mb01OhUplo,
    trans: Mb01OhTrans,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    h: &[f64],
    ldh: usize,
    a: &[f64],
    lda: usize,
) -> i32 {
    if !matches!(uplo, Mb01OhUplo::Upper | Mb01OhUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01OhTrans::NoTrans | Mb01OhTrans::Trans) {
        return -2;
    }
    if ldr < n.max(1) || ldh < n.max(1) || lda < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01OhUplo::Upper && i <= j) || (uplo == Mb01OhUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let h_mat = DMatrix::from_fn(n, n, |i, j| if i <= j + 1 { h[i + j * ldh] } else { 0.0 });
    let a_mat = DMatrix::from_fn(n, n, |i, j| if i <= j + 1 { a[i + j * lda] } else { 0.0 });
    let update = match trans {
        Mb01OhTrans::NoTrans => &h_mat * a_mat.transpose() + &a_mat * h_mat.transpose(),
        Mb01OhTrans::Trans => h_mat.transpose() * &a_mat + a_mat.transpose() * &h_mat,
    };

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01OhUplo::Upper && i <= j) || (uplo == Mb01OhUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01OhUplo::Upper && i <= j) || (uplo == Mb01OhUplo::Lower && i >= j) {
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
    fn test_mb01oh_upper_notrans() {
        let n = 2;
        let mut r = vec![0.0, 0.0, 0.0, 0.0];
        let h = [1.0, 0.0, 0.0, 1.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        assert_eq!(
            mb01oh(
                Mb01OhUplo::Upper,
                Mb01OhTrans::NoTrans,
                n,
                0.0,
                1.0,
                &mut r,
                2,
                &h,
                2,
                &a,
                2,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
