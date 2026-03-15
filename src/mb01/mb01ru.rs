//! MB01RU — R = alpha*R + beta*op(A)*X*op(A)' (SLICOT MB01RU)
//
// R, X symmetric (triangular stored); A general. op(A) = A or A'.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RuUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RuTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R with alpha*R + beta*op(A)*X*op(A)'.
/// R is M×M, X is N×N, A is M×N for NoTrans else N×M. dwork length >= M*N when beta != 0.
pub fn mb01ru(
    uplo: Mb01RuUplo,
    trans: Mb01RuTrans,
    m: usize,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    a: &[f64],
    lda: usize,
    x: &[f64],
    ldx: usize,
    dwork: &mut [f64],
) -> i32 {
    if !matches!(uplo, Mb01RuUplo::Upper | Mb01RuUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01RuTrans::NoTrans | Mb01RuTrans::Trans) {
        return -2;
    }
    if ldr < m.max(1) {
        return -8;
    }
    let (a_rows, _a_cols) = match trans {
        Mb01RuTrans::NoTrans => (m, n),
        Mb01RuTrans::Trans => (n, m),
    };
    if lda < a_rows.max(1) {
        return -10;
    }
    if ldx < n.max(1) {
        return -12;
    }
    if beta != 0.0 && dwork.len() < m * n {
        return -14;
    }
    if m == 0 {
        return 0;
    }
    if beta == 0.0 || n == 0 {
        if alpha == 0.0 {
            for i in 0..m {
                for j in 0..m {
                    r[i + j * ldr] = 0.0;
                }
            }
        } else if alpha != 1.0 {
            for j in 0..m {
                for i in 0..m {
                    r[i + j * ldr] *= alpha;
                }
            }
        }
        return 0;
    }

    // Build full R and X from triangular (column-major)
    let mut r_full = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            r_full[(i, j)] = if (uplo == Mb01RuUplo::Upper && i <= j) || (uplo == Mb01RuUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            x_full[(i, j)] = if (uplo == Mb01RuUplo::Upper && i <= j) || (uplo == Mb01RuUplo::Lower && i >= j) {
                x[i + j * ldx]
            } else {
                x[j + i * ldx]
            };
        }
    }

    let a_mat = match trans {
        Mb01RuTrans::NoTrans => DMatrix::from_fn(m, n, |i, j| a[i + j * lda]),
        Mb01RuTrans::Trans => DMatrix::from_fn(m, n, |i, j| a[j + i * lda]),
    };

    let update = &a_mat * &x_full * a_mat.transpose();
    r_full = alpha * r_full + beta * update;

    for i in 0..m {
        for j in 0..m {
            if (uplo == Mb01RuUplo::Upper && i <= j) || (uplo == Mb01RuUplo::Lower && i >= j) {
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
    fn test_mb01ru_upper_notrans() {
        let m = 2;
        let n = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0]; // identity upper
        let a = [1.0, 0.0, 0.0, 1.0]; // identity 2x2
        let x = [1.0, 0.0, 0.0, 1.0]; // identity upper
        let mut dwork = vec![0.0; m * n];
        assert_eq!(
            mb01ru(
                Mb01RuUplo::Upper,
                Mb01RuTrans::NoTrans,
                m,
                n,
                1.0,
                1.0,
                &mut r,
                2,
                &a,
                2,
                &x,
                2,
                &mut dwork,
            ),
            0
        );
        // R = I + A*I*A' = I + I = 2*I (upper: r[0]=(0,0), r[2]=(0,1), r[3]=(1,1))
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[2] - 0.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
