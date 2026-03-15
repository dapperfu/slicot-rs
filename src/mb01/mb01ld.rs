//! MB01LD — R = alpha*R + beta*op(A)*X*op(A)' (SLICOT MB01LD)
// R and X are skew-symmetric (strictly upper triangle stored).

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01LdUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01LdTrans {
    NoTrans,
    Trans,
}

/// Overwrites the strictly triangular part of R. R and X are skew-symmetric.
pub fn mb01ld(
    uplo: Mb01LdUplo,
    trans: Mb01LdTrans,
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
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(uplo, Mb01LdUplo::Upper | Mb01LdUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01LdTrans::NoTrans | Mb01LdTrans::Trans) {
        return -2;
    }
    if m == 0 {
        return 0;
    }
    if m <= 1 || n <= 1 {
        return 0;
    }
    if ldr < m.max(1) || lda < m.max(1) || ldx < n.max(1) {
        return -8;
    }
    if beta == 0.0 {
        if alpha != 1.0 {
            for j in 0..m {
                for i in 0..m {
                    if (uplo == Mb01LdUplo::Upper && i < j) || (uplo == Mb01LdUplo::Lower && i > j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            if i < j {
                x_full[(i, j)] = x[i + j * ldx];
                x_full[(j, i)] = -x[i + j * ldx];
            }
        }
    }
    let a_mat = match trans {
        Mb01LdTrans::NoTrans => DMatrix::from_fn(m, n, |i, j| a[i + j * lda]),
        Mb01LdTrans::Trans => DMatrix::from_fn(m, n, |i, j| a[j + i * lda]),
    };
    let update_full = &a_mat * &x_full * a_mat.transpose();
    let update_skew = (&update_full - update_full.transpose()) * 0.5;

    let mut r_full = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            if (uplo == Mb01LdUplo::Upper && i < j) || (uplo == Mb01LdUplo::Lower && i > j) {
                r_full[(i, j)] = r[i + j * ldr];
            } else if (uplo == Mb01LdUplo::Upper && i > j) || (uplo == Mb01LdUplo::Lower && i < j) {
                r_full[(i, j)] = -r[j + i * ldr];
            }
            // diag stays 0
        }
    }
    r_full = alpha * r_full + beta * update_skew;
    for i in 0..m {
        for j in 0..m {
            if (uplo == Mb01LdUplo::Upper && i < j) || (uplo == Mb01LdUplo::Lower && i > j) {
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
    fn test_mb01ld_upper_notrans() {
        let m = 2;
        let n = 2;
        let mut r = vec![0.0, 0.0, 0.0, 0.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        let x = [0.0, 0.0, 1.0, 0.0]; // skew: upper X(0,1)=1, so X(1,0)=-1; column-major x[0+1*2]=1
        assert_eq!(
            mb01ld(
                Mb01LdUplo::Upper,
                Mb01LdTrans::NoTrans,
                m,
                n,
                0.0,
                1.0,
                &mut r,
                2,
                &a,
                2,
                &x,
                2,
                &mut [],
            ),
            0
        );
        // A*X*A' with A=I, X skew [[0,1],[-1,0]] gives [[0,1],[-1,0]]; upper (0,1) at r[0+1*2]=r[2]
        assert!((r[2].abs() - 1.0).abs() < 1e-14);
    }
}
