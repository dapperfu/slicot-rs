//! MB01RW — A := op(Z)*A*op(Z)' (SLICOT MB01RW)
// Congruence transform. A is N×N (input), result is M×M. Z is M×N (trans N) or N×M (trans T).

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RwUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01RwTrans {
    NoTrans, // op(Z)=Z, Z is M×N, result Z*A*Z'
    Trans,   // op(Z)=Z', Z is N×M, result Z'*A*Z
}

/// Overwrites the leading M×M part of A with op(Z)*A*op(Z)'. A on entry is N×N (triangular).
pub fn mb01rw(
    uplo: Mb01RwUplo,
    trans: Mb01RwTrans,
    m: usize,
    n: usize,
    a: &mut [f64],
    lda: usize,
    z: &[f64],
    ldz: usize,
    _dwork: &mut [f64],
) -> i32 {
    if !matches!(uplo, Mb01RwUplo::Upper | Mb01RwUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01RwTrans::NoTrans | Mb01RwTrans::Trans) {
        return -2;
    }
    if lda < m.max(n).max(1) {
        return -6;
    }
    let ldz_need = if trans == Mb01RwTrans::NoTrans { m } else { n };
    if ldz < ldz_need.max(1) {
        return -8;
    }
    if n == 0 || m == 0 {
        return 0;
    }

    let mut a_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            a_full[(i, j)] = if (uplo == Mb01RwUplo::Upper && i <= j) || (uplo == Mb01RwUplo::Lower && i >= j) {
                a[i + j * lda]
            } else {
                a[j + i * lda]
            };
        }
    }
    let z_mat = match trans {
        Mb01RwTrans::NoTrans => DMatrix::from_fn(m, n, |i, j| z[i + j * ldz]),
        Mb01RwTrans::Trans => DMatrix::from_fn(m, n, |i, j| z[j + i * ldz]),
    };
    let result = &z_mat * &a_full * z_mat.transpose();
    for i in 0..m {
        for j in 0..m {
            if (uplo == Mb01RwUplo::Upper && i <= j) || (uplo == Mb01RwUplo::Lower && i >= j) {
                a[i + j * lda] = result[(i, j)];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01rw_upper_notrans() {
        let m = 2;
        let n = 2;
        let mut a = vec![1.0, 0.0, 0.0, 1.0];
        let z = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; n];
        assert_eq!(
            mb01rw(
                Mb01RwUplo::Upper,
                Mb01RwTrans::NoTrans,
                m,
                n,
                &mut a,
                2,
                &z,
                2,
                &mut dwork,
            ),
            0
        );
        assert!((a[0] - 1.0).abs() < 1e-14);
        assert!((a[3] - 1.0).abs() < 1e-14);
    }
}
