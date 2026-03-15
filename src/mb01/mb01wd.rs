//! MB01WD — R = alpha*(op(A)'*op(T)'*op(T)*op(A) ± ...) + beta*R (SLICOT MB01WD)
// Continuous: R = alpha*(op(A)'*op(T)'*op(T) + op(T)'*op(T)*op(A)) + beta*R
// Discrete:   R = alpha*(op(A)'*op(T)'*op(T)*op(A) - op(T)'*op(T)) + beta*R

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01WdDico {
    Continuous,
    Discrete,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01WdUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01WdTrans {
    NoTrans,
    Trans,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01WdHess {
    Full,
    Hessenberg,
}

/// Overwrites the triangle of R. We use full matrix path (HESS ignored for structure).
pub fn mb01wd(
    dico: Mb01WdDico,
    uplo: Mb01WdUplo,
    trans: Mb01WdTrans,
    _hess: Mb01WdHess,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    a: &[f64],
    lda: usize,
    t: &[f64],
    ldt: usize,
) -> i32 {
    if !matches!(dico, Mb01WdDico::Continuous | Mb01WdDico::Discrete) {
        return -1;
    }
    if !matches!(uplo, Mb01WdUplo::Upper | Mb01WdUplo::Lower) {
        return -2;
    }
    if !matches!(trans, Mb01WdTrans::NoTrans | Mb01WdTrans::Trans) {
        return -3;
    }
    if ldr < n.max(1) || lda < n.max(1) || ldt < n.max(1) {
        return -9;
    }
    if n == 0 {
        return 0;
    }
    if alpha == 0.0 {
        if beta == 0.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01WdUplo::Upper && i <= j) || (uplo == Mb01WdUplo::Lower && i >= j) {
                        r[i + j * ldr] = 0.0;
                    }
                }
            }
        } else if beta != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01WdUplo::Upper && i <= j) || (uplo == Mb01WdUplo::Lower && i >= j) {
                        r[i + j * ldr] *= beta;
                    }
                }
            }
        }
        return 0;
    }

    let t_mat = DMatrix::from_fn(n, n, |i, j| {
        if (uplo == Mb01WdUplo::Upper && i <= j) || (uplo == Mb01WdUplo::Lower && i >= j) {
            t[i + j * ldt]
        } else {
            0.0
        }
    });
    let a_mat = DMatrix::from_fn(n, n, |i, j| a[i + j * lda]);

    let (op_a, op_t) = match trans {
        Mb01WdTrans::NoTrans => (a_mat.clone(), t_mat.clone()),
        Mb01WdTrans::Trans => (a_mat.transpose(), t_mat.transpose()),
    };

    let ttt = op_t.transpose() * &op_t; // op(T)'*op(T)

    let update = match dico {
        Mb01WdDico::Continuous => {
            // R = alpha*(op(A)'*op(T)'*op(T) + op(T)'*op(T)*op(A)) + beta*R
            &op_a.transpose() * &ttt + &ttt * &op_a
        }
        Mb01WdDico::Discrete => {
            // R = alpha*(op(A)'*op(T)'*op(T)*op(A) - op(T)'*op(T)) + beta*R
            &op_a.transpose() * &ttt * &op_a - &ttt
        }
    };

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01WdUplo::Upper && i <= j) || (uplo == Mb01WdUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * update + beta * r_full;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01WdUplo::Upper && i <= j) || (uplo == Mb01WdUplo::Lower && i >= j) {
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
    fn test_mb01wd_discrete_upper_notrans() {
        let n = 2;
        let mut r = vec![0.0, 0.0, 0.0, 0.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        let t = [1.0, 0.0, 0.0, 1.0];
        assert_eq!(
            mb01wd(
                Mb01WdDico::Discrete,
                Mb01WdUplo::Upper,
                Mb01WdTrans::NoTrans,
                Mb01WdHess::Full,
                n,
                1.0,
                0.0,
                &mut r,
                2,
                &a,
                2,
                &t,
                2,
            ),
            0
        );
        // R = A'*T'*T*A - T'*T = I*I*I*I - I = 0
        assert!((r[0].abs()) < 1e-14);
        assert!((r[3].abs()) < 1e-14);
    }
}
