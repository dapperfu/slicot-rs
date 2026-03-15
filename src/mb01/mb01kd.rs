//! MB01KD — Skew-symmetric rank 2k update (SLICOT MB01KD)
//!
//! C := alpha*A*B' - alpha*B*A' + beta*C  (trans='N') or
//! C := alpha*A'*B - alpha*B'*A + beta*C  (trans='T').
//! Only the strictly upper or strictly lower triangle of C is referenced.

use nalgebra::DMatrix;

/// Upper or lower triangle of C.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01KdUplo {
    /// Strictly upper triangular part.
    Upper,
    /// Strictly lower triangular part.
    Lower,
}

/// Transpose flag.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01KdTrans {
    /// C := alpha*A*B' - alpha*B*A' + beta*C (A, B are N×K).
    NoTrans,
    /// C := alpha*A'*B - alpha*B'*A + beta*C (A, B are K×N).
    Trans,
}

/// Performs the skew-symmetric rank 2k update. C is N×N (only strict upper or strict lower is used).
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn mb01kd(
    uplo: Mb01KdUplo,
    trans: Mb01KdTrans,
    alpha: f64,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    beta: f64,
    c: &mut DMatrix<f64>,
) -> i32 {
    let n = c.nrows();
    if c.ncols() != n {
        return -12;
    }
    let k = match trans {
        Mb01KdTrans::NoTrans => {
            if a.nrows() != n || b.nrows() != n {
                return -7;
            }
            let k = a.ncols();
            if b.ncols() != k {
                return -9;
            }
            k
        }
        Mb01KdTrans::Trans => {
            if a.ncols() != n || b.ncols() != n {
                return -7;
            }
            let k = a.nrows();
            if b.nrows() != k {
                return -9;
            }
            k
        }
    };
    if n <= 1 || (alpha == 0.0 || k == 0) && beta == 1.0 {
        return 0;
    }
    if alpha == 0.0 {
        if beta == 0.0 {
            match uplo {
                Mb01KdUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            c[(i, j)] = 0.0;
                        }
                    }
                }
                Mb01KdUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            c[(i, j)] = 0.0;
                        }
                    }
                }
            }
        } else {
            match uplo {
                Mb01KdUplo::Upper => {
                    for j in 1..n {
                        for i in 0..j {
                            c[(i, j)] *= beta;
                        }
                    }
                }
                Mb01KdUplo::Lower => {
                    for j in 0..n - 1 {
                        for i in j + 1..n {
                            c[(i, j)] *= beta;
                        }
                    }
                }
            }
        }
        return 0;
    }
    match (trans, uplo) {
        (Mb01KdTrans::NoTrans, Mb01KdUplo::Upper) => {
            for j in 1..n {
                if beta == 0.0 {
                    for i in 0..j {
                        c[(i, j)] = 0.0;
                    }
                } else if beta != 1.0 {
                    for i in 0..j {
                        c[(i, j)] *= beta;
                    }
                }
                for l in 0..k {
                    let temp1 = alpha * b[(j, l)];
                    let temp2 = alpha * a[(j, l)];
                    for i in 0..j {
                        c[(i, j)] += a[(i, l)] * temp1 - b[(i, l)] * temp2;
                    }
                }
            }
        }
        (Mb01KdTrans::NoTrans, Mb01KdUplo::Lower) => {
            for j in 0..n - 1 {
                if beta == 0.0 {
                    for i in j + 1..n {
                        c[(i, j)] = 0.0;
                    }
                } else if beta != 1.0 {
                    for i in j + 1..n {
                        c[(i, j)] *= beta;
                    }
                }
                for l in 0..k {
                    let temp1 = alpha * b[(j, l)];
                    let temp2 = alpha * a[(j, l)];
                    for i in j + 1..n {
                        c[(i, j)] += a[(i, l)] * temp1 - b[(i, l)] * temp2;
                    }
                }
            }
        }
        (Mb01KdTrans::Trans, Mb01KdUplo::Upper) => {
            for j in 1..n {
                for i in 0..j {
                    let mut temp1 = 0.0_f64;
                    let mut temp2 = 0.0_f64;
                    for l in 0..k {
                        temp1 += a[(l, i)] * b[(l, j)];
                        temp2 += b[(l, i)] * a[(l, j)];
                    }
                    c[(i, j)] = if beta == 0.0 {
                        alpha * temp1 - alpha * temp2
                    } else {
                        beta * c[(i, j)] + alpha * temp1 - alpha * temp2
                    };
                }
            }
        }
        (Mb01KdTrans::Trans, Mb01KdUplo::Lower) => {
            for j in 0..n - 1 {
                for i in j + 1..n {
                    let mut temp1 = 0.0_f64;
                    let mut temp2 = 0.0_f64;
                    for l in 0..k {
                        temp1 += a[(l, i)] * b[(l, j)];
                        temp2 += b[(l, i)] * a[(l, j)];
                    }
                    c[(i, j)] = if beta == 0.0 {
                        alpha * temp1 - alpha * temp2
                    } else {
                        beta * c[(i, j)] + alpha * temp1 - alpha * temp2
                    };
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01kd_upper_notrans() {
        let a = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let b = DMatrix::from_row_slice(2, 1, &[0.0, 1.0]);
        let mut c = DMatrix::zeros(2, 2);
        assert_eq!(
            mb01kd(
                Mb01KdUplo::Upper,
                Mb01KdTrans::NoTrans,
                1.0,
                &a,
                &b,
                0.0,
                &mut c
            ),
            0
        );
        // C(0,1) = A(0,0)*B(1,0) - B(0,0)*A(1,0) = 1*1 - 0*0 = 1
        assert!((c[(0, 1)] - 1.0).abs() < 1e-10);
    }
}
