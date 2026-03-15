//! MB01MD — Skew-symmetric matrix-vector product (SLICOT MB01MD)
//
// y := alpha*A*x + beta*y, where A is n×n skew-symmetric (only triangular part stored).

use nalgebra::{DMatrix, DVector};

/// Which triangle of A is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01MdUplo {
    /// Strictly upper triangular part is stored.
    Upper,
    /// Strictly lower triangular part is stored.
    Lower,
}

/// Computes y := alpha*A*x + beta*y. A is skew-symmetric; only the specified triangle is read.
///
/// # Arguments
/// * `uplo` — which triangle of `a` is stored (strictly upper or strictly lower)
/// * `alpha` — scalar
/// * `a` — n×n matrix (only strictly upper or strictly lower part used)
/// * `x` — length n
/// * `beta` — scalar
/// * `y` — overwritten with result (length n)
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn mb01md(
    uplo: Mb01MdUplo,
    alpha: f64,
    a: &DMatrix<f64>,
    x: &DVector<f64>,
    beta: f64,
    y: &mut DVector<f64>,
) -> i32 {
    let n = a.nrows();
    if a.ncols() != n {
        return -4;
    }
    if x.nrows() != n || y.nrows() != n {
        return -8;
    }
    if n == 0 {
        return 0;
    }
    // y := beta*y
    if beta == 0.0 {
        for i in 0..n {
            y[i] = 0.0;
        }
    } else if beta != 1.0 {
        for i in 0..n {
            y[i] *= beta;
        }
    }
    if alpha == 0.0 {
        return 0;
    }
    match uplo {
        Mb01MdUplo::Upper => {
            for j in 1..n {
                let temp1 = alpha * x[j];
                let mut temp2 = 0.0;
                for i in 0..j {
                    y[i] += temp1 * a[(i, j)];
                    temp2 += a[(i, j)] * x[i];
                }
                y[j] -= alpha * temp2;
            }
        }
        Mb01MdUplo::Lower => {
            for j in 0..n - 1 {
                let temp1 = alpha * x[j];
                let mut temp2 = 0.0;
                for i in j + 1..n {
                    y[i] += temp1 * a[(i, j)];
                    temp2 += a[(i, j)] * x[i];
                }
                y[j] -= alpha * temp2;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01md_upper_zero_alpha() {
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, -1.0, 0.0]);
        let x = DVector::from_row_slice(&[1.0, 0.0]);
        let mut y = DVector::from_row_slice(&[3.0, 4.0]);
        assert_eq!(mb01md(Mb01MdUplo::Upper, 0.0, &a, &x, 1.0, &mut y), 0);
        assert_eq!(y[0], 3.0);
        assert_eq!(y[1], 4.0);
    }

    #[test]
    fn test_mb01md_upper_ax() {
        // A = [0 1; -1 0], A*x for x=[1,2] = [2, -1]
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, -1.0, 0.0]);
        let x = DVector::from_row_slice(&[1.0, 2.0]);
        let mut y = DVector::zeros(2);
        assert_eq!(mb01md(Mb01MdUplo::Upper, 1.0, &a, &x, 0.0, &mut y), 0);
        assert!((y[0] - 2.0).abs() < 1e-15);
        assert!((y[1] - (-1.0)).abs() < 1e-15);
    }
}
