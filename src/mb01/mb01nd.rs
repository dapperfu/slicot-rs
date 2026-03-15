//! MB01ND — Skew-symmetric rank-2 update (SLICOT MB01ND)
//
// A := alpha*x*y' - alpha*y*x' + A  (only triangular part of A updated).

use nalgebra::{DMatrix, DVector};

/// Which triangle of A is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01NdUplo {
    Upper,
    Lower,
}

/// Performs A := alpha*x*y' - alpha*y*x' + A. Only the specified triangle of A is read/written.
pub fn mb01nd(
    uplo: Mb01NdUplo,
    alpha: f64,
    x: &DVector<f64>,
    y: &DVector<f64>,
    a: &mut DMatrix<f64>,
) -> i32 {
    let n = a.nrows();
    if a.ncols() != n || x.nrows() != n || y.nrows() != n {
        return -1;
    }
    if n == 0 || alpha == 0.0 {
        return 0;
    }
    match uplo {
        Mb01NdUplo::Upper => {
            for j in 1..n {
                let temp1 = alpha * y[j];
                let temp2 = alpha * x[j];
                for i in 0..j {
                    a[(i, j)] += x[i] * temp1 - y[i] * temp2;
                }
            }
        }
        Mb01NdUplo::Lower => {
            for j in 0..n - 1 {
                let temp1 = alpha * y[j];
                let temp2 = alpha * x[j];
                for i in (j + 1)..n {
                    a[(i, j)] += x[i] * temp1 - y[i] * temp2;
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
    fn test_mb01nd_upper_zero_alpha() {
        let mut a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, -1.0, 0.0]);
        let x = DVector::from_row_slice(&[1.0, 0.0]);
        let y = DVector::from_row_slice(&[0.0, 1.0]);
        assert_eq!(mb01nd(Mb01NdUplo::Upper, 0.0, &x, &y, &mut a), 0);
        assert_eq!(a[(0, 1)], 1.0);
    }
}
