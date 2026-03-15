//! MB01XY — Compute U'*U or L*L' in place (unblocked) (SLICOT MB01XY)

use nalgebra::DMatrix;

/// Which triangle is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01XyUplo {
    Upper, // U'*U
    Lower, // L*L'
}

/// Overwrites the triangular part of A with U'*U (Upper) or L*L' (Lower).
pub fn mb01xy(uplo: Mb01XyUplo, a: &mut DMatrix<f64>) -> i32 {
    let n = a.nrows();
    if a.ncols() != n {
        return -3;
    }
    if n == 0 {
        return 0;
    }
    let a_copy = a.clone();
    match uplo {
        Mb01XyUplo::Upper => {
            // (U'*U)(i,j) = sum_k U(k,i)*U(k,j); for upper U, U(k,i)=0 for k>i.
            for j in 0..n {
                for i in 0..=j {
                    let mut s = 0.0;
                    for k in 0..=i {
                        s += a_copy[(k, i)] * a_copy[(k, j)];
                    }
                    a[(i, j)] = s;
                }
            }
        }
        Mb01XyUplo::Lower => {
            for i in 0..n {
                for j in 0..=i {
                    let mut s = 0.0;
                    for k in j..n {
                        s += a_copy[(i, k)] * a_copy[(j, k)];
                    }
                    a[(i, j)] = s;
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
    fn test_mb01xy_upper() {
        // U = [[1,2],[0,3]] in row-major: row0=[1,2], row1=[0,3]
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
        assert_eq!(mb01xy(Mb01XyUplo::Upper, &mut a), 0);
        // U'*U = [[1,2],[2,13]]
        assert!((a[(0, 0)] - 1.0).abs() < 1e-15);
        assert!((a[(0, 1)] - 2.0).abs() < 1e-15);
        assert!((a[(1, 1)] - 13.0).abs() < 1e-15);
    }
}
