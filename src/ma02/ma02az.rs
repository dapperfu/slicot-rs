//! MA02AZ — (Conjugate) transpose into another matrix (SLICOT MA02AZ)
//
// Complex in SLICOT; this is the real version (transpose = conjugate transpose).

use nalgebra::DMatrix;

/// Part of the matrix to transpose.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02AzJob {
    Upper,
    Lower,
    All,
}

/// Transposes the specified part of A into B. B must be N×M. Real version (trans and conj trans are the same).
pub fn ma02az(_trans: char, job: Ma02AzJob, a: &DMatrix<f64>, b: &mut DMatrix<f64>) -> i32 {
    let m = a.nrows();
    let n = a.ncols();
    if b.nrows() != n || b.ncols() != m {
        return -6;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    match job {
        Ma02AzJob::Upper => {
            for j in 0..n {
                for i in 0..m.min(j + 1) {
                    b[(j, i)] = a[(i, j)];
                }
            }
        }
        Ma02AzJob::Lower => {
            for j in 0..n {
                for i in j..m {
                    b[(j, i)] = a[(i, j)];
                }
            }
        }
        Ma02AzJob::All => {
            for j in 0..n {
                for i in 0..m {
                    b[(j, i)] = a[(i, j)];
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
    fn test_ma02az_all() {
        let a = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let mut b = DMatrix::zeros(3, 2);
        assert_eq!(ma02az('T', Ma02AzJob::All, &a, &mut b), 0);
        assert_eq!(b[(0, 0)], 1.0);
        assert_eq!(b[(2, 1)], 6.0);
    }
}
