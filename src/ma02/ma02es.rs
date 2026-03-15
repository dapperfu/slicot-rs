//! MA02ES — Store by skew-symmetry (copy triangle with sign flip), skew-symmetric matrix (SLICOT MA02ES)

use nalgebra::DMatrix;

/// Which triangle is given (the other is filled by skew-symmetry; diagonal set to 0).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02EsUplo {
    /// Upper triangle given; lower is filled as -upper^T; diagonal = 0.
    Upper,
    /// Lower triangle given; upper is filled as -lower^T; diagonal = 0.
    Lower,
}

/// Fills the other triangle by skew-symmetry (A(i,j) = -A(j,i)) and sets diagonal to 0.
pub fn ma02es(uplo: Ma02EsUplo, a: &mut DMatrix<f64>) -> i32 {
    let n = a.nrows();
    if a.ncols() != n {
        return -3;
    }
    if n == 0 {
        return 0;
    }
    match uplo {
        Ma02EsUplo::Lower => {
            for i in 0..n {
                a[(i, i)] = 0.0;
                for j in (i + 1)..n {
                    a[(i, j)] = -a[(j, i)];
                }
            }
        }
        Ma02EsUplo::Upper => {
            for i in 0..n {
                a[(i, i)] = 0.0;
                for j in (i + 1)..n {
                    a[(j, i)] = -a[(i, j)];
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
    fn test_ma02es_upper() {
        let mut a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 0.0, 0.0]);
        assert_eq!(ma02es(Ma02EsUplo::Upper, &mut a), 0);
        assert_eq!(a[(0, 0)], 0.0);
        assert_eq!(a[(1, 0)], -1.0);
    }
}
