//! MA02ED — Store by symmetry (copy triangle to the other), symmetric matrix (SLICOT MA02ED)

use nalgebra::DMatrix;

/// Which triangle is given (the other is filled by symmetry).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02EdUplo {
    /// Upper triangle given; lower is filled.
    Upper,
    /// Lower triangle given; upper is filled.
    Lower,
}

/// Fills the other triangle of A by symmetry. A must be n×n.
pub fn ma02ed(uplo: Ma02EdUplo, a: &mut DMatrix<f64>) -> i32 {
    let n = a.nrows();
    if a.ncols() != n {
        return -3;
    }
    if n == 0 {
        return 0;
    }
    match uplo {
        Ma02EdUplo::Lower => {
            for j in 1..n {
                for i in 0..j {
                    a[(i, j)] = a[(j, i)];
                }
            }
        }
        Ma02EdUplo::Upper => {
            for j in 1..n {
                for i in 0..j {
                    a[(j, i)] = a[(i, j)];
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
    fn test_ma02ed_upper() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
        assert_eq!(ma02ed(Ma02EdUplo::Upper, &mut a), 0);
        assert_eq!(a[(1, 0)], 2.0);
    }
}
