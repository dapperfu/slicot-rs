//! MA02DD — Pack or unpack upper/lower triangle of a symmetric matrix (SLICOT MA02DD)

use nalgebra::DMatrix;

/// Pack or unpack.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02DdJob {
    /// Pack A into AP.
    Pack,
    /// Unpack AP into A.
    Unpack,
}

/// Upper or lower triangle.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02DdUplo {
    Upper,
    Lower,
}

/// Packs or unpacks the upper or lower triangle of a symmetric matrix. AP has length n*(n+1)/2;
/// column-wise order: U: (1,1),(1,2),(2,2),...,(1,n),...,(n,n); L: (1,1),(2,1),...,(n,1),(2,2),...
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn ma02dd(
    job: Ma02DdJob,
    uplo: Ma02DdUplo,
    a: &mut DMatrix<f64>,
    ap: &mut [f64],
) -> i32 {
    let n = a.nrows();
    if a.ncols() != n {
        return -4;
    }
    let len = n * (n + 1) / 2;
    if ap.len() < len {
        return -6;
    }
    if n == 0 {
        return 0;
    }
    match (job, uplo) {
        (Ma02DdJob::Pack, Ma02DdUplo::Upper) => {
            let mut idx = 0;
            for j in 0..n {
                for i in 0..=j {
                    ap[idx] = a[(i, j)];
                    idx += 1;
                }
            }
        }
        (Ma02DdJob::Pack, Ma02DdUplo::Lower) => {
            let mut idx = 0;
            for j in 0..n {
                for i in j..n {
                    ap[idx] = a[(i, j)];
                    idx += 1;
                }
            }
        }
        (Ma02DdJob::Unpack, Ma02DdUplo::Upper) => {
            let mut idx = 0;
            for j in 0..n {
                for i in 0..=j {
                    a[(i, j)] = ap[idx];
                    idx += 1;
                }
            }
        }
        (Ma02DdJob::Unpack, Ma02DdUplo::Lower) => {
            let mut idx = 0;
            for j in 0..n {
                for i in j..n {
                    a[(i, j)] = ap[idx];
                    idx += 1;
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
    fn test_ma02dd_pack_unpack_upper() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 2.0, 4.0]);
        let mut ap = vec![0.0; 3];
        assert_eq!(ma02dd(Ma02DdJob::Pack, Ma02DdUplo::Upper, &mut a, &mut ap), 0);
        assert_eq!(ap[0], 1.0);
        assert_eq!(ap[1], 2.0);
        assert_eq!(ap[2], 4.0);
        a.fill(0.0);
        assert_eq!(ma02dd(Ma02DdJob::Unpack, Ma02DdUplo::Upper, &mut a, &mut ap), 0);
        assert_eq!(a[(0, 0)], 1.0);
        assert_eq!(a[(0, 1)], 2.0);
        assert_eq!(a[(1, 1)], 4.0);
    }
}
