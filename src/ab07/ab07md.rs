//! AB07MD — Find the dual of a given state-space representation (SLICOT AB07MD).
//!
//! Dual of (A,B,C,D) is (A', C', B', D').

use nalgebra::DMatrix;
use std::cmp::{max, min};

/// JobD: whether direct transmission matrix D is present.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobD {
    /// D is present.
    D,
    /// D is assumed zero (not referenced).
    Z,
}

/// Find the dual of (A,B,C,D): output (A', C', B', D').
///
/// # Arguments
/// * `jobd` — D present ('D') or zero ('Z').
/// * `n` — order of A (rows/cols of A).
/// * `m` — number of inputs (columns of B).
/// * `p` — number of outputs (rows of C).
/// * `a` — state matrix A (n×n), overwritten by A'.
/// * `b` — input matrix B (n×m), overwritten by C' (n×p). Must have at least max(m,p) columns.
/// * `c` — output matrix C (p×n), overwritten by B' (m×n). Must have at least max(m,p) rows.
/// * `d` — if JobD::D, direct matrix D (p×m), overwritten by D' (m×p). Must have at least max(m,p) rows and cols.
///
/// # Returns
/// * 0 success, < 0 invalid argument (-i = i-th argument).
pub fn ab07md(
    jobd: JobD,
    n: usize,
    m: usize,
    p: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: Option<&mut DMatrix<f64>>,
) -> i32 {
    let mplim = max(m, p);
    let minmp = min(m, p);

    // Argument checks
    if n > 0 && (a.nrows() != n || a.ncols() != n) {
        return -6;
    }
    if n > 0 && (b.nrows() != n || b.ncols() < mplim) {
        return -8;
    }
    if (n > 0 && (c.nrows() < mplim || c.ncols() != n)) || (n == 0 && c.nrows() < 1) {
        return -10;
    }
    if matches!(jobd, JobD::D) && (d.is_none() || d.as_ref().map(|x| x.nrows() < mplim || x.ncols() < mplim).unwrap_or(true)) {
        return -12;
    }

    if n == 0 && minmp == 0 {
        return 0;
    }

    if n > 0 {
        // Transpose A in place
        for i in 0..n {
            for j in (i + 1)..n {
                let ai = a[(i, j)];
                let aj = a[(j, i)];
                a[(i, j)] = aj;
                a[(j, i)] = ai;
            }
        }

        // Replace B by C' and C by B'
        for j in 0..mplim {
            if j < minmp {
                // Swap column j of B with row j of C
                for i in 0..n {
                    let bi = b[(i, j)];
                    let cj = c[(j, i)];
                    b[(i, j)] = cj;
                    c[(j, i)] = bi;
                }
            } else if j >= p {
                // j > p: copy B(1:n,j) into C(j,1:n)
                for i in 0..n {
                    c[(j, i)] = b[(i, j)];
                }
            } else {
                // j >= m and j < p: copy C(j,1:n) into B(1:n,j)
                for i in 0..n {
                    b[(i, j)] = c[(j, i)];
                }
            }
        }
    }

    if matches!(jobd, JobD::D) && minmp > 0 {
        if let Some(d) = d {
            // D is p×m on entry, output D' is m×p. Copy to tmp then write transpose.
            let mut tmp = DMatrix::zeros(p, m);
            for i in 0..p {
                for j in 0..m {
                    tmp[(i, j)] = d[(i, j)];
                }
            }
            for i in 0..m {
                for j in 0..p {
                    d[(i, j)] = tmp[(j, i)];
                }
            }
        }
    }

    0
}

/// Convenience wrapper for benchmarking: (n, m) -> INFO. Uses JobD::Z, p=m, zero matrices.
#[inline]
pub fn ab07md_nm(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mplim = max(m, p);
    let mut a = DMatrix::zeros(n, n);
    let mut b = DMatrix::zeros(n, mplim);
    let mut c = DMatrix::zeros(mplim, n);
    ab07md(JobD::Z, n, m, p, &mut a, &mut b, &mut c, None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab07md_trivial() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 1);
        let mut c = DMatrix::zeros(1, 0);
        assert_eq!(
            ab07md(JobD::Z, 0, 0, 0, &mut a, &mut b, &mut c, None),
            0
        );
    }

    #[test]
    fn test_ab07md_dual_1x1() {
        // (A,B,C,D) = (1, 1, 1, 0). Dual is (1, 1, 1, 0).
        let mut a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut d = DMatrix::from_row_slice(1, 1, &[0.0]);
        assert_eq!(
            ab07md(
                JobD::D,
                1,
                1,
                1,
                &mut a,
                &mut b,
                &mut c,
                Some(&mut d),
            ),
            0
        );
        assert!((a[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((b[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((c[(0, 0)] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ab07md_transpose_a() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let mut b = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut c = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        assert_eq!(
            ab07md(JobD::Z, 2, 2, 2, &mut a, &mut b, &mut c, None),
            0
        );
        // A should be transposed: [[1,3],[2,4]]
        assert!((a[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((a[(0, 1)] - 3.0).abs() < 1e-10);
        assert!((a[(1, 0)] - 2.0).abs() < 1e-10);
        assert!((a[(1, 1)] - 4.0).abs() < 1e-10);
    }
}
