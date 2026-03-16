//! TB01ID — Balancing system matrix (A,B,C) by diagonal similarity (SLICOT TB01ID)
//!
//! Reduces the 1-norm of S = [A B; C 0] (or variants) by inv(D)*A*D, inv(D)*B, C*D.

use nalgebra::DMatrix;

/// Which matrices are involved in balancing.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tb01IdJob {
    /// All matrices (A, B, C) involved.
    All,
    /// B and A involved.
    B,
    /// C and A involved.
    C,
    /// Only A (B and C not involved).
    N,
}

/// Balances the system triplet (A,B,C) by diagonal similarity.
///
/// * `job` - Which matrices to include in row/column norms (All, B, C, N).
/// * `a` - N×N state matrix; overwritten with inv(D)*A*D.
/// * `b` - N×M input matrix; overwritten with inv(D)*B (if M>0 and job includes B).
/// * `c` - P×N output matrix; overwritten with C*D (if P>0 and job includes C).
/// * `scale` - Output length N; scale(j) = D(j) applied to row/column j.
/// * `maxred` - On entry: max allowed reduction per step (<=0 => 10.0). On exit: ratio of original 1-norm of S to balanced 1-norm.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument had an illegal value.
pub fn tb01id(
    job: Tb01IdJob,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    scale: &mut [f64],
    maxred: &mut f64,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n {
        return -6;
    }
    if b.nrows() != n {
        return -8;
    }
    if c.ncols() != n {
        return -10;
    }
    if scale.len() < n {
        return -11;
    }
    if n == 0 {
        *maxred = 1.0;
        return 0;
    }
    let mut max_red_in = *maxred;
    if max_red_in <= 0.0 {
        max_red_in = 10.0;
    }
    let use_b = job == Tb01IdJob::All || job == Tb01IdJob::B;
    let use_c = job == Tb01IdJob::All || job == Tb01IdJob::C;

    // Initial scale = 1
    for i in 0..n {
        scale[i] = 1.0;
    }

    fn one_norm_s(n: usize, m: usize, p: usize, a: &DMatrix<f64>, b: &DMatrix<f64>, c: &DMatrix<f64>, use_b: bool, use_c: bool) -> f64 {
        let mut s = 0.0;
        for i in 0..n {
            for j in 0..n {
                s += a[(i, j)].abs();
            }
            if use_b {
                for j in 0..m {
                    s += b[(i, j)].abs();
                }
            }
        }
        if use_c {
            for i in 0..p {
                for j in 0..n {
                    s += c[(i, j)].abs();
                }
            }
        }
        s
    }

    let norm0 = one_norm_s(n, m, p, a, b, c, use_b, use_c);
    if norm0 == 0.0 {
        *maxred = 1.0;
        return 0;
    }

    const MAX_ITER: usize = 80;
    for _ in 0..MAX_ITER {
        let mut changed = false;
        for i in 0..n {
            let mut row = 0.0;
            for j in 0..n {
                row += a[(i, j)].abs();
            }
            if use_b {
                for j in 0..m {
                    row += b[(i, j)].abs();
                }
            }
            let mut col = 0.0;
            for j in 0..n {
                col += a[(j, i)].abs();
            }
            if use_c {
                for j in 0..p {
                    col += c[(j, i)].abs();
                }
            }
            if row == 0.0 && col == 0.0 {
                continue;
            }
            if row == 0.0 || col == 0.0 {
                continue; // no scaling
            }
            let f = (row / col).sqrt();
            if f < 1.0 {
                if f < 1.0 / max_red_in {
                    continue;
                }
            } else if f > max_red_in {
                continue;
            }
            scale[i] *= f;
            changed = true;
            // Apply f to row i (divide) and column i (multiply)
            for j in 0..n {
                a[(i, j)] /= f;
                a[(j, i)] *= f;
            }
            if use_b {
                for j in 0..m {
                    b[(i, j)] /= f;
                }
            }
            if use_c {
                for j in 0..p {
                    c[(j, i)] *= f;
                }
            }
        }
        if !changed {
            break;
        }
    }

    let norm1 = one_norm_s(n, m, p, a, b, c, use_b, use_c);
    *maxred = if norm1 > 0.0 { norm0 / norm1 } else { 1.0 };
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb01id_n_job() {
        let n = 2;
        let m = 1;
        let p = 1;
        let mut a = DMatrix::from_row_slice(n, n, &[1.0, 1e6, 1e-6, 1.0]);
        let mut b = DMatrix::zeros(n, m);
        let mut c = DMatrix::zeros(p, n);
        let mut scale = vec![0.0; n];
        let mut maxred = 0.0;
        let info = tb01id(Tb01IdJob::N, &mut a, &mut b, &mut c, &mut scale, &mut maxred);
        assert_eq!(info, 0);
        assert!(maxred >= 1.0 || (maxred - 1.0).abs() < 1e-10);
    }
}
