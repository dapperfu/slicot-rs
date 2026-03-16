//! TB01IZ — Balancing system matrix (A,B,C) complex case (SLICOT TB01IZ)
//!
//! Same as TB01ID for complex A, B, C; SCALE is real.

use nalgebra::DMatrix;
use num_complex::Complex64;

/// Which matrices are involved in balancing (same as TB01ID).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tb01IzJob {
    All,
    B,
    C,
    N,
}

/// Balances the complex system triplet (A,B,C) by real diagonal similarity.
///
/// # Returns
/// 0 on success; < 0 invalid argument.
pub fn tb01iz(
    job: Tb01IzJob,
    a: &mut DMatrix<Complex64>,
    b: &mut DMatrix<Complex64>,
    c: &mut DMatrix<Complex64>,
    scale: &mut [f64],
    maxred: &mut f64,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -6;
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
    let use_b = job == Tb01IzJob::All || job == Tb01IzJob::B;
    let use_c = job == Tb01IzJob::All || job == Tb01IzJob::C;
    for i in 0..n {
        scale[i] = 1.0;
    }
    fn norm_s(n: usize, m: usize, p: usize, a: &DMatrix<Complex64>, b: &DMatrix<Complex64>, c: &DMatrix<Complex64>, use_b: bool, use_c: bool) -> f64 {
        let mut s = 0.0;
        for i in 0..n {
            for j in 0..n {
                s += a[(i, j)].norm();
            }
            if use_b {
                for j in 0..m {
                    s += b[(i, j)].norm();
                }
            }
        }
        if use_c {
            for i in 0..p {
                for j in 0..n {
                    s += c[(i, j)].norm();
                }
            }
        }
        s
    }
    let norm0 = norm_s(n, m, p, a, b, c, use_b, use_c);
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
                row += a[(i, j)].norm();
            }
            if use_b {
                for j in 0..m {
                    row += b[(i, j)].norm();
                }
            }
            let mut col = 0.0;
            for j in 0..n {
                col += a[(j, i)].norm();
            }
            if use_c {
                for j in 0..p {
                    col += c[(j, i)].norm();
                }
            }
            if row == 0.0 && col == 0.0 {
                continue;
            }
            if row == 0.0 || col == 0.0 {
                continue;
            }
            let f = (row / col).sqrt();
            if f < 1.0 / max_red_in || f > max_red_in {
                continue;
            }
            scale[i] *= f;
            changed = true;
            let inv_f = 1.0 / f;
            for j in 0..n {
                a[(i, j)] *= inv_f;
                a[(j, i)] *= f;
            }
            if use_b {
                for j in 0..m {
                    b[(i, j)] *= inv_f;
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
    let norm1 = norm_s(n, m, p, a, b, c, use_b, use_c);
    *maxred = if norm1 > 0.0 { norm0 / norm1 } else { 1.0 };
    0
}
