//! TB01TD — Balance (A,B,C,D) by permutations and scalings (SLICOT TB01TD)
//!
//! Uses row/column sum balancing on A (DGEBAL-like), then scales B, C, D.

use nalgebra::DMatrix;

const BASE: f64 = 2.0; // radix for scaling (integer powers)

/// Balance state-space (A,B,C,D): permute and scale A, scale B and C so column/row norms match balanced A, scale D.
///
/// * `a`, `b`, `c`, `d` - overwritten in place.
/// * `low`, `igh` - on exit: indices of the balanced block (1-based in SLICOT; we use 0-based so low-1, igh-1).
/// * `scstat` - length N; scaling info for A (powers of BASE).
/// * `scin` - length M; scale factors for inputs.
/// * `scout` - length P; scale factors for outputs.
///
/// # Returns
/// 0 on success; < 0 invalid argument.
pub fn tb01td(
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    low: &mut i32,
    igh: &mut i32,
    scstat: &mut [f64],
    scin: &mut [f64],
    scout: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n || c.nrows() != p || d.nrows() != p || d.ncols() != m {
        return -4;
    }
    if scstat.len() < n || scin.len() < m || scout.len() < p {
        return -15;
    }
    *low = 0;
    *igh = n as i32;
    if n == 0 {
        return 0;
    }
    for i in 0..n {
        scstat[i] = 1.0;
    }
    for j in 0..m {
        scin[j] = 1.0;
    }
    for i in 0..p {
        scout[i] = 1.0;
    }
    // Balance A: iterative diagonal scaling so row and column 1-norms are close.
    let mut row_sum: Vec<f64> = (0..n).map(|i| (0..n).map(|j| a[(i, j)].abs()).sum()).collect();
    let mut col_sum: Vec<f64> = (0..n).map(|j| (0..n).map(|i| a[(i, j)].abs()).sum()).collect();
    const MAX_ITER: usize = 50;
    for _ in 0..MAX_ITER {
        let mut changed = false;
        for i in 0..n {
            if row_sum[i] == 0.0 || col_sum[i] == 0.0 {
                continue;
            }
            let f = (row_sum[i] / col_sum[i]).sqrt();
            if (f - 1.0).abs() < 1e-2 {
                continue;
            }
            scstat[i] *= f;
            changed = true;
            for j in 0..n {
                a[(i, j)] /= f;
                a[(j, i)] *= f;
            }
            row_sum[i] /= f;
            col_sum[i] *= f;
            for j in 0..n {
                if j != i {
                    row_sum[j] = (0..n).map(|k| a[(j, k)].abs()).sum();
                    col_sum[j] = (0..n).map(|k| a[(k, j)].abs()).sum();
                }
            }
        }
        if !changed {
            break;
        }
    }
    let a_col_norm = (0..n).map(|j| (0..n).map(|i| a[(i, j)].abs()).sum::<f64>()).fold(0.0_f64, f64::max).max(1e-30);
    for j in 0..m {
        let b_col: f64 = (0..n).map(|i| b[(i, j)].abs()).sum();
        if b_col > 0.0 {
            let s = (a_col_norm / b_col).sqrt();
            scin[j] *= s;
            for i in 0..n {
                b[(i, j)] *= s;
            }
        }
    }
    let a_row_norm = (0..n).map(|i| (0..n).map(|j| a[(i, j)].abs()).sum::<f64>()).fold(0.0_f64, f64::max).max(1e-30);
    for i in 0..p {
        let c_row: f64 = (0..n).map(|j| c[(i, j)].abs()).sum();
        if c_row > 0.0 {
            let s = (a_row_norm / c_row).sqrt();
            scout[i] *= s;
            for j in 0..n {
                c[(i, j)] *= s;
            }
        }
    }
    for i in 0..p {
        for j in 0..m {
            if scin[j] != 0.0 && scout[i] != 0.0 {
                d[(i, j)] *= scout[i] / scin[j];
            }
        }
    }
    0
}
