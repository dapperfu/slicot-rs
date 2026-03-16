//! TB01TY — Balance (A,B,C,D) complex case (SLICOT TB01TY)
//!
//! Complex version of TB01TD: permutations and scalings for (A,B,C,D).

use nalgebra::DMatrix;
use num_complex::Complex64;

const BASE: f64 = 2.0;

/// Balance complex state-space (A,B,C,D). Outputs low, igh, scstat, scin, scout.
///
/// # Returns
/// 0 on success; < 0 invalid argument.
pub fn tb01ty(
    a: &mut DMatrix<Complex64>,
    b: &mut DMatrix<Complex64>,
    c: &mut DMatrix<Complex64>,
    d: &mut DMatrix<Complex64>,
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
    for i in 0..n {
        scstat[i] = 1.0;
    }
    for j in 0..m {
        scin[j] = 1.0;
    }
    for i in 0..p {
        scout[i] = 1.0;
    }
    if n == 0 {
        return 0;
    }
    let mut row_sum: Vec<f64> = (0..n).map(|i| (0..n).map(|j| a[(i, j)].norm()).sum()).collect();
    let mut col_sum: Vec<f64> = (0..n).map(|j| (0..n).map(|i| a[(i, j)].norm()).sum()).collect();
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
            let inv_f = 1.0 / f;
            for j in 0..n {
                a[(i, j)] *= inv_f;
                a[(j, i)] *= f;
            }
            row_sum[i] *= inv_f;
            col_sum[i] *= f;
            for j in 0..n {
                if j != i {
                    row_sum[j] = (0..n).map(|k| a[(j, k)].norm()).sum();
                    col_sum[j] = (0..n).map(|k| a[(k, j)].norm()).sum();
                }
            }
        }
        if !changed {
            break;
        }
    }
    let a_col_norm = (0..n).map(|j| (0..n).map(|i| a[(i, j)].norm()).sum::<f64>()).fold(0.0_f64, f64::max).max(1e-30);
    for j in 0..m {
        let b_col: f64 = (0..n).map(|i| b[(i, j)].norm()).sum();
        if b_col > 0.0 {
            let s = (a_col_norm / b_col).sqrt();
            scin[j] *= s;
            for i in 0..n {
                b[(i, j)] *= s;
            }
        }
    }
    let a_row_norm = (0..n).map(|i| (0..n).map(|j| a[(i, j)].norm()).sum::<f64>()).fold(0.0_f64, f64::max).max(1e-30);
    for i in 0..p {
        let c_row: f64 = (0..n).map(|j| c[(i, j)].norm()).sum();
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
                d[(i, j)] *= Complex64::new(scout[i] / scin[j], 0.0);
            }
        }
    }
    0
}
