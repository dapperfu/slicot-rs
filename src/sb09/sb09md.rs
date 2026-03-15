//! SB09MD — Evaluation of closeness of two multivariable sequences (SLICOT SB09MD)
//!
//! Compares two sequences M1(k) and M2(k) for k = 1..N (each NC×NB), outputs
//! SS (sum of squares), SE (quadratic error), PRE (percentage relative error).

use nalgebra::DMatrix;

/// Compares two multivariable sequences H1 (M1) and H2 (M2), fills SS, SE, PRE.
///
/// - `n`: number of parameters (blocks), N >= 0.
/// - `nc`: rows of each M1(k), M2(k), NC >= 0.
/// - `nb`: columns of each M1(k), M2(k), NB >= 0.
/// - `h1`: NC×N*NB, M1(k) stored in columns (k-1)*NB..k*NB (0-based: k*nb..(k+1)*nb).
/// - `h2`: same layout for M2.
/// - `ss`, `se`, `pre`: output NC×NB matrices (overwritten).
/// - `tol`: tolerance; if tol < EPS then EPS is used.
///
/// Returns INFO: 0 = success, < 0 = invalid argument (-i = i-th argument).
pub fn sb09md(
    n: usize,
    nc: usize,
    nb: usize,
    h1: &DMatrix<f64>,
    h2: &DMatrix<f64>,
    ss: &mut DMatrix<f64>,
    se: &mut DMatrix<f64>,
    pre: &mut DMatrix<f64>,
    tol: f64,
) -> i32 {
    let eps = f64::EPSILON;
    let tol_use = if tol < eps { eps } else { tol };
    let bound = 1.0 / tol_use;

    if h1.nrows() != nc || h1.ncols() != n * nb {
        return -4;
    }
    if h2.nrows() != nc || h2.ncols() != n * nb {
        return -6;
    }
    if ss.nrows() != nc || ss.ncols() != nb {
        return -8;
    }
    if se.nrows() != nc || se.ncols() != nb {
        return -10;
    }
    if pre.nrows() != nc || pre.ncols() != nb {
        return -12;
    }

    if n == 0 || nc == 0 || nb == 0 {
        return 0;
    }

    for i in 0..nc {
        for j in 0..nb {
            let mut sum_ss = 0.0;
            let mut sum_se = 0.0;
            let mut overflow = false;
            for k in 0..n {
                let idx = k * nb + j;
                let m1 = h1[(i, idx)];
                let m2 = h2[(i, idx)];
                let diff = m1 - m2;
                if m1.abs() > bound || diff.abs() > bound {
                    overflow = true;
                    break;
                }
                sum_ss += m1 * m1;
                sum_se += diff * diff;
            }
            let (ss_ij, se_ij, pre_ij) = if overflow {
                (bound, bound, 1.0)
            } else if sum_ss <= tol_use {
                (sum_ss, sum_se, 100.0)
            } else {
                let pre_val = 100.0 * (sum_se / sum_ss).sqrt();
                (sum_ss, sum_se, pre_val)
            };
            ss[(i, j)] = ss_ij;
            se[(i, j)] = se_ij;
            pre[(i, j)] = pre_ij;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb09md_trivial() {
        let nc = 0_usize;
        let nb = 0_usize;
        let n = 0_usize;
        let h1 = DMatrix::<f64>::zeros(nc, n * nb);
        let h2 = DMatrix::<f64>::zeros(nc, n * nb);
        let mut ss = DMatrix::<f64>::zeros(nc, nb);
        let mut se = DMatrix::<f64>::zeros(nc, nb);
        let mut pre = DMatrix::<f64>::zeros(nc, nb);
        assert_eq!(sb09md(n, nc, nb, &h1, &h2, &mut ss, &mut se, &mut pre, 0.0), 0);
    }

    #[test]
    fn test_sb09md_example() {
        // SLICOT SB09MD example: N=2, NC=2, NB=2. H1/H2 stored as NC×(N*NB).
        let h1 = DMatrix::from_row_slice(
            2,
            4,
            &[
                1.3373, 0.1205, 0.6618, -0.3372,
                -0.4062, 1.6120, 0.9299, 0.7429,
            ],
        );
        let h2 = DMatrix::from_row_slice(
            2,
            4,
            &[
                1.1480, -0.1837, 0.8843, -0.4947,
                -0.4616, 1.4674, 0.6028, 0.9524,
            ],
        );
        let mut ss = DMatrix::zeros(2, 2);
        let mut se = DMatrix::zeros(2, 2);
        let mut pre = DMatrix::zeros(2, 2);
        let info = sb09md(2, 2, 2, &h1, &h2, &mut ss, &mut se, &mut pre, 0.0);
        assert_eq!(info, 0);
        // SS, SE >= 0; PRE = 100*sqrt(SE/SS) when SS > tol
        for i in 0..2 {
            for j in 0..2 {
                assert!(ss[(i, j)] >= 0.0);
                assert!(se[(i, j)] >= 0.0);
                assert!(pre[(i, j)] >= 0.0 && pre[(i, j)] <= 100.01);
                if ss[(i, j)] > f64::EPSILON {
                    let expected_pre = 100.0 * (se[(i, j)] / ss[(i, j)]).sqrt();
                    assert!((pre[(i, j)] - expected_pre).abs() < 1e-6);
                }
            }
        }
    }
}
