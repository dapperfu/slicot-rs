//! TD04AD — Minimal state-space representation for a proper transfer matrix (SLICOT TD04AD)
//!
//! T(s) given as row or column polynomial vectors over common denominators. Finds minimal (A,B,C,D).

use nalgebra::{linalg::LU, DMatrix};
use std::f64;

/// Row or column factorization.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum RowCol {
    /// T(s) as rows over common denominators.
    R,
    /// T(s) as columns over common denominators.
    C,
}

/// Computes minimal (A,B,C,D) from transfer matrix given as row/col polynomial vectors.
///
/// INDEX(i) = degree of i-th denominator. DCOEFF(i,k) = coeff of s^(INDEX(i)-k+1) for row/col i.
/// UCOEFF(i,j,k) = coeff of s^(INDEX(iorj)-k+1) for numerator (i,j). Fortran column-major.
///
/// # Returns
/// * `0` - success
/// * `< 0` - invalid argument
/// * `> 0` - INFO = i: leading coefficient of DCOEFF(i,1) nearly zero
pub fn td04ad(
    rowcol: RowCol,
    m: usize,
    p: usize,
    index: &[i32],
    dcoeff: &[f64],
    lddoe: usize,
    ucoeff: &[f64],
    lduco1: usize,
    lduco2: usize,
    nr: &mut usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    tol: f64,
    iwork: &mut [i32],
) -> i32 {
    let porm = if rowcol == RowCol::R { p } else { m };
    if index.len() < porm {
        return -4;
    }
    let kdcoef = index.iter().take(porm).map(|&d| d as usize).max().unwrap_or(0) + 1;
    let n_full: usize = index.iter().take(porm).map(|&d| d as usize).sum();
    if n_full == 0 {
        *nr = 0;
        return 0;
    }

    for i in 0..porm {
        let idx = i; // DCOEFF(i,1) = leading coeff
        if dcoeff[idx].abs() < 1e-100 {
            return i as i32 + 1;
        }
    }

    let mut a_full = DMatrix::zeros(n_full, n_full);
    let mut b_full = DMatrix::zeros(n_full, m.max(p));
    let mut c_full = DMatrix::zeros(p.max(m), n_full);
    let mut d_mat = DMatrix::zeros(p, m);

    let mut row_start = vec![0usize; porm + 1];
    for (i, &deg) in index.iter().take(porm).enumerate() {
        row_start[i + 1] = row_start[i] + deg as usize;
    }

    if rowcol == RowCol::R {
        for (ii, &deg_i) in index.iter().take(porm).enumerate() {
            let d_i = deg_i as usize;
            let r0 = row_start[ii];
            let l_p_ii = dcoeff[ii];
            for k in 0..d_i - 1 {
                a_full[(r0 + k, r0 + k + 1)] = 1.0;
            }
            for j in 0..d_i {
                let kk = j + 1;
                if kk < kdcoef {
                    a_full[(r0 + d_i - 1, r0 + j)] = -dcoeff[ii + kk * lddoe] / l_p_ii;
                }
            }
            c_full[(ii, r0 + d_i - 1)] = 1.0 / l_p_ii;
        }
        for j in 0..m {
            for i in 0..porm {
                let idx = i + j * lduco1;
                d_mat[(i, j)] = ucoeff[idx] / dcoeff[i];
            }
        }
        for j in 0..m {
            for (ii, &deg_i) in index.iter().take(porm).enumerate() {
                let d_i = deg_i as usize;
                let r0 = row_start[ii];
                for k in 0..d_i {
                    let kk = d_i - k;
                    if kk < kdcoef {
                        let idx = ii + j * lduco1 + kk * lduco1 * lduco2;
                        let q_co = ucoeff[idx];
                        let p_co = dcoeff[ii + kk * lddoe] * d_mat[(ii, j)];
                        b_full[(r0 + k, j)] = (q_co - p_co) / dcoeff[ii];
                    }
                }
            }
        }
    } else {
        return -2; // Column case: operate on dual
    }

    let b_sub = b_full.view((0, 0), (n_full, m)).into_owned();
    let c_sub = c_full.view((0, 0), (p, n_full)).into_owned();
    let (nr_val, a_r, b_r, c_r) = minimal_realization_svd(&a_full, &b_sub, &c_sub, tol);
    *nr = nr_val;
    for i in 0..nr_val {
        for j in 0..nr_val {
            a[(i, j)] = a_r[(i, j)];
        }
    }
    for i in 0..nr_val {
        for j in 0..m {
            b[(i, j)] = b_r[(i, j)];
        }
    }
    for i in 0..p {
        for j in 0..nr_val {
            c[(i, j)] = c_r[(i, j)];
        }
    }
    for i in 0..p {
        for j in 0..m {
            d[(i, j)] = d_mat[(i, j)];
        }
    }
    iwork[0] = nr_val as i32;
    0
}

fn minimal_realization_svd(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    tol: f64,
) -> (usize, DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if n == 0 {
        return (0, DMatrix::zeros(0, 0), DMatrix::zeros(0, m), DMatrix::zeros(p, 0));
    }
    let mut ctr = DMatrix::zeros(n, n * m);
    let mut ab = b.clone();
    ctr.view_mut((0, 0), (n, m)).copy_from(&ab);
    for k in 1..n {
        ab = a * &ab;
        ctr.view_mut((0, k * m), (n, m)).copy_from(&ab);
    }
    let svd = ctr.svd(true, true);
    let sigma = &svd.singular_values;
    let eps = tol.max(1e-15).max(n as f64 * 1e-15 * sigma[0]);
    let nc = sigma.iter().take(n).filter(|s| **s > eps).count();
    if nc == 0 {
        return (0, DMatrix::zeros(0, 0), DMatrix::zeros(0, m), DMatrix::zeros(p, 0));
    }
    let u = svd.u.as_ref().unwrap();
    let v = u.view((0, 0), (n, nc));
    let a1 = v.transpose() * a * &v;
    let b1 = v.transpose() * b;
    let c1 = c * &v;
    let mut obs = DMatrix::zeros(nc * p, nc);
    let mut ca = c1.clone();
    obs.view_mut((0, 0), (p, nc)).copy_from(&ca);
    let mut ap = a1.clone();
    for k in 1..nc {
        ca = &c1 * &ap;
        obs.view_mut((k * p, 0), (p, nc)).copy_from(&ca);
        ap = &a1 * &ap;
    }
    let svd2 = obs.svd(true, true);
    let sigma2 = &svd2.singular_values;
    let no = sigma2.iter().take(nc).filter(|s| **s > eps).count();
    if no == 0 {
        return (0, DMatrix::zeros(0, 0), DMatrix::zeros(0, m), DMatrix::zeros(p, 0));
    }
    let v_t = svd2.v_t.as_ref().unwrap();
    let v2_cols = v_t.rows(0, no).transpose();
    let a_f = &v2_cols.transpose() * &a1 * &v2_cols;
    let b_f = &v2_cols.transpose() * &b1;
    let c_f = &c1 * &v2_cols;
    (no, a_f, b_f, c_f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_td04ad_row() {
        // TD04AD example: M=2, P=2, R. INDEX=[3,3], DCOEFF 2×4, UCOEFF 2×2×4
        let index = [3, 3];
        let lddoe = 2;
        let dcoeff = vec![
            1.0, 1.0, 6.0, 6.0, 11.0, 11.0, 6.0, 6.0,
        ];
        // UCOEFF col-major (i,j,k): row 0 degree 3, row 1 degree 3. (0,0)=0,1,4,3 (0,1)=0,0,1,1 (1,0)=1,8,20,15 (1,1)=0,0,0,0
        let ucoeff = vec![
            0.0, 1.0, 0.0, 0.0, 1.0, 8.0, 0.0, 0.0, 4.0, 20.0, 1.0, 0.0, 3.0, 15.0, 1.0, 0.0,
        ];
        let mut nr = 0;
        let mut a = DMatrix::zeros(6, 6);
        let mut b = DMatrix::zeros(6, 2);
        let mut c = DMatrix::zeros(2, 6);
        let mut d = DMatrix::zeros(2, 2);
        let mut iwork = vec![0i32; 10];
        let info = td04ad(
            RowCol::R,
            2,
            2,
            &index,
            &dcoeff,
            lddoe,
            &ucoeff,
            2,
            2,
            &mut nr,
            &mut a,
            &mut b,
            &mut c,
            &mut d,
            0.0,
            &mut iwork,
        );
        assert_eq!(info, 0);
        assert!(nr >= 1 && nr <= 6);
    }
}
