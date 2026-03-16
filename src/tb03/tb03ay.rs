//! TB03AY — Polynomial matrix representation from state-space (SLICOT TB03AY)
//!
//! Builds polynomial matrix V(s) block-by-block from (A,B,C,D). Uses TB04AD to obtain
//! transfer matrix in row form, then forms P(s) and Q(s) so that T(s) = inv(P)*Q.

use nalgebra::DMatrix;

use crate::tb04::tb04ad::{tb04ad, RowCol};

/// Computes polynomial matrix representation (P(s), Q(s)) from (A,B,C,D) so that T(s) = inv(P)*Q.
///
/// Builds block-by-block using TB04AD for the transfer matrix, then fills P and Q coefficient arrays.
///
/// # Returns
/// 0 success; < 0 invalid argument; > 0 leading coefficient nearly zero (index) or singular matrix.
pub fn tb03ay(
    n: usize,
    m: usize,
    p: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &DMatrix<f64>,
    nr: &mut usize,
    indexp: &mut [i32],
    pcoeff: &mut [f64],
    ldpco1: usize,
    ldpco2: usize,
    qcoeff: &mut [f64],
    ldqco1: usize,
    ldqco2: usize,
    tol: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
) -> i32 {
    if n == 0 {
        *nr = 0;
        for i in 0..p {
            indexp[i] = 0;
        }
        return 0;
    }
    let lddcoe = p;
    let lduco1 = p;
    let lduco2 = m;
    let kdcoef = n + 1;
    if indexp.len() < p {
        return -14;
    }
    if pcoeff.len() < ldpco1 * ldpco2 * kdcoef {
        return -16;
    }
    if qcoeff.len() < ldqco1 * ldqco2 * kdcoef {
        return -20;
    }

    let mut index = vec![0i32; p];
    let mut dcoeff = vec![0.0; lddcoe * kdcoef];
    let mut ucoeff = vec![0.0; lduco1 * lduco2 * kdcoef];

    let info = tb04ad(
        RowCol::R,
        a,
        b,
        c,
        d,
        nr,
        &mut index,
        &mut dcoeff,
        lddcoe,
        &mut ucoeff,
        lduco1,
        lduco2,
        tol,
        tol,
        iwork,
        dwork,
    );
    if info != 0 {
        return info;
    }

    let tol_use = if tol > 0.0 { tol } else { 1e-10 };
    for i in 0..p {
        indexp[i] = index[i];
        if dcoeff[i].abs() < tol_use {
            return (i + 1) as i32;
        }
    }

    for i in 0..p {
        for j in 0..p {
            for k in 0..kdcoef {
                let idx = i + j * ldpco1 + k * ldpco1 * ldpco2;
                if idx < pcoeff.len() {
                    pcoeff[idx] = if i == j && (i + k * lddcoe) < dcoeff.len() {
                        dcoeff[i + k * lddcoe]
                    } else {
                        0.0
                    };
                }
            }
        }
    }
    for i in 0..p {
        for j in 0..m {
            for k in 0..kdcoef {
                let u_idx = i + j * lduco1 + k * lduco1 * lduco2;
                let q_idx = i + j * ldqco1 + k * ldqco1 * ldqco2;
                if q_idx < qcoeff.len() && u_idx < ucoeff.len() {
                    qcoeff[q_idx] = ucoeff[u_idx];
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
    fn tb03ay_n0() {
        let mut a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 1);
        let mut c = DMatrix::<f64>::zeros(1, 0);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut nr = 1;
        let mut indexp = [0i32; 1];
        let mut pcoeff = [0.0; 2];
        let mut qcoeff = [0.0; 2];
        let mut iwork = [0i32; 2];
        let mut dwork = [0.0; 2];
        let info = tb03ay(
            0, 1, 1,
            &mut a, &mut b, &mut c, &d,
            &mut nr, &mut indexp,
            &mut pcoeff, 1, 2,
            &mut qcoeff, 1, 1,
            0.0, &mut iwork, &mut dwork,
        );
        assert_eq!(info, 0);
        assert_eq!(nr, 0);
        assert_eq!(indexp[0], 0);
    }

    #[test]
    fn tb03ay_smoke() {
        let mut a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let mut c = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut nr = 0;
        let mut indexp = [0i32; 1];
        let mut pcoeff = vec![0.0; 1 * 2 * 3];
        let mut qcoeff = vec![0.0; 1 * 1 * 3];
        let mut iwork = vec![0i32; 2 + 1];
        let mut dwork = vec![0.0; 100];
        let info = tb03ay(
            2, 1, 1,
            &mut a, &mut b, &mut c, &d,
            &mut nr, &mut indexp,
            &mut pcoeff, 1, 2,
            &mut qcoeff, 1, 1,
            0.0, &mut iwork, &mut dwork,
        );
        assert_eq!(info, 0);
        assert!(nr >= 1);
    }
}
