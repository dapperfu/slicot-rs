//! SB10AD — H-infinity optimal controller (continuous-time) using Glover-Doyle formulas.
//! Computes minimal gamma and controller via bisection/scan; uses two CARE.

use nalgebra::DMatrix;

use crate::sb10::sb10fd::sb10fd;
use crate::sb10::sb10ld::sb10ld;

/// JOB: 1=bisection, 2=scan, 3=bisection then scan, 4=suboptimal only.
pub fn sb10ad(
    job: i32,
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    gamma: &mut f64,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    ak: &mut DMatrix<f64>,
    bk: &mut DMatrix<f64>,
    ck: &mut DMatrix<f64>,
    dk: &mut DMatrix<f64>,
    ac: &mut DMatrix<f64>,
    bc: &mut DMatrix<f64>,
    cc: &mut DMatrix<f64>,
    dc: &mut DMatrix<f64>,
    rcond: &mut [f64],
    gtol: f64,
    actol: f64,
    tol: f64,
) -> i32 {
    if job == 4 {
        return sb10fd(n, m, np, ncon, nmeas, *gamma, a, b, c, d, ak, bk, ck, dk, rcond, tol);
    }
    let mut g_lo = 0.0_f64;
    let mut g_hi = *gamma;
    let gtol_use = if gtol > 0.0 { gtol } else { 1e-8 };
    let mut best_info = 0i32;
    for _ in 0..60 {
        let g_mid = (g_lo + g_hi) / 2.0;
        *gamma = g_mid;
        let mut a_c = DMatrix::zeros(2 * n, 2 * n);
        let mut b_c = DMatrix::zeros(2 * n, m - ncon);
        let mut c_c = DMatrix::zeros(np - nmeas, 2 * n);
        let mut d_c = DMatrix::zeros(np - nmeas, m - ncon);
        let info_fd = sb10fd(n, m, np, ncon, nmeas, *gamma, a, b, c, d, ak, bk, ck, dk, rcond, tol);
        if info_fd != 0 {
            g_lo = g_mid;
            best_info = info_fd;
            continue;
        }
        let info_ld = sb10ld(n, m, np, ncon, nmeas, a, b, c, d, ak, bk, ck, dk, &mut a_c, &mut b_c, &mut c_c, &mut d_c);
        if info_ld != 0 {
            g_lo = g_mid;
            continue;
        }
        let schur = match a_c.clone().try_schur(1e-14, 100) {
            Some(s) => s,
            None => {
                g_lo = g_mid;
                continue;
            }
        };
        let eigs = schur.complex_eigenvalues();
        let max_re = eigs.iter().map(|c| c.re).fold(f64::NEG_INFINITY, f64::max);
        if max_re < actol {
            g_hi = g_mid;
            ac.copy_from(&a_c);
            bc.copy_from(&b_c);
            cc.copy_from(&c_c);
            dc.copy_from(&d_c);
            best_info = 0;
        } else {
            g_lo = g_mid;
        }
        if g_hi - g_lo < gtol_use * (1.0 + g_hi) {
            break;
        }
    }
    *gamma = g_hi;
    if best_info != 0 {
        let _ = sb10fd(n, m, np, ncon, nmeas, *gamma, a, b, c, d, ak, bk, ck, dk, rcond, tol);
    }
    best_info
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10ad_n0() {
        let mut gamma = 10.0;
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let d = DMatrix::zeros(2, 2);
        let mut ak = DMatrix::zeros(0, 0);
        let mut bk = DMatrix::zeros(0, 1);
        let mut ck = DMatrix::zeros(1, 0);
        let mut dk = DMatrix::zeros(1, 1);
        let mut ac = DMatrix::zeros(0, 0);
        let mut bc = DMatrix::zeros(0, 1);
        let mut cc = DMatrix::zeros(1, 0);
        let mut dc = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 4];
        assert_eq!(
            sb10ad(4, 0, 2, 2, 1, 1, &mut gamma, &a, &b, &c, &d, &mut ak, &mut bk, &mut ck, &mut dk, &mut ac, &mut bc, &mut cc, &mut dc, &mut rcond, 1e-8, -0.1, 1e-10),
            0
        );
    }
}
