//! AB09CD — Hankel-norm approximation with D (SLICOT).
//!
//! Full port: TB01ID (optional), TB01WD (Schur), AB09CX.

use nalgebra::DMatrix;
use std::cmp::max;

use crate::ab09::ab09cx::ab09cx_full;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};
use crate::tb01::tb01wd::tb01wd;

/// Full AB09CD: Hankel-norm with D. EQUIL 'S'/'N', DICO, ORDSEL. LDWORK >= max(N*N+3*N, ldwork_cx).
pub fn ab09cd_full(
    dico: u8,
    equil: u8,
    ordsel: u8,
    n: usize,
    m: usize,
    p: usize,
    nr: &mut usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    hsv: &mut [f64],
    tol1: f64,
    tol2: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
    ldwork: i32,
    iwarn: &mut i32,
) -> i32 {
    *iwarn = 0;
    if n == 0 && m == 0 && p == 0 {
        *nr = 0;
        return 0;
    }

    if equil == b'S' || equil == b's' {
        let mut scale = vec![0.0; n];
        let mut maxred = 100.0;
        let info_id = tb01id(Tb01IdJob::All, a, b, c, &mut scale, &mut maxred);
        if info_id != 0 {
            return info_id;
        }
    }

    let mut u = DMatrix::zeros(n, n);
    let mut wr = vec![0.0; n];
    let mut wi = vec![0.0; n];
    let info_wd = tb01wd(a, b, c, &mut u, &mut wr, &mut wi);
    if info_wd != 0 {
        return 1;
    }

    let ldw_cx1 = n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2;
    let ldw_cx2 = n * (m + p + 2) + 2 * m * p + n.min(m) + max(3 * m + 1, n.min(m) + p);
    let ldwork_cx = max(ldw_cx1, ldw_cx2);
    if (ldwork as usize) < ldwork_cx {
        return -20;
    }
    let ierr = ab09cx_full(
        dico,
        ordsel,
        n,
        m,
        p,
        nr,
        a,
        b,
        c,
        d,
        hsv,
        tol1,
        tol2,
        iwork,
        dwork,
        ldwork,
        iwarn,
    );
    if ierr != 0 {
        return ierr + 1;
    }
    if !dwork.is_empty() {
        dwork[0] = (n * n + 3 * n).max(ldwork_cx) as f64;
    }
    0
}

#[inline]
pub fn ab09cd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mut a = DMatrix::zeros(n.max(1), n.max(1));
    let mut b = DMatrix::zeros(n.max(1), m.max(1));
    let mut c = DMatrix::zeros(p.max(1), n.max(1));
    let mut d = DMatrix::zeros(p.max(1), m.max(1));
    let mut nr = 0_usize;
    let mut hsv = vec![0.0; n.max(1)];
    let mut iwarn = 0i32;
    let ldw1 = n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2;
    let ldw2 = n * (m + p + 2) + 2 * m * p + n.min(m) + max(3 * m + 1, n.min(m) + p);
    let ldwork = max(ldw1, ldw2).max(1);
    let mut dwork = vec![0.0; ldwork];
    let mut iwork = vec![0i32; max(1, max(n, m))];
    ab09cd_full(
        b'C',
        b'N',
        b'F',
        n,
        m,
        p,
        &mut nr,
        &mut a,
        &mut b,
        &mut c,
        &mut d,
        &mut hsv,
        0.0,
        0.0,
        &mut iwork,
        &mut dwork,
        ldwork as i32,
        &mut iwarn,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09cd_trivial() {
        assert_eq!(ab09cd(0, 0), 0);
    }
}
