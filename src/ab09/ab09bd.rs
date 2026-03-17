//! AB09BD — Singular Perturbation Approximation with D (SLICOT).
//!
//! Full port: TB01ID (optional), TB01WD (Schur), AB09BX.

use nalgebra::DMatrix;
use std::cmp::max;

use crate::ab09::ab09bx::ab09bx_core;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};
use crate::tb01::tb01wd::tb01wd;

/// Full AB09BD: SPA with D. EQUIL 'S'/'N', DICO, JOB, ORDSEL. LDWORK >= 2*N*N + N*(MAX(N,M,P)+5) + N*(N+1)/2.
pub fn ab09bd_full(
    dico: u8,
    job: u8,
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
    if n == 0 && m == 0 && p == 0 {
        *nr = 0;
        return 0;
    }
    let nn = n * n;
    let ldwork_min = (2 * nn + n * (max(n, max(m, p)) + 5) + (n * (n + 1)) / 2) as i32;
    if ldwork < ldwork_min {
        return -22;
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

    let lda = n.max(1);
    let ldb = n.max(1);
    let ldc = p.max(1);
    let ldd = p.max(1);
    let ldt = n.max(1);
    let ldti = n.max(1);
    let kt = 0;
    let kti = nn;
    let kw = 2 * nn;
    let (t_sl, rest) = dwork.split_at_mut(nn);
    let (ti_sl, ab09bx_dwork) = rest.split_at_mut(nn);

    let ierr = ab09bx_core(
        dico,
        job,
        ordsel,
        n,
        m,
        p,
        nr,
        a.as_mut_slice(),
        lda,
        b.as_mut_slice(),
        ldb,
        c.as_mut_slice(),
        ldc,
        d.as_mut_slice(),
        ldd,
        hsv,
        t_sl,
        ldt,
        ti_sl,
        ldti,
        tol1,
        tol2,
        iwork,
        ab09bx_dwork,
        iwarn,
    );
    if ierr != 0 {
        return ierr + 1;
    }
    if !dwork.is_empty() {
        dwork[0] = (kw + 5 * n) as f64;
    }
    0
}

#[inline]
pub fn ab09bd(n: usize, m: usize) -> i32 {
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
    let ldwork = 2 * n * n + n * (max(n, max(m, p)) + 5) + (n * (n + 1)) / 2;
    let mut dwork = vec![0.0; ldwork.max(1)];
    let mut iwork = vec![0i32; (2 * n).max(1)];
    ab09bd_full(
        b'C',
        b'B',
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
    fn test_ab09bd_trivial() {
        assert_eq!(ab09bd(0, 0), 0);
    }
}
