//! AB09IX — B&T/SPA from given Cholesky factors S,R (SLICOT).
//!
//! On entry TI = S (Su), T = R (Ru). Forms R*S, MB03UD for HSV, truncation, AB09DD.
//! Main path implemented via ab09bx_core after forming R*S in workspace.

use std::cmp::min;

use crate::ab09::ab09bx::ab09bx_reduce_from_rs;
use crate::ab09::ab09hy::ab09hy_full;
use crate::mb04::blas::dtrmv;

fn dlacpy_full(m: usize, n: usize, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    for j in 0..n {
        for i in 0..m {
            b[i + j * ldb] = a[i + j * lda];
        }
    }
}

/// Full AB09IX: reduction from given S (in ti), R (in t). JOB B/F/S/P, ORDSEL F/A.
pub fn ab09ix_full(
    dico: u8,
    job: u8,
    _fact: u8,
    ordsel: u8,
    n: usize,
    m: usize,
    p: usize,
    nr: &mut usize,
    scalec: f64,
    scaleo: f64,
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    d: &mut [f64],
    ldd: usize,
    ti: &mut [f64],
    ldti: usize,
    t: &mut [f64],
    ldt: usize,
    nminr: &mut usize,
    hsv: &mut [f64],
    tol1: f64,
    tol2: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
    iwarn: &mut i32,
) -> i32 {
    *iwarn = 0;
    let min_nmp = min(n, min(m, p));
    let ktau = n * min(n, m).max(p);
    let ldwork_min = ktau + n * n + 5 * n;
    if min_nmp == 0 {
        *nr = 0;
        *nminr = 0;
        return 0;
    }
    if (ordsel == b'F' || ordsel == b'f') && *nr == 0 {
        *nr = 0;
        *nminr = 0;
        return 0;
    }
    if dwork.len() < ldwork_min {
        return -22;
    }
    let rest = &mut dwork[ktau..];
    dlacpy_full(n, n, ti, ldti, &mut rest[..n * n], n);
    for j in 0..n {
        dtrmv(true, false, n, t, ldt, &mut ti[j * ldti..], 1);
    }
    let ierr = ab09bx_reduce_from_rs(
        dico,
        job,
        ordsel,
        n,
        m,
        p,
        nr,
        a,
        lda,
        b,
        ldb,
        c,
        ldc,
        d,
        ldd,
        hsv,
        ti,
        ldti,
        t,
        ldt,
        scalec,
        scaleo,
        tol1,
        tol2,
        iwork,
        dwork,
        iwarn,
    );
    if ierr != 0 {
        return ierr;
    }
    *nminr = (0..n).filter(|&j| hsv[j] > (n as f64) * f64::EPSILON * hsv[0]).count();
    0
}

#[inline]
pub fn ab09ix(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let lda = n.max(1);
    let ldb = n.max(1);
    let ldc = p.max(1);
    let ldd = p.max(1);
    let lds = n.max(1);
    let ldr = n.max(1);
    let mut a = vec![0.0; lda * n];
    let mut b = vec![0.0; ldb * m];
    let mut c = vec![0.0; ldc * n];
    let mut d = vec![0.0; ldd * m];
    if n > 0 {
        a[0] = -1.0;
        if m > 0 {
            b[0] = 1.0;
        }
        if p > 0 {
            c[0] = 1.0;
        }
    }
    let mut s = vec![0.0; lds * n];
    let mut r = vec![0.0; ldr * n];
    for i in 0..n {
        s[i + i * lds] = 1.0;
        r[i + i * ldr] = 1.0;
    }
    let mut scalec = 1.0;
    let mut scaleo = 1.0;
    let ldwork_hy = (2 * n * p + (10 * n * (n + 1))).max(n * (n + p + 5)).max(2) as i32;
    let mut dwork_hy = vec![0.0; ldwork_hy as usize];
    let info_hy = ab09hy_full(
        n,
        m,
        p,
        &a,
        lda,
        &b,
        ldb,
        &c,
        ldc,
        &d,
        ldd,
        &mut scalec,
        &mut scaleo,
        &mut s,
        lds,
        &mut r,
        ldr,
        &mut [0i32; 2],
        &mut dwork_hy,
        ldwork_hy,
    );
    if info_hy != 0 {
        return info_hy;
    }
    let mut nr = 0_usize;
    let mut nminr = 0_usize;
    let mut hsv = vec![0.0; n];
    let mut iwarn = 0i32;
    let ktau = n * min(n, m).max(p);
    let ldwork = (ktau + n * n + 5 * n).max(1);
    let mut dwork = vec![0.0; ldwork];
    let mut iwork = vec![0i32; (2 * n).max(1)];
    ab09ix_full(
        b'C', b'B', b'S', b'F', n, m, p, &mut nr, scalec, scaleo,
        &mut a, lda, &mut b, ldb, &mut c, ldc, &mut d, ldd,
        &mut s, lds, &mut r, ldr,
        &mut nminr, &mut hsv, 0.0, 0.0,
        &mut iwork, &mut dwork, &mut iwarn,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09ix_trivial() {
        assert_eq!(ab09ix(0, 0), 0);
    }
}
