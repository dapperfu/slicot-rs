//! AB09ED — Hankel-norm for alpha-stable part (SLICOT).
//!
//! Full port: TB01ID (optional), TB01KD (Unstable), AB09CX on stable part.

use nalgebra::DMatrix;
use std::cmp::max;

use crate::ab09::ab09cx::ab09cx_full;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};
use crate::tb01::tb01kd::{tb01kd, Dico as Tb01Dico, JobA, StDom};

/// Full AB09ED: Hankel-norm for alpha-stable part. EQUIL 'S'/'N', DICO, ORDSEL, ALPHA.
pub fn ab09ed_full(
    dico: u8,
    equil: u8,
    ordsel: u8,
    n: usize,
    m: usize,
    p: usize,
    nr: &mut usize,
    alpha: f64,
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

    let tdico = if dico == b'D' || dico == b'd' {
        Tb01Dico::Discrete
    } else {
        Tb01Dico::Continuous
    };
    let mut alpha_work = alpha;
    if tdico == Tb01Dico::Discrete && (alpha - 1.0).abs() < 1e-15 {
        alpha_work = 1.0 - f64::EPSILON.sqrt();
    } else if tdico == Tb01Dico::Continuous && alpha.abs() < 1e-15 {
        alpha_work = -f64::EPSILON.sqrt();
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
    let mut nu = 0_usize;
    let ierr = tb01kd(
        tdico,
        StDom::Unstable,
        JobA::General,
        a,
        b,
        c,
        alpha_work,
        &mut nu,
        &mut u,
        &mut wr,
        &mut wi,
    );
    if ierr != 0 {
        return if ierr == 3 { 2 } else { 1 };
    }
    let ns = n - nu;
    if ns == 0 {
        *nr = nu;
        if !dwork.is_empty() {
            dwork[0] = (n * (n + 2) + 3 * n) as f64;
        }
        return 0;
    }

    let nu1 = nu;
    let mut nra = if ordsel == b'F' || ordsel == b'f' {
        (*nr).saturating_sub(nu)
    } else {
        0
    };
    if ordsel == b'F' && *nr < nu {
        *iwarn = 2;
    }

    let mut a_stable = DMatrix::zeros(ns, ns);
    let mut b_stable = DMatrix::zeros(ns, m);
    let mut c_stable = DMatrix::zeros(p, ns);
    for i in 0..ns {
        for j in 0..ns {
            a_stable[(i, j)] = a[(nu1 + i, nu1 + j)];
        }
        for j in 0..m {
            b_stable[(i, j)] = b[(nu1 + i, j)];
        }
    }
    for i in 0..p {
        for j in 0..ns {
            c_stable[(i, j)] = c[(i, nu1 + j)];
        }
    }

    let ldw_cx1 = ns * (2 * ns + max(ns, max(m, p)) + 5) + (ns * (ns + 1)) / 2;
    let ldw_cx2 = ns * (m + p + 2) + 2 * m * p + ns.min(m) + max(3 * m + 1, ns.min(m) + p);
    let ldwork_cx = max(ldw_cx1, ldw_cx2);
    let kw = n * (n + 2) + 3 * n;
    if (ldwork as usize) < kw + ldwork_cx {
        return -20;
    }
    let cx_len = dwork.len() - kw;
    let ierr = ab09cx_full(
        dico,
        ordsel,
        ns,
        m,
        p,
        &mut nra,
        &mut a_stable,
        &mut b_stable,
        &mut c_stable,
        d,
        hsv,
        tol1,
        tol2,
        iwork,
        &mut dwork[kw..],
        cx_len as i32,
        iwarn,
    );
    if ierr != 0 {
        return ierr + 2;
    }
    *nr = nu + nra;
    for i in 0..nra {
        for j in 0..nra {
            a[(nu1 + i, nu1 + j)] = a_stable[(i, j)];
        }
        for j in 0..m {
            b[(nu1 + i, j)] = b_stable[(i, j)];
        }
    }
    for i in 0..p {
        for j in 0..nra {
            c[(i, nu1 + j)] = c_stable[(i, j)];
        }
    }
    if !dwork.is_empty() {
        dwork[0] = (kw + ldwork_cx) as f64;
    }
    0
}

#[inline]
pub fn ab09ed(n: usize, m: usize) -> i32 {
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
    let kw = n * (n + 2) + 3 * n;
    let ldw_cx1 = n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2;
    let ldw_cx2 = n * (m + p + 2) + 2 * m * p + n.min(m) + max(3 * m + 1, n.min(m) + p);
    let ldwork = (kw + max(ldw_cx1, ldw_cx2)).max(1);
    let mut dwork = vec![0.0; ldwork];
    let mut iwork = vec![0i32; max(1, max(n, m))];
    ab09ed_full(
        b'C',
        b'N',
        b'F',
        n,
        m,
        p,
        &mut nr,
        0.0,
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
    fn test_ab09ed_trivial() {
        assert_eq!(ab09ed(0, 0), 0);
    }
}
