//! AB09ND — Singular Perturbation Approximation for alpha-stable part (SLICOT).
//!
//! Full port: TB01ID (optional), TB01KD (Unstable), AB09BX on stable part.

use nalgebra::DMatrix;
use std::cmp::max;

use crate::ab09::ab09bx::ab09bx_core;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};
use crate::tb01::tb01kd::{tb01kd, Dico as Tb01Dico, JobA, StDom};

/// Full AB09ND: SPA for alpha-stable part with D. EQUIL 'S'/'N', DICO, JOB, ORDSEL, ALPHA.
/// LDWORK >= 2*N*N + N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2.
pub fn ab09nd_full(
    dico: u8,
    job: u8,
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
    ns: &mut usize,
    hsv: &mut [f64],
    tol1: f64,
    tol2: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
    ldwork: i32,
    iwarn: &mut i32,
) -> i32 {
    *iwarn = 0;
    if n == 0 {
        *nr = 0;
        *ns = 0;
        return 0;
    }
    if m == 0 && p == 0 {
        *nr = 0;
        *ns = 0;
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
    *ns = n - nu;
    if *ns == 0 {
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

    let mut a_stable = DMatrix::zeros(*ns, *ns);
    let mut b_stable = DMatrix::zeros(*ns, m);
    let mut c_stable = DMatrix::zeros(p, *ns);
    for i in 0..*ns {
        for j in 0..*ns {
            a_stable[(i, j)] = a[(nu1 + i, nu1 + j)];
        }
        for j in 0..m {
            b_stable[(i, j)] = b[(nu1 + i, j)];
        }
    }
    for i in 0..p {
        for j in 0..*ns {
            c_stable[(i, j)] = c[(i, nu1 + j)];
        }
    }

    let nn = n * n;
    let ldwork_bx = (*ns) * (2 * (*ns) + max(*ns, max(m, p)) + 5) + (*ns) * (*ns + 1) / 2;
    let kt = 0;
    let kti = nn;
    let kw = 2 * nn;
    if (ldwork as usize) < kw + ldwork_bx {
        return -22;
    }
    let (t_sl, rest) = dwork.split_at_mut(nn);
    let (ti_sl, bx_work) = rest.split_at_mut(nn);

    let lda = (*ns).max(1);
    let ldb = (*ns).max(1);
    let ldc = p.max(1);
    let ldd = p.max(1);
    let ldt = (*ns).max(1);
    let ldti = (*ns).max(1);
    let ierr = ab09bx_core(
        dico,
        job,
        ordsel,
        *ns,
        m,
        p,
        &mut nra,
        a_stable.as_mut_slice(),
        lda,
        b_stable.as_mut_slice(),
        ldb,
        c_stable.as_mut_slice(),
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
        bx_work,
        iwarn,
    );
    if ierr != 0 {
        return ierr + 1;
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
        dwork[0] = (kw + ldwork_bx) as f64;
    }
    0
}

#[inline]
pub fn ab09nd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mut a = DMatrix::zeros(n.max(1), n.max(1));
    let mut b = DMatrix::zeros(n.max(1), m.max(1));
    let mut c = DMatrix::zeros(p.max(1), n.max(1));
    let mut d = DMatrix::zeros(p.max(1), m.max(1));
    let mut nr = 0_usize;
    let mut ns = 0_usize;
    let mut hsv = vec![0.0; n.max(1)];
    let mut iwarn = 0i32;
    let ldwork = (2 * n * n
        + n * (2 * n + max(n, max(m, p)) + 5)
        + n * (n + 1) / 2)
        .max(1);
    let mut dwork = vec![0.0; ldwork];
    let mut iwork = vec![0i32; (2 * n).max(1)];
    ab09nd_full(
        b'C',
        b'B',
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
        &mut ns,
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
    fn test_ab09nd_trivial() {
        assert_eq!(ab09nd(0, 0), 0);
    }

    #[test]
    fn test_ab09nd_n1() {
        let info = ab09nd(1, 1);
        assert!(info == 0 || info > 0);
    }
}
