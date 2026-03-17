//! AB09AD — Balanced truncation model reduction (SLICOT).
//!
//! Full port: TB01ID (optional equilibration), TB01WD (Schur), AB09AX (balance & truncate).

use nalgebra::DMatrix;
use std::cmp::max;

use crate::ab09::ab09ax::ab09ax_full;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};
use crate::tb01::tb01wd::tb01wd;

/// Full AB09AD: balanced truncation (DICO, JOB, EQUIL, ORDSEL, N, M, P, NR, A, B, C, HSV, TOL, ...).
///
/// * DICO: 'C' continuous, 'D' discrete
/// * JOB: 'B' sqrt B&T, 'N' balancing-free
/// * EQUIL: 'S' scale (TB01ID), 'N' no equilibration
/// * ORDSEL: 'F' fixed order, 'A' automatic from TOL
/// * LDWORK >= N*(2*N+MAX(N,M,P)+5)+N*(N+1)/2
/// * Returns 0 success; 1 Schur failed; 2 unstable; 3 HSV failed; <0 invalid argument.
pub fn ab09ad_full(
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
    hsv: &mut [f64],
    tol: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
    ldwork: i32,
    iwarn: &mut i32,
) -> i32 {
    let fixord = ordsel == b'F' || ordsel == b'f';
    *iwarn = 0;

    if dico != b'C' && dico != b'c' && dico != b'D' && dico != b'd' {
        return -1;
    }
    if job != b'B' && job != b'b' && job != b'N' && job != b'n' {
        return -2;
    }
    if equil != b'S' && equil != b's' && equil != b'N' && equil != b'n' {
        return -3;
    }
    if !fixord && ordsel != b'A' && ordsel != b'a' {
        return -4;
    }
    if n == 0 && (m != 0 || p != 0) {
        return -5;
    }
    if fixord && (*nr > n) {
        return -8;
    }
    if n > 0 && (a.nrows() != n || a.ncols() != n) {
        return -10;
    }
    if n > 0 && m > 0 && (b.nrows() != n || b.ncols() != m) {
        return -12;
    }
    if p > 0 && n > 0 && (c.nrows() != p || c.ncols() != n) {
        return -14;
    }
    let ldwork_min = if n > 0 {
        n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2
    } else {
        1
    };
    if ldwork < 0 || (n > 0 && (ldwork as usize) < ldwork_min) {
        return -19;
    }
    if hsv.len() < n {
        return -15;
    }
    let liwork = if job == b'N' || job == b'n' { n } else { 0 };
    if iwork.len() < liwork {
        return -17;
    }

    if n == 0 || (fixord && *nr == 0) {
        *nr = 0;
        if !dwork.is_empty() {
            dwork[0] = 1.0;
        }
        return 0;
    }
    if m == 0 && p == 0 {
        *nr = 0;
        if !dwork.is_empty() {
            dwork[0] = 1.0;
        }
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

    let ierr = ab09ax_full(
        dico,
        job,
        ordsel,
        n,
        m,
        p,
        nr,
        a,
        b,
        c,
        hsv,
        tol,
        iwork,
        dwork,
        iwarn,
    );
    if ierr != 0 {
        return ierr + 1;
    }

    if !dwork.is_empty() {
        let wrkopt = (n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2) as f64;
        dwork[0] = wrkopt;
    }
    0
}

/// Compatibility: (n, m) -> INFO. P = m. Continuous, sqrt B&T, no equil, fixed order 0.
#[inline]
pub fn ab09ad(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mut a = DMatrix::zeros(n.max(1), n.max(1));
    let mut b = DMatrix::zeros(n.max(1), m.max(1));
    let mut c = DMatrix::zeros(p.max(1), n.max(1));
    let mut nr = 0_usize;
    let mut hsv = vec![0.0; n.max(1)];
    let ldwork = (n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2).max(1);
    let mut dwork = vec![0.0; ldwork];
    let mut iwork = vec![0i32; n.max(1)];
    let mut iwarn = 0i32;
    ab09ad_full(
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
        &mut hsv,
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
    fn test_ab09ad_trivial() {
        assert_eq!(ab09ad(0, 0), 0);
    }

    #[test]
    fn test_ab09ad_simple() {
        let n = 1;
        let m = 1;
        let p = 1;
        let mut a = DMatrix::from_row_slice(n, n, &[-1.0]);
        let mut b = DMatrix::from_row_slice(n, m, &[1.0]);
        let mut c = DMatrix::from_row_slice(p, n, &[1.0]);
        let mut nr = 1_usize;
        let mut hsv = vec![0.0; n];
        let ldwork = n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2;
        let mut dwork = vec![0.0; ldwork];
        let mut iwork = vec![0i32; n];
        let mut iwarn = 0i32;
        let info = ab09ad_full(
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
            &mut hsv,
            0.0,
            &mut iwork,
            &mut dwork,
            ldwork as i32,
            &mut iwarn,
        );
        assert_eq!(info, 0);
        assert_eq!(nr, 1);
    }
}
