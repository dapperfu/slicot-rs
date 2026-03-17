//! AB09CX — Optimal Hankel-norm approximation (SLICOT).
//!
//! Full port: AB09AX (balanced minimal), order selection, TB01WD or Hankel step (permute, pinv, TB01KD, AB04MD).

use nalgebra::DMatrix;
use std::cmp::{max, min};

use crate::ab04::ab04md::{ab04md, BilinearType};
use crate::ab09::ab09ax::ab09ax_full;
use crate::mb01::mb01sd::{mb01sd, Mb01SdJobs};
use crate::mb04::blas::dgemm;
use crate::tb01::tb01kd::{tb01kd, Dico as Tb01Dico, JobA, StDom};
use crate::tb01::tb01wd::tb01wd;

const ONE: f64 = 1.0;
const ZERO: f64 = 0.0;

fn dlacpy_full(m: usize, n: usize, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    for j in 0..n {
        for i in 0..m {
            b[i + j * ldb] = a[i + j * lda];
        }
    }
}

fn ma02ad_full(m: usize, n: usize, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    for j in 0..n {
        for i in 0..m {
            b[j + i * ldb] = a[i + j * lda];
        }
    }
}

/// Full AB09CX: Hankel-norm approximation. DICO 'C'/'D', ORDSEL 'F'/'A'.
/// LDWORK >= max( N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2, N*(M+P+2) + 2*M*P + min(N,M) + max(3*M+1, min(N,M)+P) ).
pub fn ab09cx_full(
    dico: u8,
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
    let discr = dico == b'D' || dico == b'd';
    let fixord = ordsel == b'F' || ordsel == b'f';

    let min_nmp = min(n, min(m, p));
    if min_nmp == 0 {
        *nr = 0;
        if !iwork.is_empty() {
            iwork[0] = 0;
        }
        if !dwork.is_empty() {
            dwork[0] = ONE;
        }
        return 0;
    }

    let ldw1 = n * (2 * n + max(n, max(m, p)) + 5) + (n * (n + 1)) / 2;
    let ldw2 = n * (m + p + 2) + 2 * m * p + min(n, m) + max(3 * m + 1, min(n, m) + p);
    let ldwork_min = max(ldw1, ldw2);
    if ldwork < 0 || (ldwork as usize) < ldwork_min {
        return -20;
    }
    let liwork = if discr { max(1, max(n, m)) } else { max(1, m) };
    if iwork.len() < liwork {
        return -18;
    }

    let mut nminr = 0_usize;
    let ierr = ab09ax_full(
        dico,
        b'B',
        b'A',
        n,
        m,
        p,
        &mut nminr,
        a,
        b,
        c,
        hsv,
        tol2,
        iwork,
        dwork,
        iwarn,
    );
    if ierr != 0 {
        return ierr;
    }

    let rtol = (n as f64) * f64::EPSILON;
    let mut atol = rtol * hsv[0];
    if fixord {
        if *nr > 0 && *nr > nminr {
            *nr = nminr;
            *iwarn = 1;
        }
    } else {
        atol = if tol1 > 0.0 { tol1.max(atol) } else { atol };
        *nr = 0;
        for i in 0..nminr {
            if hsv[i] <= atol {
                break;
            }
            *nr += 1;
        }
    }

    if *nr == nminr {
        iwork[0] = nminr as i32;
        if !dwork.is_empty() {
            dwork[0] = (n * (n + 2) + 1) as f64;
        }
        let mut u = DMatrix::zeros(nminr, nminr);
        let mut wr = vec![0.0; nminr];
        let mut wi = vec![0.0; nminr];
        let mut a_sub = a.view_range(0..nminr, 0..nminr).into_owned();
        let mut b_sub = b.view_range(0..nminr, 0..m).into_owned();
        let mut c_sub = c.view_range(0..p, 0..nminr).into_owned();
        let ierr_wd = tb01wd(&mut a_sub, &mut b_sub, &mut c_sub, &mut u, &mut wr, &mut wi);
        if ierr_wd != 0 {
            return 3;
        }
        a.view_range_mut(0..nminr, 0..nminr).copy_from(&a_sub);
        b.view_range_mut(0..nminr, 0..m).copy_from(&b_sub);
        c.view_range_mut(0..p, 0..nminr).copy_from(&c_sub);
        return 0;
    }

    let skp = hsv[*nr];
    let srrtol = (rtol).sqrt();
    while *nr > 0 && (hsv[*nr - 1] - skp).abs() <= srrtol * skp {
        *nr -= 1;
    }
    let mut kr = 1_usize;
    for i in (*nr + 2)..nminr {
        if (hsv[i - 1] - skp).abs() > srrtol * skp {
            break;
        }
        kr += 1;
    }
    let nr1 = *nr + 1;
    let nkr1 = min(nminr, nr1 + kr);
    let nu = nminr - *nr - kr;
    let na = *nr + nu;

    if discr {
        let mut a_d = a.view_range(0..nminr, 0..nminr).into_owned();
        let mut b_d = b.view_range(0..nminr, 0..m).into_owned();
        let mut c_d = c.view_range(0..p, 0..nminr).into_owned();
        let mut d_d = d.clone();
        let info_ab04 = ab04md(BilinearType::D2C, ONE, ONE, &mut a_d, &mut b_d, &mut c_d, &mut d_d);
        if info_ab04 != 0 {
            return 1;
        }
        a.view_range_mut(0..nminr, 0..nminr).copy_from(&a_d);
        b.view_range_mut(0..nminr, 0..m).copy_from(&b_d);
        c.view_range_mut(0..p, 0..nminr).copy_from(&c_d);
        *d = d_d;
    }

    let lda = n.max(1);
    let ldb = n.max(1);
    let ldc = p.max(1);
    let ldd = p.max(1);
    let a_sl = a.as_mut_slice();
    let b_sl = b.as_mut_slice();
    let c_sl = c.as_mut_slice();
    let d_sl = d.as_mut_slice();

    if *nr > 0 {
        for j in 0..nu {
            for i in 0..nminr {
                a_sl[i + (nr1 - 1 + j) * lda] = a_sl[i + (nkr1 - 1 + j) * lda];
            }
        }
        for j in 0..na {
            for i in 0..nu {
                a_sl[(nr1 - 1 + i) + j * lda] = a_sl[(nkr1 - 1 + i) + j * lda];
            }
        }
        for j in 0..m {
            for i in 0..nu {
                b_sl[(nr1 - 1 + i) + j * ldb] = b_sl[(nkr1 - 1 + i) + j * ldb];
            }
        }
        for j in 0..nu {
            for i in 0..p {
                c_sl[i + (nr1 - 1 + j) * ldc] = c_sl[i + (nkr1 - 1 + j) * ldc];
            }
        }
    }

    let mut b2_vec = vec![0.0; kr * m];
    let mut c2t_vec = vec![0.0; max(kr, m) * p];
    dlacpy_full(kr, m, &b_sl[(nr1 - 1)..], ldb, &mut b2_vec, kr);
    ma02ad_full(p, kr, &c_sl[(nr1 - 1) * ldc..], ldc, &mut c2t_vec, max(kr, m));

    let b2_mat = DMatrix::from_fn(kr, m, |i, j| b2_vec[i + j * kr]);
    let b2t = b2_mat.transpose();
    let pinv_b2t = match b2t.clone().pseudo_inverse(rtol) {
        Ok(pinv) => pinv,
        Err(_) => return 2,
    };
    let c2t_cols = DMatrix::from_fn(kr, p, |i, j| c2t_vec[i + j * max(kr, m)]);
    let x = &pinv_b2t * &c2t_cols.transpose();
    let mut u_vec = vec![0.0; p * m];
    for j in 0..p {
        for i in 0..m {
            u_vec[j + i * p] = x[(i, j)];
        }
    }

    for j in 0..m {
        for i in 0..p {
            d_sl[i + j * ldd] += skp * u_vec[i + j * p];
        }
    }

    if *nr > 0 {
        let skp2 = skp * skp;
        let mut hsvp = vec![0.0; na];
        for i in 0..*nr {
            hsvp[i] = hsv[i];
        }
        for i in 0..nu {
            hsvp[*nr + i] = hsv[nkr1 - 1 + i];
        }
        let mut hsvp2 = vec![0.0; na];
        for i in 0..na {
            hsvp2[i] = ONE / (hsvp[i] * hsvp[i] - skp2);
        }
        let mut b1_vec = vec![0.0; na * m];
        let mut c1_vec = vec![0.0; p * na];
        dlacpy_full(na, m, b_sl, ldb, &mut b1_vec, na);
        dlacpy_full(p, na, c_sl, ldc, &mut c1_vec, p);

        let mut c_mat = DMatrix::from_fn(p, na, |i, j| c_sl[i + j * ldc]);
        let mut r_p = vec![ONE; p];
        mb01sd(Mb01SdJobs::Column, &mut c_mat, &r_p, &hsvp);
        let mut b1t = vec![0.0; m * na];
        for j in 0..na {
            for k in 0..m {
                b1t[k + j * m] = b1_vec[j + k * na];
            }
        }
        dgemm(p, na, m, -skp, &u_vec, p, &b1t, m, ONE, c_mat.as_mut_slice(), p);
        for i in 0..p {
            for j in 0..na {
                c_sl[i + j * ldc] = c_mat[(i, j)];
            }
        }

        let mut b_mat = DMatrix::from_fn(na, m, |i, j| b_sl[i + j * ldb]);
        let c_na = vec![ONE; m];
        mb01sd(Mb01SdJobs::Row, &mut b_mat, &hsvp, &c_na);
        let mut c1t = vec![0.0; na * p];
        for i in 0..na {
            for j in 0..p {
                c1t[i + j * na] = c1_vec[j + i * p];
            }
        }
        dgemm(na, m, p, -skp, &c1t, na, &u_vec, p, ONE, b_mat.as_mut_slice(), na);
        mb01sd(Mb01SdJobs::Row, &mut b_mat, &hsvp2, &c_na);
        for i in 0..na {
            for j in 0..m {
                b_sl[i + j * ldb] = b_mat[(i, j)];
            }
        }

        let mut a_mat = DMatrix::from_fn(na, na, |i, j| a_sl[i + j * lda]);
        for j in 1..na {
            for i in 0..j {
                let t = a_mat[(i, j)];
                a_mat[(i, j)] = a_mat[(j, i)];
                a_mat[(j, i)] = t;
            }
        }
        dgemm(na, na, m, -ONE, b_mat.as_slice(), na, &b1t, m, -ONE, a_mat.as_mut_slice(), na);
        for i in 0..na {
            for j in 0..na {
                a_sl[i + j * lda] = a_mat[(i, j)];
            }
        }

        let mut a_na = DMatrix::from_fn(na, na, |i, j| a_sl[i + j * lda]);
        let mut b_na = DMatrix::from_fn(na, m, |i, j| b_sl[i + j * ldb]);
        let mut c_na = DMatrix::from_fn(p, na, |i, j| c_sl[i + j * ldc]);
        let mut u_na = DMatrix::zeros(na, na);
        let mut wr_na = vec![0.0; na];
        let mut wi_na = vec![0.0; na];
        let mut ndim = 0_usize;
        let ierr_kd = tb01kd(
            Tb01Dico::Continuous,
            StDom::Stable,
            JobA::General,
            &mut a_na,
            &mut b_na,
            &mut c_na,
            ZERO,
            &mut ndim,
            &mut u_na,
            &mut wr_na,
            &mut wi_na,
        );
        if ierr_kd != 0 {
            return 3;
        }
        if ndim != *nr {
            return 4;
        }
        for i in 0..*nr {
            for j in 0..*nr {
                a_sl[i + j * lda] = a_na[(i, j)];
            }
        }
        for i in 0..*nr {
            for j in 0..m {
                b_sl[i + j * ldb] = b_na[(i, j)];
            }
        }
        for i in 0..p {
            for j in 0..*nr {
                c_sl[i + j * ldc] = c_na[(i, j)];
            }
        }

        if discr {
            let mut a_r = a.view_range(0..*nr, 0..*nr).into_owned();
            let mut b_r = b.view_range(0..*nr, 0..m).into_owned();
            let mut c_r = c.view_range(0..p, 0..*nr).into_owned();
            let mut d_r = d.clone();
            let _ = ab04md(BilinearType::C2D, ONE, ONE, &mut a_r, &mut b_r, &mut c_r, &mut d_r);
            a.view_range_mut(0..*nr, 0..*nr).copy_from(&a_r);
            b.view_range_mut(0..*nr, 0..m).copy_from(&b_r);
            c.view_range_mut(0..p, 0..*nr).copy_from(&c_r);
            *d = d_r;
        }
    }

    iwork[0] = nminr as i32;
    if !dwork.is_empty() {
        dwork[0] = ldwork_min as f64;
    }
    0
}

#[inline]
pub fn ab09cx(n: usize, m: usize) -> i32 {
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
    let ldw2 = n * (m + p + 2) + 2 * m * p + min(n, m) + max(3 * m + 1, min(n, m) + p);
    let ldwork = max(ldw1, ldw2).max(1);
    let mut dwork = vec![0.0; ldwork];
    let mut iwork = vec![0i32; max(1, max(n, m))];
    ab09cx_full(
        b'C',
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
    fn test_ab09cx_trivial() {
        assert_eq!(ab09cx(0, 0), 0);
    }
}
