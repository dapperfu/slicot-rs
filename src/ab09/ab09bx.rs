//! AB09BX — Singular Perturbation Approximation (SPA). Full port from SLICOT AB09BX.f.
//!
//! Same flow as AB09AX (SB03OU×2, MB03UD, truncation) then AB09DD for (Ar,Br,Cr,Dr) with D.

use nalgebra::DMatrix;
use std::cmp::{max, min};

use crate::ab09::ab09dd::ab09dd_full;
use crate::mb03::mb03ud::{mb03ud, Mb03udJobP, Mb03udJobQ};
use crate::mb04::blas::{dgemm, dgemv, dscal, dtrmm_left, dtrmv, dtpmv_upper};
use crate::sb03::sb03ou::sb03ou;

const ONE: f64 = 1.0;
const ZERO: f64 = 0.0;

fn dlacpy_full(m: usize, n: usize, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    for j in 0..n {
        for i in 0..m {
            b[i + j * ldb] = a[i + j * lda];
        }
    }
}

fn ma02dd_pack_upper(n: usize, a: &[f64], lda: usize, ap: &mut [f64]) {
    let mut idx = 0;
    for j in 0..n {
        for i in 0..=j {
            ap[idx] = a[i + j * lda];
            idx += 1;
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

/// SPA model reduction (slice API). DICO 'C'/'D', JOB 'B'/'N', ORDSEL 'F'/'A'. D and T, TI outputs.
/// LDWORK >= N*(MAX(N,M,P)+5) + N*(N+1)/2. IWORK >= max(1, 2*N).
pub fn ab09bx_core(
    dico: u8,
    job: u8,
    ordsel: u8,
    n: usize,
    m: usize,
    p: usize,
    nr: &mut usize,
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    d: &mut [f64],
    ldd: usize,
    hsv: &mut [f64],
    t: &mut [f64],
    ldt: usize,
    ti: &mut [f64],
    ldti: usize,
    tol1: f64,
    tol2: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
    iwarn: &mut i32,
) -> i32 {
    let discr = dico == b'D' || dico == b'd';
    let bal = job == b'B' || job == b'b';
    let fixord = ordsel == b'F' || ordsel == b'f';

    let min_nmp = min(n, min(m, p));
    let ldwork_min = n * (min(n, m).max(p) + 5) + (n * (n + 1)) / 2;
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
    if fixord && *nr == 0 {
        *nr = 0;
        if !iwork.is_empty() {
            iwork[0] = 0;
        }
        if !dwork.is_empty() {
            dwork[0] = ONE;
        }
        return 0;
    }
    if dwork.len() < ldwork_min {
        return -25;
    }
    if iwork.len() < max(1, 2 * n) {
        return -23;
    }

    let ku = 0;
    let n_mp = n * min(n, m).max(p);
    let ktau = ku + n_mp;
    let mut scalec = ONE;
    let mut scaleo = ONE;

    let (buf_b, rest) = dwork.split_at_mut(ktau);
    let (tau_sl, work_sl) = rest.split_at_mut(n);

    dlacpy_full(n, m, b, ldb, &mut buf_b[ku..], n);
    let mut ierr = sb03ou(
        discr,
        true,
        n,
        m,
        a,
        lda,
        &mut buf_b[ku..],
        n,
        tau_sl,
        ti,
        ldti,
        &mut scalec,
        work_sl,
    );
    if ierr != 0 && ierr != 1 {
        return 1;
    }

    dlacpy_full(p, n, c, ldc, &mut buf_b[..], p);
    ierr = sb03ou(
        discr,
        false,
        n,
        p,
        a,
        lda,
        &mut buf_b[..],
        p,
        tau_sl,
        t,
        ldt,
        &mut scaleo,
        work_sl,
    );
    if ierr != 0 && ierr != 1 {
        return 1;
    }

    let kv = ktau;
    let packed = rest.len() < n * (n + 5);
    let (kw_next, packed_size) = if packed {
        ma02dd_pack_upper(n, ti, ldti, &mut rest[..]);
        (kv + (n * (n + 1)) / 2, (n * (n + 1)) / 2)
    } else {
        for j in 0..n {
            for i in 0..=j {
                rest[i + j * n] = ti[i + j * ldti];
            }
        }
        (kv + n * n, n * n)
    };

    for j in 0..n {
        dtrmv(true, false, j + 1, &t[..], ldt, &mut ti[j * ldti..], 1);
    }

    ierr = mb03ud(
        Mb03udJobQ::Compute,
        Mb03udJobP::Compute,
        n,
        ti,
        ldti,
        &mut buf_b[ku..],
        n,
        hsv,
        &mut rest[packed_size..],
    );
    if ierr != 0 {
        return 2;
    }

    dscal(n, ONE / scalec / scaleo, hsv, 1);

    let rtol = (n as f64) * f64::EPSILON;
    let mut atol = rtol * hsv[0];
    if fixord {
        if *nr > 0 && hsv[*nr - 1] <= atol {
            *nr = 0;
            *iwarn = 1;
        }
    } else {
        atol = if tol1 > 0.0 { tol1.max(atol) } else { atol };
        *nr = 0;
        for j in 0..n {
            if hsv[j] <= atol {
                break;
            }
            *nr += 1;
        }
    }

    if *nr == 0 {
        let mut rcond = 0.0;
        let _ns = 0;
        let ab09dd_dwork_len = 4 * n;
        let (dd_dwork, _) = dwork.split_at_mut(ab09dd_dwork_len.min(dwork.len()));
        let (dd_iwork, _) = iwork.split_at_mut((2 * n).min(iwork.len()));
        let _ = ab09dd_full(
            dico,
            n,
            m,
            p,
            0,
            a,
            lda,
            b,
            ldb,
            c,
            ldc,
            d,
            ldd,
            &mut rcond,
            dd_iwork,
            dd_dwork,
        );
        if !iwork.is_empty() {
            iwork[0] = 0;
        }
        dwork[0] = (kw_next + 5 * n) as f64;
        return 0;
    }

    let nr1 = *nr + 1;
    let mut nminr = *nr;
    let atol2 = if tol2 > 0.0 {
        tol2.max(rtol * hsv[0])
    } else {
        rtol * hsv[0]
    };
    for j in nr1..n {
        if hsv[j] <= atol2 {
            break;
        }
        nminr += 1;
    }
    let ns = nminr - *nr;

    dtrmm_left(true, n, nminr, ONE, t, ldt, &mut buf_b[ku..], n);
    ma02ad_full(nminr, n, ti, ldti, t, ldt);
    if packed {
        for j in 0..nminr {
            dtpmv_upper(n, &rest[..], &mut t[j * ldt..], 1);
        }
    } else {
        dtrmm_left(false, n, nminr, ONE, &rest[..], n, t, ldt);
    }

    if bal {
        for j in 0..*nr {
            let temp = ONE / hsv[j].sqrt();
            dscal(n, temp, &mut t[j * ldt..], 1);
            dscal(n, temp, &mut buf_b[ku + j * n..], 1);
        }
    } else {
        let t1 = DMatrix::from_fn(n, *nr, |i, j| t[i + j * ldt]);
        let qr_t = t1.qr();
        let q_t = qr_t.q();
        for j in 0..*nr {
            for i in 0..n {
                t[i + j * ldt] = q_t[(i, j)];
            }
        }
        let u1 = DMatrix::from_fn(n, *nr, |i, j| buf_b[ku + i + j * n]);
        let qr_u = u1.qr();
        let q_u = qr_u.q();
        for j in 0..*nr {
            for i in 0..n {
                buf_b[ku + i + j * n] = q_u[(i, j)];
            }
        }
        if ns > 0 {
            let t2 = DMatrix::from_fn(n, ns, |i, j| t[i + (nr1 - 1 + j) * ldt]);
            let qr_t2 = t2.qr();
            let q_t2 = qr_t2.q();
            for j in 0..ns {
                for i in 0..n {
                    t[i + (nr1 - 1 + j) * ldt] = q_t2[(i, j)];
                }
            }
            let u2 = DMatrix::from_fn(n, ns, |i, j| buf_b[ku + i + (*nr + j) * n]);
            let qr_u2 = u2.qr();
            let q_u2 = qr_u2.q();
            for j in 0..ns {
                for i in 0..n {
                    buf_b[ku + i + (*nr + j) * n] = q_u2[(i, j)];
                }
            }
        }
    }

    ma02ad_full(n, nminr, &buf_b[ku..], n, ti, ldti);

    if !bal {
        let w1 = DMatrix::from_fn(*nr, *nr, |i, j| {
            (0..n).map(|k| ti[i + k * ldti] * t[k + j * ldt]).sum()
        });
        let lu1 = w1.lu();
        let ti1 = DMatrix::from_fn(*nr, n, |i, j| ti[i + j * ldti]);
        if let Some(sol) = lu1.solve(&ti1) {
            for j in 0..n {
                for i in 0..*nr {
                    ti[i + j * ldti] = sol[(i, j)];
                }
            }
        } else {
            return 2;
        }
        if ns > 0 {
            let w2 = DMatrix::from_fn(ns, ns, |i, j| {
                (0..n).map(|k| ti[nr1 - 1 + i + k * ldti] * t[k + (nr1 - 1 + j) * ldt]).sum()
            });
            let lu2 = w2.lu();
            let ti2 = DMatrix::from_fn(ns, n, |i, j| ti[nr1 - 1 + i + j * ldti]);
            if let Some(sol) = lu2.solve(&ti2) {
                for j in 0..n {
                    for i in 0..ns {
                        ti[nr1 - 1 + i + j * ldti] = sol[(i, j)];
                    }
                }
            } else {
                return 2;
            }
        }
    }

    for j in 0..n {
        let k = min(j + 1, n);
        dgemv(
            false,
            nminr,
            k,
            ONE,
            ti,
            ldti,
            &a[j * lda..],
            1,
            ZERO,
            &mut buf_b[ku + j * n..],
            1,
        );
    }
    dgemm(
        nminr,
        nminr,
        n,
        ONE,
        &buf_b[ku..],
        n,
        t,
        ldt,
        ZERO,
        a,
        lda,
    );

    dlacpy_full(n, m, b, ldb, &mut buf_b[ku..], n);
    dgemm(nminr, m, n, ONE, ti, ldti, &buf_b[ku..], n, ZERO, b, ldb);

    dlacpy_full(p, n, c, ldc, &mut buf_b[ku..], p);
    dgemm(p, nminr, n, ONE, &buf_b[ku..], p, t, ldt, ZERO, c, ldc);

    let mut rcond = 0.0;
    let ab09dd_dwork_len = 4 * n;
    let ab09dd_iwork_len = 2 * n;
    let (dd_dwork, _) = dwork.split_at_mut(ab09dd_dwork_len.min(dwork.len()));
    let (dd_iwork, _) = iwork.split_at_mut(ab09dd_iwork_len.min(iwork.len()));
    ierr = ab09dd_full(
        dico,
        nminr,
        m,
        p,
        *nr,
        a,
        lda,
        b,
        ldb,
        c,
        ldc,
        d,
        ldd,
        &mut rcond,
        dd_iwork,
        dd_dwork,
    );
    if ierr != 0 {
        return 2;
    }

    if !iwork.is_empty() {
        iwork[0] = nminr as i32;
    }
    dwork[0] = (kw_next + 5 * n) as f64;
    0
}

/// Reduction from preformed R*S (in ti) and S in dwork[ktau..ktau+n*n]. Used by AB09IX.
/// Caller must set ti = R*S and dwork[ktau..ktau+n*n] = S. scalec, scaleo are inputs.
pub fn ab09bx_reduce_from_rs(
    dico: u8,
    job: u8,
    ordsel: u8,
    n: usize,
    m: usize,
    p: usize,
    nr: &mut usize,
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    d: &mut [f64],
    ldd: usize,
    hsv: &mut [f64],
    ti: &mut [f64],
    ldti: usize,
    t: &mut [f64],
    ldt: usize,
    scalec: f64,
    scaleo: f64,
    tol1: f64,
    tol2: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
    iwarn: &mut i32,
) -> i32 {
    let bal = job == b'B' || job == b'b';
    let fixord = ordsel == b'F' || ordsel == b'f';
    let min_nmp = min(n, min(m, p));
    let ku = 0;
    let n_mp = n * min(n, m).max(p);
    let ktau = ku + n_mp;
    let packed_size = n * n;
    let ldwork_min = ktau + packed_size + 5 * n;
    if min_nmp == 0 {
        *nr = 0;
        return 0;
    }
    if fixord && *nr == 0 {
        *nr = 0;
        return 0;
    }
    if dwork.len() < ldwork_min {
        return -25;
    }
    let (buf_b, rest) = dwork.split_at_mut(ktau);
    let (rest_s, mb03_work) = rest.split_at_mut(packed_size);
    let kw_next = ktau + packed_size;

    let ierr = mb03ud(
        Mb03udJobQ::Compute,
        Mb03udJobP::Compute,
        n,
        ti,
        ldti,
        &mut buf_b[ku..],
        n,
        hsv,
        mb03_work,
    );
    if ierr != 0 {
        return 2;
    }

    dscal(n, ONE / scalec / scaleo, hsv, 1);

    let rtol = (n as f64) * f64::EPSILON;
    let mut atol = rtol * hsv[0];
    if fixord {
        if *nr > 0 && hsv[*nr - 1] <= atol {
            *nr = 0;
            *iwarn = 1;
        }
    } else {
        atol = if tol1 > 0.0 { tol1.max(atol) } else { atol };
        *nr = 0;
        for j in 0..n {
            if hsv[j] <= atol {
                break;
            }
            *nr += 1;
        }
    }

    if *nr == 0 {
        let mut rcond = 0.0;
        let ab09dd_dwork_len = 4 * n;
        let (dd_dwork, _) = dwork.split_at_mut(ab09dd_dwork_len.min(dwork.len()));
        let (dd_iwork, _) = iwork.split_at_mut((2 * n).min(iwork.len()));
        let _ = ab09dd_full(
            dico,
            n,
            m,
            p,
            0,
            a,
            lda,
            b,
            ldb,
            c,
            ldc,
            d,
            ldd,
            &mut rcond,
            dd_iwork,
            dd_dwork,
        );
        if !iwork.is_empty() {
            iwork[0] = 0;
        }
        dwork[0] = (kw_next + 5 * n) as f64;
        return 0;
    }

    let nr1 = *nr + 1;
    let mut nminr = *nr;
    let atol2 = if tol2 > 0.0 {
        tol2.max(rtol * hsv[0])
    } else {
        rtol * hsv[0]
    };
    for j in nr1..n {
        if hsv[j] <= atol2 {
            break;
        }
        nminr += 1;
    }
    let ns = nminr - *nr;

    dtrmm_left(true, n, nminr, ONE, t, ldt, &mut buf_b[ku..], n);
    ma02ad_full(nminr, n, ti, ldti, t, ldt);
    dtrmm_left(false, n, nminr, ONE, &*rest_s, n, t, ldt);

    if bal {
        for j in 0..*nr {
            let temp = ONE / hsv[j].sqrt();
            dscal(n, temp, &mut t[j * ldt..], 1);
            dscal(n, temp, &mut buf_b[ku + j * n..], 1);
        }
    } else {
        let t1 = DMatrix::from_fn(n, *nr, |i, j| t[i + j * ldt]);
        let qr_t = t1.qr();
        let q_t = qr_t.q();
        for j in 0..*nr {
            for i in 0..n {
                t[i + j * ldt] = q_t[(i, j)];
            }
        }
        let u1 = DMatrix::from_fn(n, *nr, |i, j| buf_b[ku + i + j * n]);
        let qr_u = u1.qr();
        let q_u = qr_u.q();
        for j in 0..*nr {
            for i in 0..n {
                buf_b[ku + i + j * n] = q_u[(i, j)];
            }
        }
        if ns > 0 {
            let t2 = DMatrix::from_fn(n, ns, |i, j| t[i + (nr1 - 1 + j) * ldt]);
            let qr_t2 = t2.qr();
            let q_t2 = qr_t2.q();
            for j in 0..ns {
                for i in 0..n {
                    t[i + (nr1 - 1 + j) * ldt] = q_t2[(i, j)];
                }
            }
            let u2 = DMatrix::from_fn(n, ns, |i, j| buf_b[ku + i + (*nr + j) * n]);
            let qr_u2 = u2.qr();
            let q_u2 = qr_u2.q();
            for j in 0..ns {
                for i in 0..n {
                    buf_b[ku + i + (*nr + j) * n] = q_u2[(i, j)];
                }
            }
        }
    }

    ma02ad_full(n, nminr, &buf_b[ku..], n, ti, ldti);

    if !bal {
        let w1 = DMatrix::from_fn(*nr, *nr, |i, j| {
            (0..n).map(|k| ti[i + k * ldti] * t[k + j * ldt]).sum()
        });
        let lu1 = w1.lu();
        let ti1 = DMatrix::from_fn(*nr, n, |i, j| ti[i + j * ldti]);
        if let Some(sol) = lu1.solve(&ti1) {
            for j in 0..n {
                for i in 0..*nr {
                    ti[i + j * ldti] = sol[(i, j)];
                }
            }
        } else {
            return 2;
        }
        if ns > 0 {
            let w2 = DMatrix::from_fn(ns, ns, |i, j| {
                (0..n)
                    .map(|k| ti[nr1 - 1 + i + k * ldti] * t[k + (nr1 - 1 + j) * ldt])
                    .sum()
            });
            let lu2 = w2.lu();
            let ti2 = DMatrix::from_fn(ns, n, |i, j| ti[nr1 - 1 + i + j * ldti]);
            if let Some(sol) = lu2.solve(&ti2) {
                for j in 0..n {
                    for i in 0..ns {
                        ti[nr1 - 1 + i + j * ldti] = sol[(i, j)];
                    }
                }
            } else {
                return 2;
            }
        }
    }

    for j in 0..n {
        let k = min(j + 1, n);
        dgemv(
            false,
            nminr,
            k,
            ONE,
            ti,
            ldti,
            &a[j * lda..],
            1,
            ZERO,
            &mut buf_b[ku + j * n..],
            1,
        );
    }
    dgemm(
        nminr,
        nminr,
        n,
        ONE,
        &buf_b[ku..],
        n,
        t,
        ldt,
        ZERO,
        a,
        lda,
    );

    dlacpy_full(n, m, b, ldb, &mut buf_b[ku..], n);
    dgemm(nminr, m, n, ONE, ti, ldti, &buf_b[ku..], n, ZERO, b, ldb);

    dlacpy_full(p, n, c, ldc, &mut buf_b[ku..], p);
    dgemm(p, nminr, n, ONE, &buf_b[ku..], p, t, ldt, ZERO, c, ldc);

    let mut rcond = 0.0;
    let ab09dd_dwork_len = 4 * n;
    let ab09dd_iwork_len = 2 * n;
    let (dd_dwork, _) = dwork.split_at_mut(ab09dd_dwork_len.min(dwork.len()));
    let (dd_iwork, _) = iwork.split_at_mut(ab09dd_iwork_len.min(iwork.len()));
    let ierr = ab09dd_full(
        dico,
        nminr,
        m,
        p,
        *nr,
        a,
        lda,
        b,
        ldb,
        c,
        ldc,
        d,
        ldd,
        &mut rcond,
        dd_iwork,
        dd_dwork,
    );
    if ierr != 0 {
        return 2;
    }

    if !iwork.is_empty() {
        iwork[0] = nminr as i32;
    }
    dwork[0] = (kw_next + 5 * n) as f64;
    0
}

/// DMatrix wrapper: fixed order (ordsel='F'). D is P×M.
pub fn ab09bx_full(
    dico: u8,
    job: u8,
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
    iwarn: &mut i32,
) -> i32 {
    if n == 0 && m == 0 && p == 0 {
        *nr = 0;
        return 0;
    }
    let lda = n.max(1);
    let ldb = n.max(1);
    let ldc = p.max(1);
    let ldd = p.max(1);
    let ldt = n.max(1);
    let ldti = n.max(1);
    let mut t = vec![0.0; n * ldt];
    let mut ti = vec![0.0; n * ldti];
    ab09bx_core(
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
        &mut t,
        ldt,
        &mut ti,
        ldti,
        tol1,
        tol2,
        iwork,
        dwork,
        iwarn,
    )
}

/// Compatibility: (n, m) -> INFO. Uses default options.
#[inline]
pub fn ab09bx(n: usize, m: usize) -> i32 {
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
    let ldwork = n * (min(n, m).max(p) + 5) + (n * (n + 1)) / 2;
    let mut dwork = vec![0.0; ldwork.max(1)];
    let mut iwork = vec![0i32; (2 * n).max(1)];
    ab09bx_full(
        b'C',
        b'B',
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
        &mut iwarn,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09bx_trivial() {
        assert_eq!(ab09bx(0, 0), 0);
    }
}
