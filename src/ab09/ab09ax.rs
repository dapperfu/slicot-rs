//! AB09AX — Balance & Truncate (square-root or balancing-free). Full port from SLICOT AB09AX.f.

use nalgebra::DMatrix;
use std::cmp::min;

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

/// Balance & Truncate. DICO: 'C'=continuous 'D'=discrete; JOB: 'B'=sqrt B&T 'N'=balancing-free;
/// ORDSEL: 'F'=fixed order 'A'=automatic. All arrays column-major. LDWORK >= N*(MAX(N,M,P)+5) + N*(N+1)/2.
pub fn ab09ax(
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
    hsv: &mut [f64],
    t: &mut [f64],
    ldt: usize,
    ti: &mut [f64],
    ldti: usize,
    tol: f64,
    _iwork: &mut [i32],
    dwork: &mut [f64],
    iwarn: &mut i32,
) -> i32 {
    let discr = dico == b'D' || dico == b'd';
    let bal = job == b'B' || job == b'b';
    let fixord = ordsel == b'F' || ordsel == b'f';

    let min_nmp = min(n, min(m, p));
    let ldwork_min = n * (min(n, m).max(p) + 5) + (n * (n + 1)) / 2;
    if min_nmp == 0 || (fixord && *nr == 0) {
        *nr = 0;
        if !dwork.is_empty() {
            dwork[0] = ONE;
        }
        return 0;
    }
    if dwork.len() < ldwork_min {
        return -22;
    }

    let ku = 0;
    let n_mp = n * min(n, m).max(p);
    let ktau = ku + n_mp;
    let _kw = ktau + n;
    let mut scalec = ONE;
    let mut scaleo = ONE;

    let (buf_b, rest) = dwork.split_at_mut(ktau);
    let (tau_sl, work_sl) = rest.split_at_mut(n);

    // Copy B to DWORK(KU), solve for Su (output in TI)
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

    // Copy C to DWORK(KU), solve for Ru (output in T)
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

    // Ru*Su in TI (column j = T(1:j,1:j) * TI(1:j,j))
    for j in 0..n {
        dtrmv(true, false, j + 1, &t[..], ldt, &mut ti[j * ldti..], 1);
    }

    // SVD of Ru*Su: HSV, V in DWORK(KU), U' in TI
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
        atol = if tol > 0.0 { tol.max(atol) } else { atol };
        *nr = 0;
        for j in 0..n {
            if hsv[j] <= atol {
                break;
            }
            *nr += 1;
        }
    }

    if *nr == 0 {
        dwork[0] = (kw_next + 5 * n) as f64;
        return 0;
    }

    // TI' = Ru'*V1 in U (DWORK(KU))
    dtrmm_left(true, n, *nr, ONE, t, ldt, &mut buf_b[ku..], n);

    // T = (U1)' then T = Su*U1
    ma02ad_full(*nr, n, ti, ldti, t, ldt);
    if packed {
        for j in 0..*nr {
            dtpmv_upper(n, &rest[..], &mut t[j * ldt..], 1);
        }
    } else {
        dtrmm_left(false, n, *nr, ONE, &rest[..], n, t, ldt);
    }

    if bal {
        for j in 0..*nr {
            let temp = ONE / hsv[j].sqrt();
            dscal(n, temp, &mut t[j * ldt..], 1);
            dscal(n, temp, &mut buf_b[ku + j * n..], 1);
        }
    } else {
        let t_mat = DMatrix::from_fn(n, *nr, |i, j| t[i + j * ldt]);
        let qr_t = t_mat.qr();
        let q_t = qr_t.q();
        for j in 0..*nr {
            for i in 0..n {
                t[i + j * ldt] = q_t[(i, j)];
            }
        }
        let u_mat = DMatrix::from_fn(n, *nr, |i, j| buf_b[ku + i + j * n]);
        let qr_u = u_mat.qr();
        let q_u = qr_u.q();
        for j in 0..*nr {
            for i in 0..n {
                buf_b[ku + i + j * n] = q_u[(i, j)];
            }
        }
    }

    // TI := (DWORK(KU))' (N×NR -> NR×N)
    ma02ad_full(n, *nr, &buf_b[ku..], n, ti, ldti);

    if !bal {
        // TI := (TI*T)^{-1} * TI
        let w = DMatrix::from_fn(*nr, *nr, |i, j| {
            (0..n).map(|k| ti[i + k * ldti] * t[k + j * ldt]).sum()
        });
        let mut w_slice = vec![0.0; *nr * *nr];
        for j in 0..*nr {
            for i in 0..*nr {
                w_slice[i + j * *nr] = w[(i, j)];
            }
        }
        let ti_mat = DMatrix::from_fn(*nr, n, |i, j| ti[i + j * ldti]);
        let lu = w.lu();
        if let Some(ti_sol) = lu.solve(&ti_mat) {
            for j in 0..n {
                for i in 0..*nr {
                    ti[i + j * ldti] = ti_sol[(i, j)];
                }
            }
        } else {
            return 2;
        }
    }

    // Ar = TI*A*T
    for j in 0..n {
        let k = min(j + 1, n);
        dgemv(
            false,
            *nr,
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
        *nr,
        *nr,
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

    // Br = TI*B, Cr = C*T
    dlacpy_full(n, m, b, ldb, &mut buf_b[ku..], n);
    dgemm(*nr, m, n, ONE, ti, ldti, &buf_b[ku..], n, ZERO, b, ldb);

    dlacpy_full(p, n, c, ldc, &mut buf_b[ku..], p);
    dgemm(p, *nr, n, ONE, &buf_b[ku..], p, t, ldt, ZERO, c, ldc);

    dwork[0] = (kw_next + 5 * n) as f64;
    0
}

/// Wrapper with DMatrix (copies to/from slices). For fixed order, ordsel='F'.
pub fn ab09ax_full(
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
    hsv: &mut [f64],
    tol: f64,
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
    let ldt = n.max(1);
    let ldti = n.max(1);
    let mut t = vec![0.0; n * ldt];
    let mut ti = vec![0.0; n * ldti];
    let info = ab09ax(
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
        hsv,
        &mut t,
        ldt,
        &mut ti,
        ldti,
        tol,
        iwork,
        dwork,
        iwarn,
    );
    info
}

#[inline]
pub fn ab09ax_simple(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mut a = DMatrix::zeros(n.max(1), n.max(1));
    let mut b = DMatrix::zeros(n.max(1), m.max(1));
    let mut c = DMatrix::zeros(p.max(1), n.max(1));
    let mut nr = 0_usize;
    let mut hsv = vec![0.0; n.max(1)];
    let mut iwarn = 0i32;
    let ldwork = n * (min(n, m).max(p) + 5) + (n * (n + 1)) / 2;
    let mut dwork = vec![0.0; ldwork.max(1)];
    let mut iwork = vec![0i32; if p > 0 { n } else { 0 }];
    ab09ax_full(
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
        &mut hsv,
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
    fn test_ab09ax_trivial() {
        assert_eq!(ab09ax_simple(0, 0), 0);
    }
}
