//! AB09HD — B&T/SPA for systems with unstable part (SLICOT). TB01ID, TB01KD, then AB09HY+IX on stable block.

use nalgebra::DMatrix;
use std::cmp::min;

use crate::ab09::ab09hy::ab09hy_full;
use crate::ab09::ab09ix::ab09ix_full;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};
use crate::tb01::tb01kd::{tb01kd, Dico, JobA, StDom};

#[inline]
pub fn ab09hd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mut a = DMatrix::from_fn(n, n, |i, j| if i == j { -1.0 } else { 0.0 });
    if n > 0 && m > 0 {
        a[(0, 0)] = -1.0;
    }
    let mut b = DMatrix::from_fn(n, m, |i, j| if i == 0 && j == 0 { 1.0 } else { 0.0 });
    let mut c = DMatrix::from_fn(p, n, |i, j| if i == 0 && j == 0 { 1.0 } else { 0.0 });
    let mut scale = vec![1.0; n];
    let mut maxred = 0.0;
    let info_id = tb01id(Tb01IdJob::All, &mut a, &mut b, &mut c, &mut scale, &mut maxred);
    if info_id != 0 {
        return info_id;
    }
    let mut ndim = 0_usize;
    let mut u = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.0 });
    let mut wr = vec![0.0; n];
    let mut wi = vec![0.0; n];
    let info_kd = tb01kd(
        Dico::Continuous,
        StDom::Stable,
        JobA::Schur,
        &mut a,
        &mut b,
        &mut c,
        0.0,
        &mut ndim,
        &mut u,
        &mut wr,
        &mut wi,
    );
    if info_kd != 0 {
        return info_kd;
    }
    if ndim == 0 {
        return 0;
    }
    let lda = ndim.max(1);
    let ldb = ndim.max(1);
    let ldc = p.max(1);
    let ldd = p.max(1);
    let lds = ndim.max(1);
    let ldr = ndim.max(1);
    let mut ar = vec![0.0; lda * ndim];
    let mut br = vec![0.0; ldb * m];
    let mut cr = vec![0.0; ldc * ndim];
    let mut dr = vec![0.0; ldd * m];
    for i in 0..ndim {
        for j in 0..ndim {
            ar[i + j * lda] = a[(i, j)];
        }
    }
    for i in 0..ndim {
        for j in 0..m {
            br[i + j * ldb] = b[(i, j)];
        }
    }
    for i in 0..p {
        for j in 0..ndim {
            cr[i + j * ldc] = c[(i, j)];
        }
    }
    let mut s = vec![0.0; lds * ndim];
    let mut r = vec![0.0; ldr * ndim];
    let mut scalec = 1.0;
    let mut scaleo = 1.0;
    let ldwork_hy = (2 * ndim * p + (10 * ndim * (ndim + 1)))
        .max(ndim * (ndim + p + 5))
        .max(2) as i32;
    let mut dwork_hy = vec![0.0; ldwork_hy as usize];
    let info_hy = ab09hy_full(
        ndim,
        m,
        p,
        &ar,
        lda,
        &br,
        ldb,
        &cr,
        ldc,
        &dr,
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
    let mut hsv = vec![0.0; ndim];
    let mut iwarn = 0i32;
    let ktau = ndim * min(ndim, m).max(p);
    let ldwork = (ktau + ndim * ndim + 5 * ndim).max(1);
    let mut dwork = vec![0.0; ldwork];
    let mut iwork = vec![0i32; (2 * ndim).max(1)];
    ab09ix_full(
        b'C',
        b'B',
        b'S',
        b'F',
        ndim,
        m,
        p,
        &mut nr,
        scalec,
        scaleo,
        &mut ar,
        lda,
        &mut br,
        ldb,
        &mut cr,
        ldc,
        &mut dr,
        ldd,
        &mut s,
        lds,
        &mut r,
        ldr,
        &mut nminr,
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
    fn test_ab09hd_trivial() {
        assert_eq!(ab09hd(0, 0), 0);
    }
}
