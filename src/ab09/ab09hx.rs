//! AB09HX — B&T/SPA with standard gramians (SLICOT). Calls AB09HY then AB09IX.

use std::cmp::min;

use crate::ab09::ab09hy::ab09hy_full;
use crate::ab09::ab09ix::ab09ix_full;

#[inline]
pub fn ab09hx(n: usize, m: usize) -> i32 {
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
        b'C',
        b'B',
        b'S',
        b'F',
        n,
        m,
        p,
        &mut nr,
        scalec,
        scaleo,
        &mut a,
        lda,
        &mut b,
        ldb,
        &mut c,
        ldc,
        &mut d,
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
    fn test_ab09hx_trivial() {
        assert_eq!(ab09hx(0, 0), 0);
    }
}
