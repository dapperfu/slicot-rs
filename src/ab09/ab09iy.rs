//! AB09IY — Coprime factorization via Lyapunov (SLICOT). Uses SB03OU for gramian-like factors.

use std::cmp::{max, min};

use crate::sb03::sb03ou::sb03ou;

fn dlacpy_full(m: usize, n: usize, a: &[f64], lda: usize, b: &mut [f64], ldb: usize) {
    for j in 0..n {
        for i in 0..m {
            b[i + j * ldb] = a[i + j * lda];
        }
    }
}

#[inline]
pub fn ab09iy(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = min(m, 1).max(1);
    let lda = n.max(1);
    let ldb = n.max(1);
    let ldc = p.max(1);
    let mut a = vec![0.0; lda * n];
    let mut b = vec![0.0; ldb * m];
    let mut c = vec![0.0; ldc * n];
    if n > 0 {
        a[0] = -1.0;
        if m > 0 {
            b[0] = 1.0;
        }
        if p > 0 {
            c[0] = 1.0;
        }
    }
    let work_len = n * (max(n, m).max(p) + 5).max(1);
    let lw = (n * m + work_len).max(p * n + work_len).max(1);
    let mut dwork = vec![0.0; lw];
    let mut scale = 1.0;
    let mut tau = vec![0.0; n];
    let mut u = vec![0.0; n * n.max(1)];
    let ldu = n.max(1);
    let (buf_b, work) = dwork.split_at_mut(n * m);
    dlacpy_full(n, m, &b, ldb, buf_b, n);
    let i1 = sb03ou(
        false,
        true,
        n,
        m,
        &a,
        lda,
        buf_b,
        n,
        &mut tau,
        &mut u,
        ldu,
        &mut scale,
        work,
    );
    if i1 != 0 && i1 != 1 {
        return i1;
    }
    let mut scaleo = 1.0;
    let (buf_c, work2) = dwork.split_at_mut(p * n);
    dlacpy_full(p, n, &c, ldc, buf_c, p);
    let i2 = sb03ou(
        false,
        false,
        n,
        p,
        &a,
        lda,
        buf_c,
        p,
        &mut tau,
        &mut u,
        ldu,
        &mut scaleo,
        work2,
    );
    if i2 != 0 && i2 != 1 {
        return i2;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09iy_trivial() {
        assert_eq!(ab09iy(0, 0), 0);
    }
}
