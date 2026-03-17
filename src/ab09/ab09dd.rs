//! AB09DD — Singular perturbation approximation formulas (SLICOT).
//!
//! Full port: residualization Ar = A11 + A12*(g*I-A22)^{-1}*A21, etc., with g=0 (C) or 1 (D).

use nalgebra::DMatrix;
use std::cmp::max;

use crate::mb04::blas::dgemm;

const ONE: f64 = 1.0;
const ZERO: f64 = 0.0;

/// Full AB09DD: reduced order by SPA formulas. DICO: 'C' continuous (g=0), 'D' discrete (g=1).
/// IWORK length >= 2*(N-NR), DWORK length >= 4*(N-NR). Column-major.
pub fn ab09dd_full(
    dico: u8,
    n: usize,
    m: usize,
    p: usize,
    nr: usize,
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    d: &mut [f64],
    ldd: usize,
    rcond: &mut f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
) -> i32 {
    let discr = dico == b'D' || dico == b'd';
    if dico != b'C' && dico != b'c' && !discr {
        return -1;
    }
    if nr > n {
        return -5;
    }
    if lda < max(1, n) || ldb < max(1, n) || ldc < max(1, p) || ldd < max(1, p) {
        return -7;
    }

    if nr == n {
        *rcond = ONE;
        return 0;
    }

    let k = nr; // 0-based start of A22
    let ns = n - nr;

    // T = -A22 (continuous) or I - A22 (discrete)
    let mut t = DMatrix::from_fn(ns, ns, |i, j| {
        let val = a[(k + i) + (k + j) * lda];
        if discr && i == j {
            ONE - val
        } else {
            -val
        }
    });

    let tnorm = t.iter().map(|x| x.abs()).fold(0.0_f64, |a, b| a + b) / (ns * ns) as f64;
    let lu = t.lu();
    let inv_t = match lu.try_inverse() {
        Some(inv) => inv,
        None => {
            *rcond = ZERO;
            return 1;
        }
    };
    let inv_norm = inv_t.iter().map(|x| x.abs()).fold(0.0_f64, |a, b| a + b) / (ns * ns) as f64;
    *rcond = if tnorm * inv_norm > 0.0 {
        1.0 / (tnorm * inv_norm)
    } else {
        ONE
    };

    let eps = f64::EPSILON;
    if *rcond <= eps {
        return 1;
    }

    // A21 = A(k:nr, 0:nr) in 0-based: rows k..n, cols 0..nr => NR columns, NS rows
    let a21 = DMatrix::from_fn(ns, nr, |i, j| a[(k + i) + j * lda]);
    let x_a21 = match lu.solve(&a21) {
        Some(x) => x,
        None => {
            *rcond = ZERO;
            return 1;
        }
    };
    for j in 0..nr {
        for i in 0..ns {
            a[(k + i) + j * lda] = x_a21[(i, j)];
        }
    }

    // B2 = B(k:n, 0:m)
    let b2 = DMatrix::from_fn(ns, m, |i, j| b[(k + i) + j * ldb]);
    let x_b2 = match lu.solve(&b2) {
        Some(x) => x,
        None => {
            *rcond = ZERO;
            return 1;
        }
    };
    for j in 0..m {
        for i in 0..ns {
            b[(k + i) + j * ldb] = x_b2[(i, j)];
        }
    }

    // Ar = A11 + A12 * X_A21 (compute product into temp then add to A11)
    let mut ar_add = vec![0.0; nr * nr];
    dgemm(
        nr,
        nr,
        ns,
        ONE,
        &a[k * lda..],
        lda,
        &a[k..],
        lda,
        ZERO,
        &mut ar_add,
        nr,
    );
    for j in 0..nr {
        for i in 0..nr {
            a[i + j * lda] += ar_add[i + j * nr];
        }
    }

    // Br = B1 + A12 * X_B2
    let mut br_add = vec![0.0; nr * m];
    dgemm(nr, m, ns, ONE, &a[k * lda..], lda, &b[k..], ldb, ZERO, &mut br_add, nr);
    for j in 0..m {
        for i in 0..nr {
            b[i + j * ldb] += br_add[i + j * nr];
        }
    }

    // Cr = C1 + C2 * X_A21
    let mut cr_add = vec![0.0; p * nr];
    dgemm(
        p,
        nr,
        ns,
        ONE,
        &c[k * ldc..],
        ldc,
        &a[k..],
        lda,
        ZERO,
        &mut cr_add,
        p,
    );
    for j in 0..nr {
        for i in 0..p {
            c[i + j * ldc] += cr_add[i + j * p];
        }
    }

    // Dr = D + C2 * X_B2
    let mut dr_add = vec![0.0; p * m];
    dgemm(p, m, ns, ONE, &c[k * ldc..], ldc, &b[k..], ldb, ZERO, &mut dr_add, p);
    for j in 0..m {
        for i in 0..p {
            d[i + j * ldd] += dr_add[i + j * p];
        }
    }

    0
}

/// Compatibility: (n, m) -> INFO. Uses NR = n (no reduction) so returns 0.
#[inline]
pub fn ab09dd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let nr = n;
    let mut a = vec![0.0; n * n.max(1)];
    let mut b = vec![0.0; n * m.max(1)];
    let mut c = vec![0.0; p * n.max(1)];
    let mut d = vec![0.0; p * m.max(1)];
    let mut rcond = 0.0;
    let mut iwork = vec![0i32; 2 * n.max(1)];
    let mut dwork = vec![0.0; 4 * n.max(1)];
    ab09dd_full(
        b'C',
        n,
        m,
        p,
        nr,
        &mut a,
        n.max(1),
        &mut b,
        n.max(1),
        &mut c,
        p.max(1),
        &mut d,
        p.max(1),
        &mut rcond,
        &mut iwork,
        &mut dwork,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09dd_trivial() {
        assert_eq!(ab09dd(0, 0), 0);
    }

    #[test]
    fn test_ab09dd_nr_eq_n() {
        let n = 2;
        let m = 1;
        let p = 1;
        let nr = 2;
        let mut a = vec![1.0, 0.0, 0.0, 2.0];
        let mut b = vec![1.0, 0.0];
        let mut c = vec![1.0, 0.0];
        let mut d = vec![0.0];
        let mut rcond = -1.0;
        let mut iwork = vec![0i32; 4];
        let mut dwork = vec![0.0; 8];
        let info = ab09dd_full(
            b'C',
            n,
            m,
            p,
            nr,
            &mut a,
            2,
            &mut b,
            2,
            &mut c,
            1,
            &mut d,
            1,
            &mut rcond,
            &mut iwork,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert_eq!(rcond, 1.0);
    }

    #[test]
    fn test_ab09dd_reduce_2_to_1() {
        let n = 2;
        let m = 1;
        let p = 1;
        let nr = 1;
        // A = [-1 0; 0 -2], B = [1; 0], C = [1 0], D = 0. (A22 = -2, g=0 => T = 2)
        let mut a = vec![-1.0, 0.0, 0.0, -2.0];
        let mut b = vec![1.0, 0.0];
        let mut c = vec![1.0, 0.0];
        let mut d = vec![0.0];
        let mut rcond = -1.0;
        let mut iwork = vec![0i32; 2];
        let mut dwork = vec![0.0; 4];
        let info = ab09dd_full(
            b'C',
            n,
            m,
            p,
            nr,
            &mut a,
            2,
            &mut b,
            2,
            &mut c,
            1,
            &mut d,
            1,
            &mut rcond,
            &mut iwork,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert!(rcond > 0.0 && rcond <= 1.0);
        // Ar = A11 + A12*inv(T)*A21 = -1 + 0*inv(2)*0 = -1
        assert!((a[0] - (-1.0)).abs() < 1e-10);
    }
}
