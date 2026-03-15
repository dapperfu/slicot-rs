//! BB01AD — Benchmark generator for continuous-time algebraic Riccati equations (CAREs).
//!
//! For NR(1)=1, NR(2)=1 returns a single fixed small example (N=2, M=1, P=2).
//! Otherwise returns INFO=1.

/// For nr == [1, 1]: sets N=2, M=1, P=2, fills A, B, C, G, Q, X with fixed small values (column-major),
/// sets vec[0..9] = true, returns 0. Otherwise returns 1.
/// Buffers: a (lda*n), b (ldb*m), c (ldc*n), g (ldg*n), q (ldq*n), x (ldx*n). vec.len() >= 10.
pub fn bb01ad(
    nr: [usize; 2],
    n: &mut usize,
    m: &mut usize,
    p: &mut usize,
    vec: &mut [bool],
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    g: &mut [f64],
    ldg: usize,
    q: &mut [f64],
    ldq: usize,
    x: &mut [f64],
    ldx: usize,
) -> i32 {
    if nr[0] != 1 || nr[1] != 1 {
        return 1;
    }
    *n = 2;
    *m = 1;
    *p = 2;
    if vec.len() < 10 {
        return -1;
    }
    for v in vec.iter_mut().take(10) {
        *v = true;
    }
    let nn = 2;
    if a.len() < lda * nn || b.len() < ldb * 1 || c.len() < ldc * nn || g.len() < ldg * nn
        || q.len() < ldq * nn || x.len() < ldx * nn
    {
        return -2;
    }
    // A 2x2 fixed small values (column-major)
    a[0] = 0.1;
    a[1] = 0.2;
    a[lda] = 0.3;
    a[lda + 1] = 0.4;
    // B 2x1
    b[0] = 0.5;
    b[1] = 0.6;
    // C 2x2
    c[0] = 0.7;
    c[1] = 0.8;
    c[ldc] = 0.9;
    c[ldc + 1] = 1.0;
    // G, Q, X 2x2
    g[0] = 0.1;
    g[1] = 0.0;
    g[ldg] = 0.0;
    g[ldg + 1] = 0.2;
    q[0] = 0.2;
    q[1] = 0.0;
    q[ldq] = 0.0;
    q[ldq + 1] = 0.3;
    x[0] = 0.5;
    x[1] = 0.0;
    x[ldx] = 0.0;
    x[ldx + 1] = 0.5;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bb01ad_nr11() {
        let mut n = 0;
        let mut m = 0;
        let mut p = 0;
        let mut vec = [false; 10];
        let mut a = [0.0_f64; 4];
        let mut b = [0.0_f64; 2];
        let mut c = [0.0_f64; 4];
        let mut g = [0.0_f64; 4];
        let mut q = [0.0_f64; 4];
        let mut x = [0.0_f64; 4];
        let info = bb01ad(
            [1, 1],
            &mut n,
            &mut m,
            &mut p,
            &mut vec,
            &mut a,
            2,
            &mut b,
            2,
            &mut c,
            2,
            &mut g,
            2,
            &mut q,
            2,
            &mut x,
            2,
        );
        assert_eq!(info, 0);
        assert_eq!(n, 2);
        assert_eq!(m, 1);
        assert_eq!(p, 2);
        assert!(vec[0] && vec[1] && vec[2]);
        assert_eq!(a[0], 0.1);
        assert_eq!(b[0], 0.5);
    }

    #[test]
    fn test_bb01ad_other_nr() {
        let mut n = 0;
        let mut m = 0;
        let mut p = 0;
        let mut vec = [false; 10];
        let mut a = [0.0_f64; 4];
        let info = bb01ad(
            [1, 2],
            &mut n,
            &mut m,
            &mut p,
            &mut vec,
            &mut a,
            2,
            &mut [0.0; 2],
            2,
            &mut [0.0; 4],
            2,
            &mut [0.0; 4],
            2,
            &mut [0.0; 4],
            2,
            &mut [0.0; 4],
            2,
        );
        assert_eq!(info, 1);
    }
}
