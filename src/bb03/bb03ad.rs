//! BB03AD — Benchmark generator for Lyapunov equations.
//!
//! For NR(1)=1, NR(2)=1 returns a single fixed small example (N=2, M=1, P=2).
//! Otherwise returns INFO=1.

/// For nr == [1, 1]: sets N=2, M=1, P=2, fills output arrays with fixed small values (column-major),
/// sets vec flags, returns 0. Otherwise returns 1.
pub fn bb03ad(
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
    x: &mut [f64],
    ldx: usize,
) -> i32 {
    if nr[0] != 1 || nr[1] != 1 {
        return 1;
    }
    *n = 2;
    *m = 1;
    *p = 2;
    if vec.len() < 8 {
        return -1;
    }
    for v in vec.iter_mut().take(8) {
        *v = true;
    }
    let nn = 2;
    if a.len() < lda * nn || b.len() < ldb * 1 || c.len() < ldc * nn || x.len() < ldx * nn {
        return -2;
    }
    a[0] = -0.1;
    a[1] = 0.2;
    a[lda] = 0.0;
    a[lda + 1] = -0.3;
    b[0] = 0.5;
    b[1] = 0.6;
    c[0] = 0.7;
    c[1] = 0.8;
    c[ldc] = 0.9;
    c[ldc + 1] = 1.0;
    x[0] = 0.25;
    x[1] = 0.0;
    x[ldx] = 0.0;
    x[ldx + 1] = 0.25;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bb03ad_nr11() {
        let mut n = 0;
        let mut m = 0;
        let mut p = 0;
        let mut vec = [false; 8];
        let mut a = [0.0_f64; 4];
        let mut b = [0.0_f64; 2];
        let mut c = [0.0_f64; 4];
        let mut x = [0.0_f64; 4];
        let info = bb03ad(
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
            &mut x,
            2,
        );
        assert_eq!(info, 0);
        assert_eq!(n, 2);
    }
}
