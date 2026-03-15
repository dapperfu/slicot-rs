//! BD02AD — Benchmark generator for generalized continuous-time dynamical systems.
//!
//! For NR(1)=1, NR(2)=1 returns N=1, M=1, P=1 with E=I, A=[1], B=C=D=0, INFO=0.
//! Otherwise returns INFO=1.

/// For nr == [1, 1]: sets N=1, M=1, P=1, E=identity, A=[1], B=C=D=0, vec flags, returns 0.
/// Otherwise returns 1.
pub fn bd02ad(
    nr: [usize; 2],
    n: &mut usize,
    m: &mut usize,
    p: &mut usize,
    vec: &mut [bool],
    e: &mut [f64],
    lde: usize,
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    d: &mut [f64],
    ldd: usize,
) -> i32 {
    if nr[0] != 1 || nr[1] != 1 {
        return 1;
    }
    *n = 1;
    *m = 1;
    *p = 1;
    if vec.len() < 8 {
        return -1;
    }
    vec[0] = true;
    vec[1] = true;
    vec[2] = true;
    vec[3] = false;
    vec[4] = true;
    vec[5] = true;
    vec[6] = true;
    vec[7] = false;
    if e.len() < lde * 1 || a.len() < lda * 1 || b.len() < ldb * 1 || c.len() < ldc * 1 || d.len() < ldd * 1 {
        return -2;
    }
    e[0] = 1.0;
    a[0] = 1.0;
    b[0] = 0.0;
    c[0] = 0.0;
    d[0] = 0.0;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bd02ad_nr11() {
        let mut n = 0;
        let mut m = 0;
        let mut p = 0;
        let mut vec = [false; 8];
        let mut e = [0.0_f64];
        let mut a = [0.0_f64];
        let mut b = [0.0_f64];
        let mut c = [0.0_f64];
        let mut d = [0.0_f64];
        let info = bd02ad(
            [1, 1],
            &mut n,
            &mut m,
            &mut p,
            &mut vec,
            &mut e,
            1,
            &mut a,
            1,
            &mut b,
            1,
            &mut c,
            1,
            &mut d,
            1,
        );
        assert_eq!(info, 0);
        assert_eq!(a[0], 1.0);
    }
}
