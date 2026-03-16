//! TG01AZ — Balancing descriptor pencil (A-lambda*E,B,C), complex case (SLICOT TG01AZ)

use nalgebra::DMatrix;
use num_complex::Complex64;

/// Job for balancing (same as TG01AD).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tg01AzJob {
    All,
    B,
    C,
    N,
}

/// Complex balancing: Dl*A*Dr, Dl*E*Dr, Dl*B, C*Dr. Magnitude = |re|+|im| for thresh.
pub fn tg01az(
    job: Tg01AzJob,
    l: usize,
    n: usize,
    m: usize,
    p: usize,
    thresh: f64,
    a: &mut DMatrix<Complex64>,
    e: &mut DMatrix<Complex64>,
    b: &mut DMatrix<Complex64>,
    c: &mut DMatrix<Complex64>,
    lscale: &mut [f64],
    rscale: &mut [f64],
) -> i32 {
    if a.nrows() != l || a.ncols() != n {
        return -7;
    }
    if e.nrows() != l || e.ncols() != n {
        return -9;
    }
    if m > 0 && (b.nrows() != l || b.ncols() != m) {
        return -11;
    }
    if p > 0 && (c.nrows() != p || c.ncols() != n) {
        return -13;
    }
    if lscale.len() < l || rscale.len() < n {
        return -14;
    }
    for i in 0..l {
        lscale[i] = 1.0;
    }
    for j in 0..n {
        rscale[j] = 1.0;
    }
    if l == 0 || n == 0 {
        return 0;
    }
    fn mag(z: Complex64) -> f64 {
        z.re.abs() + z.im.abs()
    }
    let use_b = job == Tg01AzJob::All || job == Tg01AzJob::B;
    let use_c = job == Tg01AzJob::All || job == Tg01AzJob::C;
    const MAX_ITER: usize = 80;
    for _ in 0..MAX_ITER {
        let mut changed = false;
        for i in 0..l {
            let mut r = 0.0;
            for j in 0..n {
                if mag(a[(i, j)]) > thresh {
                    r += mag(a[(i, j)]);
                }
                if mag(e[(i, j)]) > thresh {
                    r += mag(e[(i, j)]);
                }
            }
            if use_b && m > 0 {
                for j in 0..m {
                    if mag(b[(i, j)]) > thresh {
                        r += mag(b[(i, j)]);
                    }
                }
            }
            if r > thresh {
                let f = 1.0 / r.sqrt();
                changed = true;
                lscale[i] *= f;
                for j in 0..n {
                    a[(i, j)] *= f;
                    e[(i, j)] *= f;
                }
                if use_b && m > 0 {
                    for j in 0..m {
                        b[(i, j)] *= f;
                    }
                }
            }
        }
        for j in 0..n {
            let mut c_norm = 0.0;
            for i in 0..l {
                if mag(a[(i, j)]) > thresh {
                    c_norm += mag(a[(i, j)]);
                }
            }
            if use_c && p > 0 {
                for i in 0..p {
                    if mag(c[(i, j)]) > thresh {
                        c_norm += mag(c[(i, j)]);
                    }
                }
            }
            if c_norm > thresh {
                let f = 1.0 / c_norm.sqrt();
                changed = true;
                rscale[j] *= f;
                for i in 0..l {
                    a[(i, j)] *= f;
                    e[(i, j)] *= f;
                }
                if use_c && p > 0 {
                    for i in 0..p {
                        c[(i, j)] *= f;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tg01az_smoke() {
        let l = 2;
        let n = 2;
        let a = DMatrix::from_fn(l, n, |i, j| Complex64::new((i + j) as f64, 0.0));
        let e = DMatrix::from_fn(l, n, |i, j| Complex64::new(if i == j { 1.0 } else { 0.0 }, 0.0));
        let b = DMatrix::from_fn(l, 1, |_, _| Complex64::new(0.0, 0.0));
        let c = DMatrix::from_fn(1, n, |_, _| Complex64::new(0.0, 0.0));
        let mut a = a;
        let mut e = e;
        let mut b = b;
        let mut c = c;
        let mut lscale = vec![0.0; l];
        let mut rscale = vec![0.0; n];
        assert_eq!(tg01az(Tg01AzJob::All, l, n, 1, 1, 0.0, &mut a, &mut e, &mut b, &mut c, &mut lscale, &mut rscale), 0);
        assert!(lscale[0].is_finite());
    }
}
