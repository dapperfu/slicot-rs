//! TG01AD — Balancing descriptor system pencil (A-lambda*E,B,C) (SLICOT TG01AD)
//!
//! Applies diagonal Dl and Dr so that diag(Dl,I)*S*diag(Dr,I) has balanced row/column 1-norms.

use nalgebra::DMatrix;

/// Which matrices are involved in balancing.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tg01AdJob {
    /// All matrices (A, E, B, C).
    All,
    /// B, A, E.
    B,
    /// C, A, E.
    C,
    /// B and C not involved (A, E only).
    N,
}

/// Balances (A,E,B,C): A,E are L×N, B is L×M, C is P×N.
/// On exit: A <- Dl*A*Dr, E <- Dl*E*Dr, B <- Dl*B, C <- C*Dr; LSCALE(1..L), RSCALE(1..N).
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tg01ad(
    job: Tg01AdJob,
    l: usize,
    n: usize,
    m: usize,
    p: usize,
    thresh: f64,
    a: &mut DMatrix<f64>,
    e: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
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
    if l == 0 || n == 0 {
        for i in 0..l {
            lscale[i] = 1.0;
        }
        for j in 0..n {
            rscale[j] = 1.0;
        }
        return 0;
    }

    let use_b = job == Tg01AdJob::All || job == Tg01AdJob::B;
    let use_c = job == Tg01AdJob::All || job == Tg01AdJob::C;

    for i in 0..l {
        lscale[i] = 1.0;
    }
    for j in 0..n {
        rscale[j] = 1.0;
    }

    fn row_norm_l(
        l: usize,
        n: usize,
        m: usize,
        a: &DMatrix<f64>,
        e: &DMatrix<f64>,
        b: &DMatrix<f64>,
        i: usize,
        thresh: f64,
        use_b: bool,
    ) -> f64 {
        let mut s = 0.0;
        for j in 0..n {
            let v = a[(i, j)].abs();
            if v > thresh {
                s += v;
            }
        }
        for j in 0..n {
            let v = e[(i, j)].abs();
            if v > thresh {
                s += v;
            }
        }
        if use_b && m > 0 {
            for j in 0..m {
                let v = b[(i, j)].abs();
                if v > thresh {
                    s += v;
                }
            }
        }
        s
    }

    fn col_norm_n(
        l: usize,
        n: usize,
        p: usize,
        a: &DMatrix<f64>,
        c: &DMatrix<f64>,
        j: usize,
        thresh: f64,
        use_c: bool,
    ) -> f64 {
        let mut s = 0.0;
        for i in 0..l {
            let v = a[(i, j)].abs();
            if v > thresh {
                s += v;
            }
        }
        if use_c && p > 0 {
            for i in 0..p {
                let v = c[(i, j)].abs();
                if v > thresh {
                    s += v;
                }
            }
        }
        s
    }

    const MAX_ITER: usize = 80;
    for _ in 0..MAX_ITER {
        let mut changed = false;
        for i in 0..l {
            let r = row_norm_l(l, n, m, a, e, b, i, thresh, use_b);
            if r > thresh {
                let f = 1.0 / r.sqrt();
                if (lscale[i] * f - lscale[i]).abs() > 1e-14 * lscale[i].abs() {
                    changed = true;
                }
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
            let c_norm = col_norm_n(l, n, p, a, c, j, thresh, use_c);
            if c_norm > thresh {
                let f = 1.0 / c_norm.sqrt();
                if (rscale[j] * f - rscale[j]).abs() > 1e-14 * rscale[j].abs() {
                    changed = true;
                }
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
    fn test_tg01ad_slicot_example() {
        let l = 4;
        let n = 4;
        let m = 2;
        let p = 2;
        let a = DMatrix::from_row_slice(l, n, &[
            -1.0, 0.0, 0.0, 0.003,
            0.0, 0.0, 0.1, 0.02,
            100.0, 10.0, 0.0, 0.4,
            0.0, 0.0, 0.0, 0.0,
        ]);
        let e = DMatrix::from_row_slice(l, n, &[
            1.0, 0.2, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.01,
            300.0, 90.0, 6.0, 0.3,
            0.0, 0.0, 20.0, 0.0,
        ]);
        let b = DMatrix::from_row_slice(l, m, &[
            0.0, 0.0,
            0.0, 0.0,
            0.0, 1000.0,
            10000.0, 10000.0,
        ]);
        let c = DMatrix::from_row_slice(p, n, &[
            -0.1, 0.0, 0.001, 0.0,
            0.0, 0.01, -0.001, 0.0001,
        ]);
        let mut a = a;
        let mut e = e;
        let mut b = b;
        let mut c = c;
        let mut lscale = vec![0.0; l];
        let mut rscale = vec![0.0; n];
        assert_eq!(tg01ad(
            Tg01AdJob::All, l, n, m, p, 0.0,
            &mut a, &mut e, &mut b, &mut c,
            &mut lscale, &mut rscale,
        ), 0);
        assert!((a[(0, 0)] - (-1.0)).abs() < 1.0, "a[(0,0)] = {}", a[(0, 0)]);
        assert!(lscale[0].abs() < 1e10 && rscale[0].abs() < 1e10);
    }
}
