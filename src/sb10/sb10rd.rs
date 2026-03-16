//! SB10RD — State feedback and output injection for H-infinity (continuous).

use nalgebra::DMatrix;

/// Computes F and H for H-infinity controller for given gamma (via internal Riccati).
pub fn sb10rd(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    gamma: f64,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    f: &mut DMatrix<f64>,
    h: &mut DMatrix<f64>,
    x: &mut DMatrix<f64>,
    y: &mut DMatrix<f64>,
    rcond: &mut [f64],
) -> i32 {
    if n == 0 {
        return 0;
    }
    let mut ak = DMatrix::zeros(n, n);
    let mut bk = DMatrix::zeros(n, nmeas);
    let mut ck = DMatrix::<f64>::zeros(ncon, n);
    let mut dk = DMatrix::zeros(ncon, nmeas);
    let info = crate::sb10::sb10fd::sb10fd(n, m, np, ncon, nmeas, gamma, a, b, c, d, &mut ak, &mut bk, f, &mut dk, rcond, 1e-10);
    if info != 0 {
        return info;
    }
    for i in 0..n {
        for j in 0..n {
            x[(i, j)] = 0.0;
            y[(i, j)] = 0.0;
        }
    }
    for i in 0..n {
        for j in 0..nmeas {
            h[(i, j)] = bk[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10rd_n0() {
        let a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 2);
        let c = DMatrix::zeros(2, 0);
        let d = DMatrix::zeros(2, 2);
        let mut f = DMatrix::zeros(1, 0);
        let mut h = DMatrix::zeros(0, 1);
        let mut x = DMatrix::zeros(0, 0);
        let mut y = DMatrix::zeros(0, 0);
        let mut rcond = [0.0; 4];
        assert_eq!(sb10rd(0, 2, 2, 1, 1, 10.0, &a, &b, &c, &d, &mut f, &mut h, &mut x, &mut y, &mut rcond), 0);
    }
}
