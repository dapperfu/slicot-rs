//! TG01CD — Reduce (A-lambda*E,B) to QR-coordinate form: Q'*E upper trapezoidal (SLICOT TG01CD)

use nalgebra::DMatrix;

/// COMPQ: N = don't compute Q; I = identity then Q; U = Q1*Q.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tg01CdCompq {
    N,
    I,
    U,
}

/// Reduce E to upper trapezoidal via QR: E = Q*R, so Q'*E = R. Apply Q' to A and B.
pub fn tg01cd(
    compq: Tg01CdCompq,
    l: usize,
    n: usize,
    m: usize,
    a: &mut DMatrix<f64>,
    e: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    q: Option<&mut DMatrix<f64>>,
) -> i32 {
    if a.nrows() != l || a.ncols() != n {
        return -5;
    }
    if e.nrows() != l || e.ncols() != n {
        return -7;
    }
    if m > 0 && (b.nrows() != l || b.ncols() != m) {
        return -9;
    }
    if matches!(compq, Tg01CdCompq::I | Tg01CdCompq::U) {
        if let Some(ref qref) = q {
            if qref.nrows() != l || qref.ncols() != l {
                return -11;
            }
        } else {
            return -11;
        }
    }
    if l == 0 {
        return 0;
    }
    let qr = e.clone().qr();
    let q_from_qr = qr.q();
    let r = qr.r();
    for i in 0..l {
        for j in 0..n {
            e[(i, j)] = r[(i, j)];
        }
    }
    let q_tr = q_from_qr.transpose();
    *a = &q_tr * a.clone();
    if m > 0 {
        *b = &q_tr * b.clone();
    }
    if let Some(ref mut qout) = q {
        if compq == Tg01CdCompq::I {
            *qout = q_from_qr;
        } else if compq == Tg01CdCompq::U {
            *qout = qout.clone() * &q_from_qr;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tg01cd_smoke() {
        let l = 4;
        let n = 4;
        let m = 2;
        let a = DMatrix::from_row_slice(l, n, &[-1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 2.0, 1.0, 1.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0]);
        let e = DMatrix::from_row_slice(l, n, &[1.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 3.0, 9.0, 6.0, 3.0, 0.0, 0.0, 2.0, 0.0]);
        let b = DMatrix::from_row_slice(l, m, &[1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        let mut a = a;
        let mut e = e;
        let mut b = b;
        let mut q = DMatrix::zeros(l, l);
        assert_eq!(tg01cd(Tg01CdCompq::I, l, n, m, &mut a, &mut e, &mut b, Some(&mut q)), 0);
        for i in 0..l {
            for j in 0..i {
                assert!(e[(i, j)].abs() < 1e-8, "E should be upper triangular");
            }
        }
    }
}
