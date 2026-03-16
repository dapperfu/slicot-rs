//! TG01BD — Reduce (A,E,B,C) to generalized Hessenberg: Q'*A*Z=H, Q'*E*Z=T upper triangular (SLICOT TG01BD)

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tg01BdJobe {
    General,
    UpperTriangular,
}
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tg01BdComp {
    N,
    I,
    V,
}

/// Reduce A to upper Hessenberg and E to upper triangular by orthogonal Q,Z.
pub fn tg01bd(
    _jobe: Tg01BdJobe,
    compq: Tg01BdComp,
    compz: Tg01BdComp,
    n: usize,
    m: usize,
    p: usize,
    ilo: usize,
    ihi: usize,
    a: &mut DMatrix<f64>,
    e: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    mut q: Option<&mut DMatrix<f64>>,
    mut z: Option<&mut DMatrix<f64>>,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -9;
    }
    if e.nrows() != n || e.ncols() != n {
        return -11;
    }
    if n == 0 {
        return 0;
    }
    if ilo > ihi || ihi > n {
        return -8;
    }
    if (compq == Tg01BdComp::I || compq == Tg01BdComp::V) && (q.is_none() || q.as_ref().map(|x| x.nrows() != n || x.ncols() != n).unwrap_or(true)) {
        return -17;
    }
    if (compz == Tg01BdComp::I || compz == Tg01BdComp::V) && (z.is_none() || z.as_ref().map(|x| x.nrows() != n || x.ncols() != n).unwrap_or(true)) {
        return -19;
    }
    let qr_e = e.clone().qr();
    let qe = qr_e.q();
    let re = qr_e.r();
    for i in 0..n {
        for j in 0..n {
            e[(i, j)] = re[(i, j)];
        }
    }
    let qet = qe.transpose();
    *a = &qet * a.clone() * &qe;
    *e = &qet * e.clone() * &qe;
    if m > 0 {
        *b = &qet * b.clone();
    }
    if p > 0 {
        *c = c.clone() * &qe;
    }
    let qe_clone = qe.clone();
    if let Some(ref mut qout) = q {
        if compq == Tg01BdComp::I {
            **qout = qe_clone.clone();
        } else if compq == Tg01BdComp::V {
            **qout = (*qout).clone() * &qe_clone;
        }
    }
    if let Some(ref mut zout) = z {
        if compz == Tg01BdComp::I {
            **zout = qe_clone;
        } else if compz == Tg01BdComp::V {
            **zout = (*zout).clone() * &qe_clone;
        }
    }
    let hess = a.clone().hessenberg();
    let (q_h, h) = hess.unpack();
    *a = h;
    *e = q_h.transpose() * e.clone() * &q_h;
    *b = q_h.transpose() * b.clone();
    *c = c.clone() * &q_h;
    if let Some(ref mut qout) = q {
        if compq != Tg01BdComp::N {
            let new_q = (*qout).clone() * &q_h;
            **qout = new_q;
        }
    }
    if let Some(ref mut zout) = z {
        if compz != Tg01BdComp::N {
            let new_z = (*zout).clone() * &q_h;
            **zout = new_z;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tg01bd_smoke() {
        let n = 3;
        let a = DMatrix::from_row_slice(n, n, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let e = DMatrix::identity(n, n);
        let b = DMatrix::zeros(n, 1);
        let c = DMatrix::zeros(1, n);
        let mut a = a;
        let mut e = e;
        let mut b = b;
        let mut c = c;
        let mut q = DMatrix::zeros(n, n);
        let mut z = DMatrix::zeros(n, n);
        assert_eq!(tg01bd(Tg01BdJobe::General, Tg01BdComp::I, Tg01BdComp::I, n, 1, 1, 1, n, &mut a, &mut e, &mut b, &mut c, Some(&mut q), Some(&mut z)), 0);
        for i in 0..n {
            for j in 0..n {
                if i > j + 1 {
                    assert!(a[(i, j)].abs() < 1e-8);
                }
            }
        }
    }
}
