//! TF01RD — Markov parameters M(k) = C*A^{k-1}*B from (A,B,C) (SLICOT TF01RD)
//!
//! M(1)=C*B, M(2)=C*A*B, ..., M(N)=C*A^{N-1}*B.

use nalgebra::DMatrix;

/// H(i, (k-1)*NB+j) = (i,j) element of M(k).
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tf01rd(
    na: usize,
    nb: usize,
    nc: usize,
    n: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    h: &mut DMatrix<f64>,
) -> i32 {
    if a.nrows() != na || a.ncols() != na {
        return -5;
    }
    if b.nrows() != na || b.ncols() != nb || c.nrows() != nc || c.ncols() != na {
        return -6;
    }
    if h.nrows() < nc || h.ncols() < n * nb {
        return -9;
    }
    if n == 0 {
        return 0;
    }

    let mut a_pow = DMatrix::identity(na, na);
    for k in 0..n {
        let m_k = c * &a_pow * b;
        for i in 0..nc {
            for j in 0..nb {
                h[(i, k * nb + j)] = m_k[(i, j)];
            }
        }
        a_pow = a * &a_pow;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tf01rd_smoke() {
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 0.0, 0.5]);
        let b = DMatrix::from_row_slice(2, 1, &[0.0, 1.0]);
        let c = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let mut h = DMatrix::zeros(1, 3);
        assert_eq!(tf01rd(2, 1, 1, 3, &a, &b, &c, &mut h), 0);
        assert!((h[(0, 0)] - 0.0).abs() < 1e-10);
        assert!((h[(0, 1)] - 1.0).abs() < 1e-10);
        assert!((h[(0, 2)] - 0.5).abs() < 1e-10);
    }
}
