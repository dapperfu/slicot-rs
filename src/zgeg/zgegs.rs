//! ZGEGS — Generalized complex Schur form (SLICOT/LAPACK auxiliary).
//!
//! Computes A = Q*S*Z^H, B = Q*T*Z^H with S, T upper triangular. When B is nonsingular,
//! uses Schur of B^{-1}*A then QR of B*Z to get Q and triangular T, S.

use nalgebra::linalg::LU;
use nalgebra::DMatrix;
use num_complex::Complex64;

/// Computes the generalized Schur form of (A, B). On exit, A is overwritten with S and B with T.
/// Optionally returns VSL = Q and VSR = Z (unitary) such that A = Q*S*Z^H, B = Q*T*Z^H.
///
/// # Returns
/// 0 = success; 1 = B singular or Schur/QR failed; < 0 = invalid argument.
pub fn zgegs(
    n: usize,
    a: &mut DMatrix<Complex64>,
    b: &mut DMatrix<Complex64>,
    mut vsl: Option<&mut DMatrix<Complex64>>,
    mut vsr: Option<&mut DMatrix<Complex64>>,
) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n {
        return -2;
    }
    if b.nrows() < n || b.ncols() < n {
        return -3;
    }
    let a_sub = a.view((0, 0), (n, n)).into_owned();
    let b_sub = b.view((0, 0), (n, n)).into_owned();
    let lu = LU::new(b_sub.clone());
    let b_inv = match lu.try_inverse() {
        Some(inv) => inv,
        None => return 1,
    };
    let c = &b_inv * &a_sub;
    let schur = match c.try_schur(1e-14, 200) {
        Some(s) => s,
        None => return 1,
    };
    let (q_schur, t_schur) = schur.unpack();
    // Z = Q_schur (right generalized eigenvectors of (A,B)); eigenvalues on diagonal of T_schur
    let z = q_schur;
    let bz = &b_sub * &z;
    let qr = bz.qr();
    let q_new = qr.q();
    let r = qr.r();
    // S = R * T_schur (upper), T = R (upper)
    let s = &r * &t_schur;
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = s[(i, j)];
            b[(i, j)] = r[(i, j)];
        }
    }
    if let Some(ref mut v) = vsl {
        if v.nrows() >= n && v.ncols() >= n {
            for i in 0..n {
                for j in 0..n {
                    v[(i, j)] = q_new[(i, j)];
                }
            }
        }
    }
    if let Some(ref mut v) = vsr {
        if v.nrows() >= n && v.ncols() >= n {
            for i in 0..n {
                for j in 0..n {
                    v[(i, j)] = z[(i, j)];
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zgegs_n_zero() {
        let mut a = DMatrix::from_element(1, 1, Complex64::new(1.0, 0.0));
        let mut b = DMatrix::from_element(1, 1, Complex64::new(1.0, 0.0));
        assert_eq!(zgegs(0, &mut a, &mut b, None, None), 0);
    }

    #[test]
    fn test_zgegs_1x1() {
        let mut a = DMatrix::from_element(1, 1, Complex64::new(2.0, 0.0));
        let mut b = DMatrix::from_element(1, 1, Complex64::new(1.0, 0.0));
        assert_eq!(zgegs(1, &mut a, &mut b, None, None), 0);
        assert!((a[(0, 0)].re - 2.0).abs() < 1e-10);
        assert!(b[(0, 0)].re.abs() - 1.0 < 1e-10);
    }
}
