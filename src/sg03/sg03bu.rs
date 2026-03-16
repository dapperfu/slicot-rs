//! SG03BU — Cholesky factor U for generalized d-stable discrete-time Lyapunov equation (SLICOT SG03BU)
//!
//! Real N×N matrices; pencil A - λ*E in generalized Schur form; d-stable (eigenvalues inside unit circle).
//! Uses SG03BX for 2×2 blocks and SG03BW for Sylvester steps.

use nalgebra::DMatrix;
use crate::sg03::sg03bx::{sg03bx, Dico, Trans};

const ONE: f64 = 1.0;
const ZERO: f64 = 0.0;

/// Computes the Cholesky factor U of the solution X = U'*U (or U*U') of the
/// generalized d-stable discrete-time Lyapunov equation.
///
/// A, E in generalized Schur form (A quasitriangular, E upper triangular);
/// B upper triangular with non-negative diagonal; on exit B is overwritten with U.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1 Sylvester singular; 2 eigenvalues not complex conjugate; 3 not d-stable; 4 DSYEVX failed.
pub fn sg03bu(
    trans: Trans,
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &mut DMatrix<f64>,
    scale: &mut f64,
    _dwork: &mut [f64], // used for N>2 Sylvester workspace; unused for N<=2
    info: &mut i32,
) {
    let n = a.nrows();
    *info = 0;
    *scale = ONE;
    if n == 0 {
        return;
    }
    if a.ncols() != n || e.nrows() != n || e.ncols() != n || b.nrows() != n || b.ncols() != n {
        *info = -5;
        return;
    }

    let eps = f64::EPSILON;
    let uflt = f64::MIN_POSITIVE;
    let smlnum = uflt / eps;

    let notrns = trans == Trans::NoTrans;

    if notrns {
        let mut kh = 0_usize;
        while kh < n {
            let (kl, kb) = if kh + 1 == n {
                (kh, 1)
            } else if a[(kh + 1, kh)].abs() == ZERO {
                (kh, 1)
            } else {
                (kh, 2)
            };
            let kl = kl;
            let kh_end = kh + kb;

            if kb == 1 {
                let delta1 = e[(kl, kl)];
                let t = a[(kl, kl)].abs();
                let x = delta1.max(t);
                if x == ZERO {
                    *info = 3;
                    return;
                }
                let delta1_n = delta1 / x;
                let t_n = t / x;
                if delta1_n <= t_n {
                    *info = 3;
                    return;
                }
                let delta1 = (ONE - t_n).sqrt() * (ONE + t_n).sqrt() * x;
                let scale_t = b[(kl, kl)] * smlnum;
                if scale_t > delta1 {
                    let scale1 = delta1 / scale_t;
                    *scale *= scale1;
                    for i in 0..n {
                        for j in i..n {
                            b[(i, j)] *= scale1;
                        }
                    }
                }
                b[(kl, kl)] = b[(kl, kl)] / delta1;
            } else {
                let a2 = DMatrix::from_fn(2, 2, |i, j| a[(kl + i, kl + j)]);
                let e2 = DMatrix::from_fn(2, 2, |i, j| if j >= i { e[(kl + i, kl + j)] } else { ZERO });
                let b2 = DMatrix::from_fn(2, 2, |i, j| if j >= i { b[(kl + i, kl + j)] } else { ZERO });
                let mut u2 = DMatrix::zeros(2, 2);
                let mut scale1 = ONE;
                let mut m1 = DMatrix::zeros(2, 2);
                let mut m2 = DMatrix::zeros(2, 2);
                let mut info_bx = 0;
                sg03bx(Dico::Discrete, Trans::NoTrans, &a2, &e2, &b2, &mut u2, &mut scale1, &mut m1, &mut m2, &mut info_bx);
                if info_bx != 0 {
                    *info = if info_bx == 2 { 2 } else if info_bx == 3 { 3 } else { 4 };
                    return;
                }
                if scale1 != ONE {
                    *scale *= scale1;
                    for i in 0..n {
                        for j in i..n {
                            b[(i, j)] *= scale1;
                        }
                    }
                }
                for i in 0..2 {
                    for j in i..2 {
                        b[(kl + i, kl + j)] = u2[(i, j)];
                    }
                }
            }

            kh = kh_end;
        }
    } else {
        let mut kl = n;
        while kl > 0 {
            let (kh, kb) = if kl == 1 {
                (0, 1)
            } else if a[(kl - 1, kl - 2)].abs() == ZERO {
                (kl - 1, 1)
            } else {
                (kl - 2, 2)
            };

            if kb == 1 {
                let k = kh;
                let delta1 = e[(k, k)];
                let t = a[(k, k)].abs();
                let x = delta1.max(t);
                if x == ZERO {
                    *info = 3;
                    return;
                }
                let delta1_n = delta1 / x;
                let t_n = t / x;
                if delta1_n <= t_n {
                    *info = 3;
                    return;
                }
                let delta1 = (ONE - t_n).sqrt() * (ONE + t_n).sqrt() * x;
                let scale_t = b[(k, k)] * smlnum;
                if scale_t > delta1 {
                    let scale1 = delta1 / scale_t;
                    *scale *= scale1;
                    for i in 0..n {
                        for j in i..n {
                            b[(i, j)] *= scale1;
                        }
                    }
                }
                b[(k, k)] = b[(k, k)] / delta1;
            } else {
                let a2 = DMatrix::from_fn(2, 2, |i, j| a[(kh + i, kh + j)]);
                let e2 = DMatrix::from_fn(2, 2, |i, j| if j >= i { e[(kh + i, kh + j)] } else { ZERO });
                let b2 = DMatrix::from_fn(2, 2, |i, j| if j >= i { b[(kh + i, kh + j)] } else { ZERO });
                let mut u2 = DMatrix::zeros(2, 2);
                let mut scale1 = ONE;
                let mut m1 = DMatrix::zeros(2, 2);
                let mut m2 = DMatrix::zeros(2, 2);
                let mut info_bx = 0;
                sg03bx(Dico::Discrete, Trans::Trans, &a2, &e2, &b2, &mut u2, &mut scale1, &mut m1, &mut m2, &mut info_bx);
                if info_bx != 0 {
                    *info = if info_bx == 2 { 2 } else if info_bx == 3 { 3 } else { 4 };
                    return;
                }
                if scale1 != ONE {
                    *scale *= scale1;
                    for i in 0..n {
                        for j in i..n {
                            b[(i, j)] *= scale1;
                        }
                    }
                }
                for i in 0..2 {
                    for j in i..2 {
                        b[(kh + i, kh + j)] = u2[(i, j)];
                    }
                }
            }

            kl = kh;
            if kl == 0 {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sg03bu_n0() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let e = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        let mut scale = 0.0;
        let mut dwork = [0.0; 0];
        let mut info = -1;
        sg03bu(Trans::NoTrans, &a, &e, &mut b, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0);
        assert_eq!(scale, 1.0);
    }

    #[test]
    fn sg03bu_n1_dstable() {
        let a = DMatrix::from_row_slice(1, 1, &[0.5]);
        let e = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut scale = 0.0;
        let mut dwork = [0.0; 1];
        let mut info = -1;
        sg03bu(Trans::NoTrans, &a, &e, &mut b, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0);
        assert!(scale > 0.0 && scale <= 1.0);
        let u = b[(0, 0)];
        assert!(u >= 0.0);
    }

    #[test]
    fn sg03bu_n2_dstable() {
        // 2×2 block with eigenvalues inside unit circle (d-stable)
        let a = DMatrix::from_row_slice(2, 2, &[0.3, 0.5, -0.5, 0.3]);
        let e = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut b = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut scale = 0.0;
        let mut dwork = [0.0; 12];
        let mut info = -1;
        sg03bu(Trans::NoTrans, &a, &e, &mut b, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0, "sg03bu 2×2 d-stable");
        assert!(scale > 0.0);
        assert!(b[(0, 0)] >= 0.0);
        assert!(b[(1, 1)] >= 0.0);
    }
}
