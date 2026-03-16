//! SG03AD — Solution of generalized discrete-time Lyapunov equation (driver) (SLICOT SG03AD)
//!
//! Driver for the generalized discrete-time Lyapunov equation. With (A,E) in generalized
//! Schur form, solves for the Cholesky factor U. Uses SG03AX.

use nalgebra::DMatrix;
use crate::sg03::sg03ax::sg03ax;
use crate::sg03::sg03bx::Trans;

/// Driver for the Cholesky factor U of X = U'*U satisfying
/// A'*X*A - E'*X*E = -scale^2*B'*B.
///
/// A, E must be in generalized Schur form. B is upper triangular; on exit overwritten with U.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1 Sylvester singular; 2 eigenvalues not complex conjugate; 3 not d-stable; 4 DSYEVX failed.
pub fn sg03ad(
    trans: Trans,
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &mut DMatrix<f64>,
    scale: &mut f64,
    dwork: &mut [f64],
    info: &mut i32,
) {
    sg03ax(trans, a, e, b, scale, dwork, info);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sg03ad_n1() {
        let a = DMatrix::from_row_slice(1, 1, &[0.5]);
        let e = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut scale = 0.0;
        let mut dwork = [0.0; 1];
        let mut info = -1;
        sg03ad(Trans::NoTrans, &a, &e, &mut b, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0);
        assert!(b[(0, 0)] >= 0.0);
    }
}
