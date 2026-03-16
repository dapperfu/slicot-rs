//! SG03AY — Solution of generalized continuous-time Lyapunov equation (SLICOT SG03AY)
//!
//! Solves the generalized continuous-time Lyapunov equation for the Cholesky factor U
//! with A, E in generalized Schur form. Uses SG03BT (which uses SG03BX and SG03BW).

use nalgebra::DMatrix;
use crate::sg03::sg03bt::sg03bt;
use crate::sg03::sg03bx::Trans;

/// Solves for the Cholesky factor U of X = U'*U satisfying
/// A'*X*E + E'*X*A = -scale^2*B'*B (or transposed form).
///
/// A, E must be in generalized Schur form (A quasitriangular, E upper triangular).
/// B is upper triangular with non-negative diagonal; on exit overwritten with U.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1 Sylvester singular; 2 eigenvalues not complex conjugate; 3 not c-stable; 4 DSYEVX failed.
pub fn sg03ay(
    trans: Trans,
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &mut DMatrix<f64>,
    scale: &mut f64,
    dwork: &mut [f64],
    info: &mut i32,
) {
    sg03bt(trans, a, e, b, scale, dwork, info);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sg03ay_n0() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let e = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        let mut scale = 0.0;
        let mut dwork = [0.0; 0];
        let mut info = -1;
        sg03ay(Trans::NoTrans, &a, &e, &mut b, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0);
        assert_eq!(scale, 1.0);
    }

    #[test]
    fn sg03ay_n1_cstable() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let e = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut scale = 0.0;
        let mut dwork = [0.0; 1];
        let mut info = -1;
        sg03ay(Trans::NoTrans, &a, &e, &mut b, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0);
        assert!(b[(0, 0)] >= 0.0);
    }
}
