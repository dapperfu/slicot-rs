//! SG03BD — Solution of generalized discrete-time Lyapunov equation (Cholesky driver) (SLICOT SG03BD)
//!
//! Cholesky-factor driver for the generalized discrete-time Lyapunov equation.
//! With (A,E) in Schur form, uses SG03BS.

use nalgebra::DMatrix;
use crate::sg03::sg03bs::sg03bs;
use crate::sg03::sg03bx::Trans;

/// Computes the Cholesky factor U of the solution of the generalized discrete-time
/// Lyapunov equation. A, E in Schur form; B overwritten with U.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1–4 from SG03BS.
pub fn sg03bd(
    trans: Trans,
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &mut DMatrix<f64>,
    scale: &mut f64,
    dwork: &mut [f64],
    info: &mut i32,
) {
    sg03bs(trans, a, e, b, scale, dwork, info);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sg03bd_n0() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let e = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        let mut scale = 0.0;
        let mut dwork = [0.0; 0];
        let mut info = -1;
        sg03bd(Trans::NoTrans, &a, &e, &mut b, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0);
    }
}
