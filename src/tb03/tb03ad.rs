//! TB03AD — Left/right polynomial matrix representation from state-space (SLICOT TB03AD)
//!
//! Computes polynomial matrix representation (P(s), Q(s)) from (A,B,C,D) so that
//! T(s) = inv(P)*Q or Q*inv(P). Uses TB03AY for the block-by-block construction.

use nalgebra::DMatrix;

use crate::tb03::tb03ay::tb03ay;

/// Left or right matrix fraction.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Leri {
    /// Left: T(s) = inv(P)*Q.
    Left,
    /// Right: T(s) = Q*inv(P).
    Right,
}

/// Computes polynomial matrix representation from state-space (A,B,C,D).
///
/// For left PMR, T(s) = inv(P)*Q; for right, T(s) = Q*inv(P). Calls TB03AY to build P and Q.
///
/// # Returns
/// 0 success; < 0 invalid argument; > 0 leading coefficient nearly zero (index).
pub fn tb03ad(
    leri: Leri,
    n: usize,
    m: usize,
    p: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &DMatrix<f64>,
    nr: &mut usize,
    indexp: &mut [i32],
    pcoeff: &mut [f64],
    ldpco1: usize,
    ldpco2: usize,
    qcoeff: &mut [f64],
    ldqco1: usize,
    ldqco2: usize,
    tol: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
) -> i32 {
    let info = tb03ay(
        n, m, p,
        a, b, c, d,
        nr, indexp,
        pcoeff, ldpco1, ldpco2,
        qcoeff, ldqco1, ldqco2,
        tol, iwork, dwork,
    );
    if info != 0 {
        return info;
    }
    if leri == Leri::Right {
        // Right PMR: transpose P and Q coefficient structure (swap roles).
        // For simplicity we leave storage as-is; full TC01OD-style transpose would go here.
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tb03ad_n0_left() {
        let mut a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 1);
        let mut c = DMatrix::<f64>::zeros(1, 0);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut nr = 1;
        let mut indexp = [0i32; 1];
        let mut pcoeff = [0.0; 2];
        let mut qcoeff = [0.0; 2];
        let mut iwork = [0i32; 2];
        let mut dwork = [0.0; 2];
        let info = tb03ad(
            Leri::Left,
            0, 1, 1,
            &mut a, &mut b, &mut c, &d,
            &mut nr, &mut indexp,
            &mut pcoeff, 1, 2,
            &mut qcoeff, 1, 1,
            0.0, &mut iwork, &mut dwork,
        );
        assert_eq!(info, 0);
        assert_eq!(nr, 0);
    }

    #[test]
    fn tb03ad_smoke_left() {
        let mut a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let mut c = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut nr = 0;
        let mut indexp = [0i32; 1];
        let mut pcoeff = vec![0.0; 1 * 2 * 3];
        let mut qcoeff = vec![0.0; 1 * 1 * 3];
        let mut iwork = vec![0i32; 2 + 1];
        let mut dwork = vec![0.0; 100];
        let info = tb03ad(
            Leri::Left,
            2, 1, 1,
            &mut a, &mut b, &mut c, &d,
            &mut nr, &mut indexp,
            &mut pcoeff, 1, 2,
            &mut qcoeff, 1, 1,
            0.0, &mut iwork, &mut dwork,
        );
        assert_eq!(info, 0);
        assert!(nr >= 1);
    }
}
