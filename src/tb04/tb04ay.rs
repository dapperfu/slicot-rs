//! TB04AY — Transfer matrix as column polynomial vectors (SLICOT TB04AY)
//!
//! Variant of TB04AD with T(s) as columns over common denominators (dual).

use nalgebra::DMatrix;

use super::tb04ad::{tb04ad, RowCol};

/// Wrapper that calls TB04AD with ROWCOL = 'C' (columns over common denominators).
/// See TB04AD for full documentation; arguments match SLICOT TB04AY.
pub fn tb04ay(
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &DMatrix<f64>,
    nr: &mut usize,
    index: &mut [i32],
    dcoeff: &mut [f64],
    lddcoe: usize,
    ucoeff: &mut [f64],
    lduco1: usize,
    lduco2: usize,
    tol1: f64,
    tol2: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
) -> i32 {
    tb04ad(
        RowCol::C,
        a,
        b,
        c,
        d,
        nr,
        index,
        dcoeff,
        lddcoe,
        ucoeff,
        lduco1,
        lduco2,
        tol1,
        tol2,
        iwork,
        dwork,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb04ay_smoke() {
        let mut a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let mut c = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut nr = 0;
        let mut index = [0i32; 1];
        let mut dcoeff = vec![0.0; 1 * 3];
        let mut ucoeff = vec![0.0; 1 * 1 * 3];
        let mut iwork = vec![0i32; 2 + 1];
        let mut dwork = vec![0.0; 50];
        let info = tb04ay(
            &mut a,
            &mut b,
            &mut c,
            &d,
            &mut nr,
            &mut index,
            &mut dcoeff,
            1,
            &mut ucoeff,
            1,
            1,
            0.0,
            0.0,
            &mut iwork,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert!(nr >= 1);
    }
}
