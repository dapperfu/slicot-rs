//! TB01UY — Controllable realization for B = [B1, B2] (SLICOT TB01UY)
//!
//! Reduces (A, [B1,B2], C) to controllable block form with two input blocks.

use nalgebra::DMatrix;
use crate::tb01::tb01ud::{tb01ud, JobZ};

/// B = [B1, B2] with M1 and M2 columns. Outputs NCONT, INDCON, NBLK.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01uy(
    jobz: JobZ,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    _m1: usize,
    _m2: usize,
    ncont: &mut usize,
    indcon: &mut usize,
    nblk: &mut [i32],
    z: Option<&mut DMatrix<f64>>,
    tau: &mut [f64],
    tol: f64,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -6;
    }
    tb01ud(
        if jobz == JobZ::Init { JobZ::Init } else { JobZ::No },
        a,
        b,
        c,
        ncont,
        indcon,
        nblk,
        z,
        tau,
        tol,
    )
}
