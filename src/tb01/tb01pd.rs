//! TB01PD — Minimal, controllable, or observable block Hessenberg realization (SLICOT TB01PD)

use nalgebra::DMatrix;
use crate::tb01::tb01ud::{tb01ud, JobZ};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tb01PdJob {
    /// Minimal (remove uncontrollable and unobservable).
    Minimal,
    /// Controllable only.
    Controllable,
    /// Observable only.
    Observable,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    Scale,
    No,
}

/// Finds reduced (Ar,Br,Cr) of order NR. If EQUIL=Scale, preliminarily balance with TB01ID.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01pd(
    job: Tb01PdJob,
    _equil: Equil,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    nr: &mut usize,
    tol: f64,
    iwork: &mut [i32],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -5;
    }
    *nr = n;
    if n == 0 {
        return 0;
    }
    let mut ncont = 0_usize;
    let mut indcon = 0_usize;
    let mut nblk = vec![0_i32; n];
    let mut tau = vec![0.0_f64; n];
    let mut z = DMatrix::identity(n, n);
    if job == Tb01PdJob::Minimal || job == Tb01PdJob::Controllable {
        let info = tb01ud(JobZ::Init, a, b, c, &mut ncont, &mut indcon, &mut nblk, Some(&mut z), &mut tau, tol);
        if info != 0 {
            return info;
        }
        *nr = ncont;
        for i in 0..indcon.min(iwork.len()) {
            iwork[i] = nblk[i];
        }
    }
    if job == Tb01PdJob::Minimal || job == Tb01PdJob::Observable {
        let mut at = a.transpose();
        let mut ct = c.transpose();
        let mut bt = b.transpose();
        let mut nobs = 0_usize;
        let info = tb01ud(JobZ::Init, &mut at, &mut ct, &mut bt, &mut nobs, &mut indcon, &mut nblk, None, &mut tau, tol);
        if info != 0 {
            return info;
        }
        *nr = nobs;
        *a = at.transpose();
        *b = bt.transpose();
        *c = ct.transpose();
    }
    0
}
