//! TB01PX — Reduced (minimal/controllable/observable) state-space variant (SLICOT TB01PX)
//!
//! Same as TB01PD but returns INFRED and block structure in IWORK.

use nalgebra::DMatrix;
use crate::tb01::tb01ud::{tb01ud, JobZ};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tb01PxJob {
    Minimal,
    Controllable,
    Observable,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    Scale,
    No,
}

/// Reduced realization; INFRED(1:4) and IWORK contain block info.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01px(
    job: Tb01PxJob,
    _equil: Equil,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    nr: &mut usize,
    infred: &mut [i32],
    tol: f64,
    iwork: &mut [i32],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -5;
    }
    if infred.len() < 4 {
        return -12;
    }
    infred[0] = -1;
    infred[1] = -1;
    infred[2] = 0;
    infred[3] = 0;
    *nr = n;
    if n == 0 {
        return 0;
    }
    let mut ncont = 0_usize;
    let mut indcon = 0_usize;
    let mut nblk = vec![0_i32; n];
    let mut tau = vec![0.0_f64; n];
    if job == Tb01PxJob::Minimal || job == Tb01PxJob::Controllable {
        let info = tb01ud(JobZ::No, a, b, c, &mut ncont, &mut indcon, &mut nblk, None, &mut tau, tol);
        if info != 0 {
            return info;
        }
        infred[0] = (n - ncont) as i32;
        *nr = ncont;
        for i in 0..indcon.min(iwork.len()) {
            iwork[i] = nblk[i];
        }
        infred[3] = indcon as i32;
    }
    if job == Tb01PxJob::Minimal || job == Tb01PxJob::Observable {
        let mut at = a.transpose();
        let mut ct = c.transpose();
        let mut bt = b.transpose();
        let mut nobs = 0_usize;
        let info = tb01ud(JobZ::No, &mut at, &mut ct, &mut bt, &mut nobs, &mut indcon, &mut nblk, None, &mut tau, tol);
        if info != 0 {
            return info;
        }
        infred[1] = (*nr - nobs) as i32;
        *nr = nobs;
        *a = at.transpose();
        *b = bt.transpose();
        *c = ct.transpose();
        infred[3] = indcon as i32;
    }
    0
}
