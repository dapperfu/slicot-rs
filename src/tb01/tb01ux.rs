//! TB01UX — Observable-unobservable decomposition (SLICOT TB01UX)
//!
//! Orthogonal Z such that Z'*A*Z = [Ano *; 0 Ao], C*Z = [0 Co]; (Ao,Bo,Co) observable.

use nalgebra::DMatrix;
use crate::tb01::tb01ud::{tb01ud, JobZ};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum CompZ {
    No,
    Init,
}

/// Decomposes into observable (Ao,Bo,Co) and unobservable part. NOBSV = order of Ao.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01ux(
    compz: CompZ,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    z: Option<&mut DMatrix<f64>>,
    nobsv: &mut usize,
    nlblck: &mut usize,
    ctau: &mut [i32],
    tol: f64,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -5;
    }
    *nobsv = 0;
    *nlblck = 0;
    if n == 0 {
        return 0;
    }
    let mut at = a.transpose();
    let mut ct = c.transpose();
    let mut bt = b.transpose();
    let mut ncont = 0_usize;
    let mut indcon = 0_usize;
    let mut nblk = vec![0_i32; n];
    let mut tau = vec![0.0_f64; n];
    let mut zt = z.map(|_| DMatrix::identity(n, n));
    let info = tb01ud(if compz == CompZ::Init { JobZ::Init } else { JobZ::No }, &mut at, &mut ct, &mut bt, &mut ncont, &mut indcon, &mut nblk, zt.as_mut(), &mut tau, tol);
    if info != 0 {
        return info;
    }
    *nobsv = ncont;
    *nlblck = indcon;
    for i in 0..indcon.min(ctau.len()) {
        ctau[i] = nblk[i];
    }
    *a = at.transpose();
    *b = bt.transpose();
    *c = ct.transpose();
    if let (Some(zout), Some(zt)) = (z, zt) {
        *zout = zt.transpose();
    }
    0
}
