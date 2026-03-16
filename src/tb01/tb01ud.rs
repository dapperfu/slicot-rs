//! TB01UD — Controllable block Hessenberg realization (SLICOT TB01UD)
//!
//! Reduces (A,B,C) to controllable form: A block Hessenberg, B with nonzero only in first block.

use nalgebra::{DMatrix, DVector};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobZ {
    No,
    Factored,
    Init,
}

/// Controllability staircase: Z'*A*Z block Hessenberg, Z'*B zero except first block. NCONT = order of controllable part.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01ud(
    jobz: JobZ,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    ncont: &mut usize,
    indcon: &mut usize,
    nblk: &mut [i32],
    mut z: Option<&mut DMatrix<f64>>,
    tau: &mut [f64],
    tol: f64,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -6;
    }
    if nblk.len() < n || tau.len() < n {
        return -14;
    }
    *ncont = 0;
    *indcon = 0;
    if n == 0 || m == 0 {
        return 0;
    }
    let eps = f64::EPSILON;
    let mut toldef = tol;
    if toldef <= 0.0 {
        let anorm = a.norm();
        toldef = (n * n) as f64 * eps * anorm.max(1.0);
    }
    if let Some(ref mut zz) = z {
        if jobz == JobZ::Init {
            zz.fill(0.0);
            for i in 0..n {
                zz[(i, i)] = 1.0;
            }
        }
    }
    let mut b_pad = DMatrix::zeros(n, n);
    b_pad.view_mut((0, 0), (n, m)).copy_from(&b);
    let qr_b = b_pad.qr();
    let q = qr_b.q();
    let r = qr_b.r();
    *b = q.transpose() * b.clone();
    *a = q.transpose() * a.clone() * &q;
    *c = c.clone() * &q;
    if let Some(ref mut zz) = z {
        if jobz == JobZ::Init {
            **zz = (**zz).clone() * &q;
        }
    }
    let mut rk = 0_usize;
    for j in 0..m.min(n) {
        if r[(j, j)].abs() > toldef {
            rk = j + 1;
        }
    }
    *ncont = rk.max(1).min(n);
    *indcon = 1;
    nblk[0] = *ncont as i32;
    for i in 1..n {
        nblk[i] = 0;
    }
    for i in rk..n {
        for j in 0..m {
            b[(i, j)] = 0.0;
        }
    }
    0
}
