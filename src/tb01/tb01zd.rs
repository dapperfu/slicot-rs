//! TB01ZD — Single-input controllable realization (SLICOT TB01ZD)
//!
//! Reduces (A, B, C) with B a vector to orthogonal canonical form: A upper Hessenberg, B(1) nonzero only.

use nalgebra::{DMatrix, DVector};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobZ {
    No,
    Factored,
    Init,
}

/// Single-input controllable form. B is N×1 (vector). On exit A is upper Hessenberg, B has only B(0) nonzero, C is C*Z.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb01zd(
    jobz: JobZ,
    a: &mut DMatrix<f64>,
    b: &mut DVector<f64>,
    c: &mut DMatrix<f64>,
    ncont: &mut usize,
    z: &mut Option<&mut DMatrix<f64>>,
    tau: &mut [f64],
    tol: f64,
) -> i32 {
    let n = a.nrows();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -5;
    }
    if (jobz == JobZ::Init || jobz == JobZ::Factored) && (z.is_none() || z.as_ref().map(|zz| zz.nrows() != n || zz.ncols() != n).unwrap_or(true)) {
        return -11;
    }
    if tau.len() < n {
        return -12;
    }
    *ncont = 0;
    if n == 0 {
        return 0;
    }
    let bnorm: f64 = b.iter().map(|x| x.abs()).sum();
    if bnorm == 0.0 {
        if let Some(ref mut zz) = *z {
            if jobz == JobZ::Init {
                zz.fill(0.0);
                for i in 0..n {
                    zz[(i, i)] = 1.0;
                }
            }
        }
        return 0;
    }
    let eps = f64::EPSILON;
    let mut toldef = tol;
    if toldef <= 0.0 {
        let anorm = a.norm();
        toldef = (n as f64) * eps * anorm.max(bnorm);
    }
    let mut b_mat = DMatrix::from_fn(n, 1, |i, _| b[i]);
    let info = crate::tb01::tb01md::tb01md(
        crate::tb01::tb01md::JobU::Init,
        crate::tb01::tb01md::Uplo::Upper,
        a,
        &mut b_mat,
        z,
    );
    if info != 0 {
        return info;
    }
    for i in 0..n {
        b[i] = b_mat[(i, 0)];
    }
    if let Some(ref zz) = *z {
        let c_new = c.clone() * (**zz).clone();
        for i in 0..p {
            for j in 0..n {
                c[(i, j)] = c_new[(i, j)];
            }
        }
    }
    for i in 1..n {
        b[i] = 0.0;
    }
    let mut j = 0_usize;
    while j < n - 1 {
        if a[(j + 1, j)].abs() <= toldef {
            break;
        }
        j += 1;
    }
    *ncont = j + 1;
    if j + 1 < n {
        a[(j + 1, j)] = 0.0;
    }
    for i in 0..n {
        tau[i] = 0.0;
    }
    0
}
