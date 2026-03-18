//! AB01ND — Find controllable realization for multi-input system (SLICOT AB01ND)
//!
//! Reduces (A,B) to orthogonal canonical form (block Hessenberg) using
//! QR-based reduction with rank estimation from R diagonal.

use nalgebra::DMatrix;

use crate::mb01::mb01pd::{mb01pd, Mb01PdScun, Mb01PdType};

/// Whether to accumulate the orthogonal transformation matrix Z.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobZ {
    No,
    Factored,
    Init,
}

/// Finds a controllable realization for the multi-input system (A,B).
/// Reduces (A,B) to block Hessenberg form (Ac, Bc) with Ac = Z'*A*Z, Bc = Z'*B.
///
/// # Returns
/// 0 on success; < 0 = invalid argument index (-i for i-th argument).
pub fn ab01nd(
    jobz: JobZ,
    n: usize,
    m: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    ncont: &mut usize,
    indcon: &mut usize,
    nblk: &mut [i32],
    mut z: Option<&mut DMatrix<f64>>,
    tol: f64,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -4;
    }
    if b.nrows() != n || b.ncols() != m {
        return -6;
    }
    let ljobi = matches!(jobz, JobZ::Init);
    let ljobf = matches!(jobz, JobZ::Factored);
    let ljobz = ljobi || ljobf;
    if let Some(ref zmat) = z {
        if ljobz && (zmat.nrows() < n || zmat.ncols() < n) {
            return -12;
        }
    } else if ljobz {
        return -12;
    }

    *ncont = 0;
    *indcon = 0;
    if n == 0 || m == 0 {
        if n > 0 && ljobi {
            if let Some(ref mut zmat) = z {
                zmat.fill(0.0);
                for i in 0..n {
                    zmat[(i, i)] = 1.0;
                }
            }
        }
        return 0;
    }

    let anorm = a.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let bnorm = b.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);

    if bnorm == 0.0 {
        if ljobi {
            if let Some(ref mut zmat) = z {
                zmat.fill(0.0);
                for i in 0..n {
                    zmat[(i, i)] = 1.0;
                }
            }
        }
        return 0;
    }

    let mut nblk_local = vec![0i32; n];
    let _ = mb01pd(
        Mb01PdScun::Scale,
        Mb01PdType::General,
        n,
        n,
        0,
        0,
        anorm,
        0,
        &nblk_local,
        a,
    );
    let _ = mb01pd(
        Mb01PdScun::Scale,
        Mb01PdType::General,
        n,
        m,
        0,
        0,
        bnorm,
        0,
        &nblk_local,
        b,
    );

    let _fnrm = b.norm();
    let toldef = if tol <= 0.0 {
        (n * n) as f64 * f64::EPSILON
    } else {
        tol
    };

    let mut ni = 0usize;
    let mut nrt = n;
    let mut mrt = m;
    let mut iqr = 0usize;

    if ljobi {
        if let Some(ref mut zmat) = z {
            zmat.fill(0.0);
            for i in 0..n {
                zmat[(i, i)] = 1.0;
            }
        }
    }

    loop {
        let nrt_i = nrt;
        let mrt_i = mrt;
        if nrt_i == 0 || mrt_i == 0 {
            break;
        }
        let b_block = b.view((iqr, 0), (nrt_i, mrt_i));
        let fnrm_block = b_block.norm();
        if fnrm_block == 0.0 {
            break;
        }
        let b_work = b_block.clone_owned();
        let qr = b_work.qr();
        let r = qr.r();
        let q_thin = qr.q();

        let mut rank = 0usize;
        let r_min = nrt_i.min(mrt_i);
        let thresh = toldef * fnrm_block.max(1e-300);
        for i in 0..r_min {
            if r[(i, i)].abs() > thresh {
                rank = i + 1;
            } else {
                break;
            }
        }

        if rank == 0 {
            break;
        }

        let q = if q_thin.ncols() < nrt_i {
            let mut full = DMatrix::identity(nrt_i, nrt_i);
            full.view_mut((0, 0), (nrt_i, q_thin.ncols()))
                .copy_from(&q_thin);
            full
        } else {
            q_thin
        };

        let nj = ni;
        ni = *ncont;
        *ncont += rank;
        *indcon += 1;
        nblk[*indcon - 1] = rank as i32;
        nblk_local[*indcon - 1] = rank as i32;

        let mut a_block = a.view_mut((ni, ni), (nrt_i, nrt_i));
        let a_bl = a_block.clone_owned();
        let tmp = q.transpose() * &a_bl;
        a_block.copy_from(&(tmp * &q));

        let a_col = a.view((0, ni), (n, nrt_i)).clone_owned();
        a.view_mut((0, ni), (n, nrt_i)).copy_from(&(a_col * &q));

        for j in 0..r.ncols().min(mrt_i) {
            for i in 0..r.nrows() {
                b[(iqr + i, j)] = r[(i, j)];
            }
        }
        for j in 0..rank {
            for i in j + 1..nrt_i {
                b[(iqr + i, j)] = 0.0;
            }
        }

        if ljobi {
            if let Some(ref mut zmat) = z {
                let z_cols = zmat.view_mut((0, ni), (n, nrt_i)).clone_owned();
                zmat.view_mut((0, ni), (n, nrt_i)).copy_from(&(z_cols * &q));
            }
        }

        if *indcon == 1 {
            iqr = rank;
        } else {
            for j in 0..mrt_i {
                for i in 0..rank {
                    a[(ni + i, nj + j)] = b[(iqr + i, j)];
                }
            }
        }

        if rank != nrt_i {
            mrt = rank;
            nrt = nrt_i - rank;
            let a_next = a.view((*ncont, ni), (nrt, mrt));
            for j in 0..mrt {
                for i in 0..nrt {
                    b[(iqr + i, j)] = a_next[(i, j)];
                }
            }
            a.view_mut((*ncont, ni), (nrt, mrt)).fill(0.0);
        } else {
            break;
        }
    }

    for i in iqr..n {
        for j in 0..m {
            b[(i, j)] = 0.0;
        }
    }

    let nbl = if *indcon < n {
        nblk_local[*indcon] = (n - *ncont) as i32;
        *indcon + 1
    } else {
        0
    };

    let _ = mb01pd(
        Mb01PdScun::Unscale,
        Mb01PdType::Hessenberg,
        n,
        n,
        0,
        0,
        anorm,
        nbl as i32,
        &nblk_local,
        a,
    );
    let _ = mb01pd(
        Mb01PdScun::Unscale,
        Mb01PdType::General,
        *ncont,
        m,
        0,
        0,
        bnorm,
        0,
        nblk,
        b,
    );

    0
}
