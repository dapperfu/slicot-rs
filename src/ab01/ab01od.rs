//! AB01OD — Staircase form for (A,B) with optional forward/backward stages (SLICOT AB01OD)
//!
//! Reduces (A,B) to upper staircase form: forward stage via AB01ND (orthogonal
//! canonical form), backward stage via RQ factorization to triangularize blocks.

use nalgebra::DMatrix;

use crate::ab01::ab01nd::{ab01nd, JobZ};

/// Reduction stages: Forward only, Backward only, or All.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Stages {
    Forward,
    Backward,
    All,
}

/// Whether to accumulate transformation matrix.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobUV {
    No,
    Init,
}

/// Reduces (A,B) to staircase form. Forward stage uses AB01ND; backward stage
/// triangularizes blocks via RQ.
///
/// # Returns
/// 0 on success; < 0 = invalid argument index.
pub fn ab01od(
    stages: Stages,
    jobu: JobUV,
    jobv: JobUV,
    n: usize,
    m: usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    mut u: Option<&mut DMatrix<f64>>,
    mut v: Option<&mut DMatrix<f64>>,
    ncont: &mut usize,
    indcon: &mut usize,
    kstair: &mut [i32],
    tol: f64,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -6;
    }
    if b.nrows() != n || b.ncols() != m {
        return -8;
    }
    let lstagb = matches!(stages, Stages::Backward);
    let lstgab = matches!(stages, Stages::All);

    if !lstagb {
        let jobz = match jobu {
            JobUV::No => JobZ::No,
            JobUV::Init => JobZ::Init,
        };
        let info = ab01nd(jobz, n, m, a, b, ncont, indcon, kstair, u.as_deref_mut(), tol);
        if info != 0 {
            return info;
        }
    }

    if !lstgab {
        return 0;
    }
    if *ncont == 0 || *indcon == 0 {
        if matches!(jobv, JobUV::Init) {
            if let Some(ref mut vmat) = v {
                if vmat.nrows() >= m && vmat.ncols() >= m {
                    vmat.fill(0.0);
                    for i in 0..m {
                        vmat[(i, i)] = 1.0;
                    }
                }
            }
        }
        return 0;
    }

    let mut mrt = kstair[*indcon - 1] as usize;
    let mut i0 = *ncont - mrt;

    for ibstep in (2..=*indcon).rev() {
        let nrt = kstair[ibstep - 2] as usize;
        let j0 = i0 - nrt;
        let mm = nrt.min(mrt);

        let a_block = a.view((i0, j0), (mrt, nrt)).clone_owned();
        let at = a_block.transpose();
        let qr_at = at.qr();
        let r = qr_at.r().transpose();
        let q = qr_at.q().transpose();
        a.view_mut((i0, j0), (mrt, nrt)).copy_from(&r);

        let jini = if ibstep > 2 {
            j0 - kstair[ibstep - 3] as usize
        } else {
            0
        };

        if ibstep == 2 {
            let b_block = b.view_mut((0, 0), (nrt, m)).clone_owned();
            b.view_mut((0, 0), (nrt, m)).copy_from(&(q.transpose() * b_block));
        }

        let a_row = a.view((j0, jini), (nrt, n - jini)).clone_owned();
        a.view_mut((j0, jini), (nrt, n - jini)).copy_from(&(q.transpose() * a_row));

        let a_col = a.view((0, j0), (i0, nrt)).clone_owned();
        a.view_mut((0, j0), (i0, nrt)).copy_from(&(a_col * &q));

        if matches!(jobu, JobUV::Init) {
            if let Some(ref mut umat) = u {
                let u_col = umat.view((0, j0), (n, nrt)).clone_owned();
                umat.view_mut((0, j0), (n, nrt)).copy_from(&(u_col * &q));
            }
        }

        for j in 0..nrt.saturating_sub(mrt) {
            for i in 0..mrt {
                a[(i0 + i, j0 + j)] = 0.0;
            }
        }
        if i0 + mrt < n && mrt > 1 {
            for j in 1..mrt {
                for i in 0..j {
                    a[(i0 + j, i0 - mrt + i)] = 0.0;
                }
            }
        }

        mrt = nrt;
        i0 = j0;
    }

    let b1 = b.view((0, 0), (mrt, m)).clone_owned();
    let bt = b1.transpose();
    let qr_bt = bt.qr();
    let rb = qr_bt.r().transpose();
    let qb = qr_bt.q().transpose();
    b.view_mut((0, 0), (mrt, m)).copy_from(&rb);

    if matches!(jobv, JobUV::Init) {
        if let Some(ref mut vmat) = v {
            if vmat.nrows() >= m && vmat.ncols() >= m {
                vmat.fill(0.0);
                for i in 0..m {
                    vmat[(i, i)] = 1.0;
                }
                let v_block = vmat.view_mut((0, 0), (m, m)).clone_owned();
                vmat.view_mut((0, 0), (m, m)).copy_from(&(v_block * qb.transpose()));
            }
        }
    }

    for j in 0..m.saturating_sub(mrt) {
        for i in 0..mrt {
            b[(i, j)] = 0.0;
        }
    }
    if mrt > 1 {
        for j in (m - mrt)..m {
            for i in 1..mrt {
                b[(i, j)] = 0.0;
            }
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab01od_n0_m0() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut ncont = 1;
        let mut indcon = 0;
        let mut kstair = [0i32; 1];
        assert_eq!(
            ab01od(
                Stages::Forward,
                JobUV::No,
                JobUV::No,
                0,
                0,
                &mut a,
                &mut b,
                None,
                None,
                &mut ncont,
                &mut indcon,
                &mut kstair,
                0.0,
            ),
            0
        );
        assert_eq!(ncont, 0);
    }

    #[test]
    fn test_ab01od_forward_small() {
        let n = 2;
        let m = 1;
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let mut ncont = 0;
        let mut indcon = 0;
        let mut kstair = [0i32; 2];
        assert_eq!(
            ab01od(
                Stages::Forward,
                JobUV::No,
                JobUV::No,
                n,
                m,
                &mut a,
                &mut b,
                None,
                None,
                &mut ncont,
                &mut indcon,
                &mut kstair,
                0.0,
            ),
            0
        );
        assert!(ncont >= 1);
    }
}
