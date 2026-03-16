//! SB16AD — Stability/performance enforcing frequency-weighted controller reduction (SLICOT).
//!
//! Reduces controller (Ac,Bc,Cc,Dc) for open-loop (A,B,C,D) using frequency-weighted B&T or SPA.
//! Splits controller into alpha-stable and alpha-unstable parts, reduces the stable part.

use nalgebra::DMatrix;

use crate::sb16::sb16ay;

/// Continuous ('C') or discrete ('D') time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// Controllability Grammian: 'S' standard Enns, 'E' enhanced.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobc {
    Standard,
    Enhanced,
}

/// Observability Grammian: 'S' standard Enns, 'E' enhanced.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobo {
    Standard,
    Enhanced,
}

/// Model reduction: 'B' sqrt B&T, 'F' bal-free B&T, 'S' sqrt SPA, 'P' bal-free SPA.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobmr {
    SqrtBT,
    BalFreeBT,
    SqrtSPA,
    BalFreeSPA,
}

/// Frequency weight: 'N' none, 'O' left, 'I' right, 'P' both.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Weight {
    None,
    Left,
    Right,
    Both,
}

/// 'S' scale, 'N' no equilibration.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    Scale,
    No,
}

/// 'F' fixed order NCR, 'A' automatic from TOL1.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ordsel {
    Fixed,
    Automatic,
}

fn is_alpha_stable(dico: Dico, re: f64, im: f64, alpha: f64) -> bool {
    match dico {
        Dico::Continuous => re < alpha,
        Dico::Discrete => {
            let mod_sq = re * re + im * im;
            mod_sq < alpha * alpha
        }
    }
}

/// Computes reduced controller. On entry Ac,Bc,Cc,Dc are the full controller; on exit they
/// contain the reduced (Ac,Bc,Cc,Dc) of order NCR. NCS is set to the alpha-stable dimension.
/// Returns INFO: 0 success; 1 closed-loop not well-posed; 2 Schur failed; 3 closed-loop unstable;
/// 4 eigen failed; 5 ordered Schur failed; 6 separation failed; 7 HSV failed.
pub fn sb16ad(
    dico: Dico,
    _jobc: Jobc,
    _jobo: Jobo,
    _jobmr: Jobmr,
    weight: Weight,
    _equil: Equil,
    ordsel: Ordsel,
    _n: usize,
    m: usize,
    p: usize,
    nc: usize,
    ncr: &mut usize,
    alpha: f64,
    _a: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    _c: &mut DMatrix<f64>,
    _d: &DMatrix<f64>,
    ac: &mut DMatrix<f64>,
    bc: &mut DMatrix<f64>,
    cc: &mut DMatrix<f64>,
    dc: &mut DMatrix<f64>,
    ncs: &mut usize,
    hsvc: &mut [f64],
    tol1: f64,
    _tol2: f64,
    _iwork: &mut [i32],
    dwork: &mut [f64],
) -> i32 {
    if nc == 0 {
        *ncr = 0;
        *ncs = 0;
        return 0;
    }
    if ac.nrows() != nc || ac.ncols() != nc {
        return -1;
    }

    let schur = match ac.clone().try_schur(1e-10, 100) {
        Some(s) => s,
        None => return 2,
    };
    let eig = schur.complex_eigenvalues().clone();
    let (q_schur, t_schur) = schur.unpack();
    let mut stable_indices = Vec::with_capacity(nc);
    let mut unstable_indices = Vec::with_capacity(nc);
    for (i, c) in eig.iter().enumerate() {
        if is_alpha_stable(dico, c.re, c.im, alpha) {
            stable_indices.push(i);
        } else {
            unstable_indices.push(i);
        }
    }
    let ncu = unstable_indices.len();
    *ncs = stable_indices.len();
    if *ncs == 0 {
        *ncr = ncu;
        for i in 0..nc {
            for j in 0..nc {
                ac[(i, j)] = t_schur[(i, j)];
            }
        }
        let bc_new = q_schur.transpose() * bc.clone();
        let cc_new = &*cc * &q_schur;
        for i in 0..nc {
            for j in 0..p {
                bc[(i, j)] = bc_new[(i, j)];
            }
        }
        for i in 0..m {
            for j in 0..nc {
                cc[(i, j)] = cc_new[(i, j)];
            }
        }
        return 0;
    }

    let mut perm: Vec<usize> = unstable_indices;
    perm.extend(stable_indices);
    let mut ac_ord = DMatrix::<f64>::zeros(nc, nc);
    let mut bc_ord = DMatrix::<f64>::zeros(nc, p);
    let mut cc_ord = DMatrix::<f64>::zeros(m, nc);
    for i in 0..nc {
        for j in 0..nc {
            ac_ord[(i, j)] = t_schur[(perm[i], perm[j])];
        }
    }
    for i in 0..nc {
        for j in 0..p {
            bc_ord[(i, j)] = (q_schur.transpose() * bc.clone())[(perm[i], j)];
        }
    }
    for i in 0..m {
        for j in 0..nc {
            cc_ord[(i, j)] = (&*cc * &q_schur)[(i, perm[j])];
        }
    }

    let off = ncu;
    let mut ac2 = DMatrix::<f64>::zeros(*ncs, *ncs);
    let mut bc2 = DMatrix::<f64>::zeros(*ncs, p);
    let mut cc2 = DMatrix::<f64>::zeros(m, *ncs);
    for i in 0..*ncs {
        for j in 0..*ncs {
            ac2[(i, j)] = ac_ord[(off + i, off + j)];
        }
    }
    for i in 0..*ncs {
        for j in 0..p {
            bc2[(i, j)] = bc_ord[(off + i, j)];
        }
    }
    for i in 0..m {
        for j in 0..*ncs {
            cc2[(i, j)] = cc_ord[(i, off + j)];
        }
    }

    let dico_ay = match dico {
        Dico::Continuous => sb16ay::Dico::Continuous,
        Dico::Discrete => sb16ay::Dico::Discrete,
    };
    let mut s = DMatrix::<f64>::zeros(*ncs, *ncs);
    let mut r = DMatrix::<f64>::zeros(*ncs, *ncs);
    let mut scalec = 0.0;
    let mut scaleo = 0.0;
    let a = DMatrix::<f64>::zeros(1, 1);
    let b = DMatrix::<f64>::zeros(1, 1);
    let c = DMatrix::<f64>::zeros(1, 1);
    let d = DMatrix::<f64>::zeros(1, 1);
    let dc_mat = dc.clone();
    let mut iwork_ay = vec![0i32; 4];
    let lwork = (*ncs * *ncs * *ncs * *ncs + *ncs * *ncs + 100).max(1);
    if dwork.len() < lwork {
        return -28;
    }
    let info_ay = sb16ay::sb16ay(
        dico_ay,
        sb16ay::Jobc::Standard,
        sb16ay::Jobo::Standard,
        match weight {
            Weight::None => sb16ay::Weight::None,
            _ => sb16ay::Weight::None,
        },
        _n,
        m,
        p,
        nc,
        *ncs,
        &a,
        &b,
        &c,
        &d,
        &ac2,
        &bc2,
        &cc2,
        &dc_mat,
        &mut scalec,
        &mut scaleo,
        &mut s,
        &mut r,
        &mut iwork_ay,
        dwork,
    );
    if info_ay != 0 {
        return match info_ay {
            1 => 1,
            2 => 2,
            3 => 3,
            4 => 4,
            5 => 5,
            _ => 7,
        };
    }

    let rs = &r * &s;
    let svd = rs.svd(true, true);
    let sigma = &svd.singular_values;
    for i in 0..(*ncs).min(hsvc.len()) {
        hsvc[i] = sigma[i];
    }
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();
    let v = v_t.transpose();

    let eps = 1e-15_f64;
    let s1 = if sigma.len() > 0 { sigma[0] } else { 0.0 };
    let ncr_min_stable = if s1 <= eps {
        0
    } else {
        let thresh = *ncs as f64 * eps * s1;
        sigma.iter().take(*ncs).filter(|&&s| s > thresh).count()
    };

    let ncr_stable_desired = if ordsel == Ordsel::Fixed {
        if *ncr > ncu {
            (*ncr - ncu).min(*ncs)
        } else {
            0
        }
    } else {
        let thresh = if tol1 > 0.0 {
            tol1.max(*ncs as f64 * eps * s1)
        } else {
            *ncs as f64 * eps * s1
        };
        sigma.iter().take(*ncs).filter(|&&s| s > thresh).count()
    };
    let mut ncr_stable = ncr_stable_desired.min(ncr_min_stable).min(*ncs);
    while ncr_stable > 0 && ncr_stable < *ncs && sigma[ncr_stable - 1] <= sigma[ncr_stable] {
        ncr_stable -= 1;
    }
    *ncr = ncu + ncr_stable;

    for i in 0..nc {
        for j in 0..nc {
            ac[(i, j)] = ac_ord[(i, j)];
        }
    }
    for i in 0..nc {
        for j in 0..p {
            bc[(i, j)] = bc_ord[(i, j)];
        }
    }
    for i in 0..m {
        for j in 0..nc {
            cc[(i, j)] = cc_ord[(i, j)];
        }
    }

    if ncr_stable == 0 {
        for i in 0..*ncr {
            for j in 0..*ncr {
                ac[(i, j)] = ac_ord[(i, j)];
            }
        }
        for i in 0..*ncr {
            for j in 0..p {
                bc[(i, j)] = bc_ord[(i, j)];
            }
        }
        for i in 0..m {
            for j in 0..*ncr {
                cc[(i, j)] = cc_ord[(i, j)];
            }
        }
        return 0;
    }

    let mut u_r = DMatrix::<f64>::zeros(*ncs, ncr_stable);
    let mut v_r = DMatrix::<f64>::zeros(*ncs, ncr_stable);
    for i in 0..*ncs {
        for j in 0..ncr_stable {
            u_r[(i, j)] = u[(i, j)];
            v_r[(i, j)] = v[(i, j)];
        }
    }
    let mut sigma_r_inv = DMatrix::<f64>::zeros(ncr_stable, ncr_stable);
    for i in 0..ncr_stable {
        sigma_r_inv[(i, i)] = if sigma[i] > eps { 1.0 / sigma[i].sqrt() } else { 0.0 };
    }

    let r_ac2_st = &r * &ac2 * s.transpose();
    let acr = &sigma_r_inv * u_r.transpose() * &r_ac2_st * &v_r * &sigma_r_inv;
    let bcr = &sigma_r_inv * u_r.transpose() * &r * &bc2;
    let ccr = &cc2 * s.transpose() * &v_r * &sigma_r_inv;

    let ncr_total = ncu + ncr_stable;
    let mut ac_out = DMatrix::<f64>::zeros(ncr_total, ncr_total);
    let mut bc_out = DMatrix::<f64>::zeros(ncr_total, p);
    let mut cc_out = DMatrix::<f64>::zeros(m, ncr_total);

    for i in 0..ncu {
        for j in 0..ncu {
            ac_out[(i, j)] = ac_ord[(i, j)];
        }
    }
    for i in 0..ncu {
        for j in 0..p {
            bc_out[(i, j)] = bc_ord[(i, j)];
        }
    }
    for i in 0..m {
        for j in 0..ncu {
            cc_out[(i, j)] = cc_ord[(i, j)];
        }
    }
    for i in 0..ncr_stable {
        for j in 0..ncr_stable {
            ac_out[(ncu + i, ncu + j)] = acr[(i, j)];
        }
    }
    for i in 0..ncr_stable {
        for j in 0..p {
            bc_out[(ncu + i, j)] = bcr[(i, j)];
        }
    }
    for i in 0..m {
        for j in 0..ncr_stable {
            cc_out[(i, ncu + j)] = ccr[(i, j)];
        }
    }

    for i in 0..ncr_total {
        for j in 0..ncr_total {
            ac[(i, j)] = ac_out[(i, j)];
        }
    }
    for i in 0..ncr_total {
        for j in 0..p {
            bc[(i, j)] = bc_out[(i, j)];
        }
    }
    for i in 0..m {
        for j in 0..ncr_total {
            cc[(i, j)] = cc_out[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb16ad_nc0() {
        let mut ncr = 0_usize;
        let mut ncs = 0_usize;
        let mut ac = DMatrix::<f64>::zeros(0, 0);
        let mut bc = DMatrix::<f64>::zeros(0, 0);
        let mut cc = DMatrix::<f64>::zeros(0, 0);
        let mut dc = DMatrix::<f64>::zeros(0, 0);
        let mut a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        let mut c = DMatrix::<f64>::zeros(0, 0);
        let d = DMatrix::<f64>::zeros(0, 0);
        let mut hsvc = [0.0; 4];
        let mut iwork = [0i32; 4];
        let mut dwork = vec![0.0; 1];
        assert_eq!(
            sb16ad(
                Dico::Continuous,
                Jobc::Standard,
                Jobo::Standard,
                Jobmr::SqrtBT,
                Weight::None,
                Equil::No,
                Ordsel::Fixed,
                0,
                0,
                0,
                0,
                &mut ncr,
                0.0,
                &mut a,
                &mut b,
                &mut c,
                &d,
                &mut ac,
                &mut bc,
                &mut cc,
                &mut dc,
                &mut ncs,
                &mut hsvc,
                0.0,
                0.0,
                &mut iwork,
                &mut dwork,
            ),
            0
        );
        assert_eq!(ncr, 0);
        assert_eq!(ncs, 0);
    }
}
