//! SB16CD — Coprime factorization based frequency-weighted state feedback controller reduction (SLICOT).
//!
//! For (A,B,C,D), F, G with A+B*F and A+G*C stable, computes reduced (Ac,Bc,Cc) using
//! frequency-weighted B&T on coprime factors. No Dc (output is strictly proper).

use nalgebra::DMatrix;

use crate::sb16::sb16cy::sb16cy;

/// Continuous ('C') or discrete ('D') time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// 'D' D present, 'Z' D zero.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobd {
    D,
    Zero,
}

/// Model reduction: 'B' square-root B&T, 'F' balancing-free square-root B&T.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobmr {
    SqrtBT,
    BalFreeBT,
}

/// 'L' left or 'R' right coprime factorization.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobcf {
    Left,
    Right,
}

/// 'F' fixed order NCR, 'A' automatic from TOL.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ordsel {
    Fixed,
    Automatic,
}

/// Computes reduced controller (Ac,Bc,Cc) from (A,B,C,D), F, G using coprime-factor
/// frequency-weighted B&T. On exit A, B, C, F, G are overwritten with Ac, Bc, Cc (in place).
/// Returns INFO: 0 success; 1 eigen failure; 2 A+G*C not stable; 3 A+B*F not stable;
/// 4 observability Lyapunov singular; 5 controllability Lyapunov singular; 6 HSV failed.
pub fn sb16cd(
    dico: Dico,
    jobd: Jobd,
    jobmr: Jobmr,
    jobcf: Jobcf,
    ordsel: Ordsel,
    n: usize,
    m: usize,
    p: usize,
    ncr: &mut usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &DMatrix<f64>,
    f: &mut DMatrix<f64>,
    g: &mut DMatrix<f64>,
    hsv: &mut [f64],
    tol: f64,
    iwork: &mut [i32],
    dwork: &mut [f64],
) -> i32 {
    if n == 0 {
        *ncr = 0;
        return 0;
    }
    let n2 = n * n;
    let min_work = (n * (n + n.max(m).max(p) + 7)).max(2 * n2 + 5 * n).max(1);
    if dwork.len() < min_work {
        return -22;
    }
    let dico_cy = match dico {
        Dico::Continuous => crate::sb16::sb16cy::Dico::Continuous,
        Dico::Discrete => crate::sb16::sb16cy::Dico::Discrete,
    };
    let jobcf_cy = match jobcf {
        Jobcf::Left => crate::sb16::sb16cy::Jobcf::Left,
        Jobcf::Right => crate::sb16::sb16cy::Jobcf::Right,
    };

    let mut s = DMatrix::<f64>::zeros(n, n);
    let mut r = DMatrix::<f64>::zeros(n, n);
    let mut scalec = 0.0;
    let mut scaleo = 0.0;
    let a_eff = if jobd == Jobd::D {
        a.clone() + &(&*g * d * &*f)
    } else {
        a.clone()
    };
    let b_cl = g.clone();
    let c_cl = f.clone();

    let info_cy = sb16cy(
        dico_cy,
        jobcf_cy,
        n,
        m,
        p,
        &a_eff,
        b,
        c,
        &c_cl,
        &b_cl,
        &mut scalec,
        &mut scaleo,
        &mut s,
        &mut r,
        dwork,
    );
    if info_cy != 0 {
        return match info_cy {
            2 => 2,
            3 => 3,
            4 => 4,
            5 => 5,
            _ => 1,
        };
    }

    let ac = &a_eff + &*b * &*f + &*g * &*c;
    let rs = &r * &s;
    let svd = rs.svd(true, true);
    let sigma = &svd.singular_values;
    for i in 0..n.min(hsv.len()) {
        hsv[i] = sigma[i];
    }
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();
    let v = v_t.transpose();

    let eps = 1e-15_f64;
    let s1 = if sigma.len() > 0 { sigma[0] } else { 0.0 };
    let ncr_min = if s1 <= eps {
        0
    } else {
        let thresh = if ordsel == Ordsel::Automatic && tol > 0.0 {
            tol.max(n as f64 * eps * s1)
        } else {
            n as f64 * eps * s1
        };
        sigma.iter().take(n).filter(|&&s| s > thresh).count()
    };

    let ncr_desired = if ordsel == Ordsel::Fixed { *ncr } else { ncr_min };
    let mut ncr_final = ncr_desired.min(n);
    if ncr_final > ncr_min && ordsel == Ordsel::Fixed {
        ncr_final = ncr_min;
    }
    while ncr_final > 0 && ncr_final < n && sigma[ncr_final - 1] <= sigma[ncr_final] {
        ncr_final -= 1;
    }
    *ncr = ncr_final;

    if ncr_final == 0 {
        *ncr = 0;
        return 0;
    }

    let mut u_r = DMatrix::<f64>::zeros(n, ncr_final);
    let mut v_r = DMatrix::<f64>::zeros(n, ncr_final);
    for i in 0..n {
        for j in 0..ncr_final {
            u_r[(i, j)] = u[(i, j)];
            v_r[(i, j)] = v[(i, j)];
        }
    }
    let mut sigma_r_inv = DMatrix::<f64>::zeros(ncr_final, ncr_final);
    for i in 0..ncr_final {
        sigma_r_inv[(i, i)] = if sigma[i] > eps { 1.0 / sigma[i].sqrt() } else { 0.0 };
    }

    let r_ac_st = &r * &ac * s.transpose();
    let ac_red = &sigma_r_inv * u_r.transpose() * &r_ac_st * &v_r * &sigma_r_inv;
    let bc_red = &sigma_r_inv * u_r.transpose() * &r * &b_cl;
    let cc_red = &c_cl * s.transpose() * &v_r * &sigma_r_inv;

    for i in 0..ncr_final {
        for j in 0..ncr_final {
            a[(i, j)] = ac_red[(i, j)];
        }
    }
    for i in 0..ncr_final {
        for j in 0..m {
            b[(i, j)] = bc_red[(i, j)];
        }
    }
    for i in 0..p {
        for j in 0..ncr_final {
            c[(i, j)] = cc_red[(i, j)];
        }
    }
    for i in 0..m {
        for j in 0..ncr_final {
            f[(i, j)] = cc_red[(i, j)];
        }
    }
    for i in 0..ncr_final {
        for j in 0..p {
            g[(i, j)] = bc_red[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb16cd_n0() {
        let mut ncr = 0_usize;
        let mut a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        let mut c = DMatrix::<f64>::zeros(0, 0);
        let d = DMatrix::<f64>::zeros(0, 0);
        let mut f = DMatrix::<f64>::zeros(0, 0);
        let mut g = DMatrix::<f64>::zeros(0, 0);
        let mut hsv = [0.0; 4];
        let mut iwork = [0i32; 4];
        let mut dwork = vec![0.0; 1];
        assert_eq!(
            sb16cd(
                Dico::Continuous,
                Jobd::Zero,
                Jobmr::SqrtBT,
                Jobcf::Right,
                Ordsel::Fixed,
                0,
                0,
                0,
                &mut ncr,
                &mut a,
                &mut b,
                &mut c,
                &d,
                &mut f,
                &mut g,
                &mut hsv,
                0.0,
                &mut iwork,
                &mut dwork,
            ),
            0
        );
        assert_eq!(ncr, 0);
    }
}
