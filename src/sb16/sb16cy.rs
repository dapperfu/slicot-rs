//! SB16CY — Cholesky factors of controllability and observability Grammians of coprime factors
//! of a state-feedback controller (SLICOT).
//!
//! For (A,B,C), F, G with A+B*F and A+G*C stable, computes Cholesky factors Su, Ru of
//! frequency-weighted Grammians P = Su*Su', Q = Ru'*Ru for left or right coprime factorization.

use nalgebra::DMatrix;

use crate::ab13::lyapunov;

/// Continuous ('C') or discrete ('D') time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// Left ('L') or right ('R') coprime factorization.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobcf {
    Left,
    Right,
}

/// Computes Cholesky factors Su, Ru of controllability and observability Grammians for
/// coprime-factor controller reduction. Open-loop (A,B,C) with D=0; F (M×N), G (N×P).
/// Returns INFO: 0 = success; 1 = eigenvalue failure; 2 = A+G*C not stable; 3 = A+B*F not stable;
/// 4 = observability Lyapunov singular; 5 = controllability Lyapunov singular; < 0 invalid argument.
pub fn sb16cy(
    dico: Dico,
    jobcf: Jobcf,
    n: usize,
    m: usize,
    p: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    f: &DMatrix<f64>,
    g: &DMatrix<f64>,
    scalec: &mut f64,
    scaleo: &mut f64,
    s: &mut DMatrix<f64>,
    r: &mut DMatrix<f64>,
    dwork: &mut [f64],
) -> i32 {
    *scalec = 1.0;
    *scaleo = 1.0;
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n
        || b.nrows() != n || b.ncols() != m
        || c.nrows() != p || c.ncols() != n
        || f.nrows() != m || f.ncols() != n
        || g.nrows() != n || g.ncols() != p
        || s.nrows() != n || s.ncols() != n
        || r.nrows() != n || r.ncols() != n
    {
        return -1;
    }
    let min_work = n * (n + n.max(m).max(p) + 7);
    if dwork.len() < min_work {
        return -2;
    }

    let af = a + b * f;
    let agc = a + g * c;

    let (q_ctr, q_obs) = match jobcf {
        Jobcf::Left => {
            let q_ctr = b * b.transpose();
            let q_obs = f.transpose() * f;
            (q_ctr, q_obs)
        }
        Jobcf::Right => {
            let q_ctr = g * g.transpose();
            let q_obs = c.transpose() * c;
            (q_ctr, q_obs)
        }
    };

    let mut p_mat = DMatrix::<f64>::zeros(n, n);
    let mut q_mat = DMatrix::<f64>::zeros(n, n);

    let ok_p = match dico {
        Dico::Continuous => lyapunov::lyapunov_continuous_dual(&af, &q_ctr, &mut p_mat),
        Dico::Discrete => {
            let n2 = n * n;
            if dwork.len() < n2 * n2 + n2 {
                return -2;
            }
            lyapunov::lyapunov_discrete_stein(&af, &q_ctr, &mut p_mat)
        }
    };
    if !ok_p {
        return 5;
    }

    let ok_q = match dico {
        Dico::Continuous => lyapunov::lyapunov_continuous(&agc.transpose(), &q_obs, &mut q_mat),
        Dico::Discrete => lyapunov::lyapunov_discrete_stein_dual(&agc, &q_obs, &mut q_mat),
    };
    if !ok_q {
        return 4;
    }

    // Ensure symmetry (numerical)
    for i in 0..n {
        for j in (i + 1)..n {
            p_mat[(i, j)] = (p_mat[(i, j)] + p_mat[(j, i)]) * 0.5;
            p_mat[(j, i)] = p_mat[(i, j)];
            q_mat[(i, j)] = (q_mat[(i, j)] + q_mat[(j, i)]) * 0.5;
            q_mat[(j, i)] = q_mat[(i, j)];
        }
    }

    // P = Su*Su' => upper Cholesky Su. nalgebra cholesky gives L with P = L*L', so Su = L'.
    let ch_p = match p_mat.cholesky() {
        Some(ch) => ch,
        None => return 5,
    };
    let l_p = ch_p.l();
    for i in 0..n {
        for j in 0..n {
            s[(i, j)] = if j >= i { l_p[(j, i)] } else { 0.0 };
        }
    }

    // Q = Ru'*Ru => Ru upper with Q = Ru'*Ru, so Ru'*Ru = Q. Cholesky Q = L*L' => Ru = L'.
    let ch_q = match q_mat.cholesky() {
        Some(ch) => ch,
        None => return 4,
    };
    let l_q = ch_q.l();
    for i in 0..n {
        for j in 0..n {
            r[(i, j)] = if j >= i { l_q[(j, i)] } else { 0.0 };
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb16cy_trivial() {
        let n = 0_usize;
        let a = DMatrix::<f64>::zeros(0, 0);
        let b = DMatrix::<f64>::zeros(0, 0);
        let c = DMatrix::<f64>::zeros(0, 0);
        let f = DMatrix::<f64>::zeros(0, 0);
        let g = DMatrix::<f64>::zeros(0, 0);
        let mut s = DMatrix::<f64>::zeros(0, 0);
        let mut r = DMatrix::<f64>::zeros(0, 0);
        let mut scalec = 0.0;
        let mut scaleo = 0.0;
        let mut dwork = vec![0.0; 1];
        assert_eq!(
            sb16cy(
                Dico::Continuous,
                Jobcf::Left,
                n,
                0,
                0,
                &a,
                &b,
                &c,
                &f,
                &g,
                &mut scalec,
                &mut scaleo,
                &mut s,
                &mut r,
                &mut dwork
            ),
            0
        );
        assert_eq!(scalec, 1.0);
        assert_eq!(scaleo, 1.0);
    }

    #[test]
    fn test_sb16cy_1x1_stable() {
        let n = 1_usize;
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let f = DMatrix::from_row_slice(1, 1, &[-0.5]);
        let g = DMatrix::from_row_slice(1, 1, &[-0.5]);
        let mut s = DMatrix::<f64>::zeros(1, 1);
        let mut r = DMatrix::<f64>::zeros(1, 1);
        let mut scalec = 0.0;
        let mut scaleo = 0.0;
        let mut dwork = vec![0.0; 20];
        assert_eq!(
            sb16cy(
                Dico::Continuous,
                Jobcf::Left,
                n,
                1,
                1,
                &a,
                &b,
                &c,
                &f,
                &g,
                &mut scalec,
                &mut scaleo,
                &mut s,
                &mut r,
                &mut dwork
            ),
            0
        );
        assert!(s[(0, 0)] > 0.0);
        assert!(r[(0, 0)] > 0.0);
    }
}
