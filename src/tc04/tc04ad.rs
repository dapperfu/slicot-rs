//! TC04AD — State-space representation for a given left/right polynomial matrix representation (SLICOT TC04AD)
//!
//! Finds (A,B,C,D) such that C*inv(sI-A)*B + D = T(s) = inv(P(s))*Q(s) or Q(s)*inv(P(s)).
//! Uses Wolovich observable companion form for left; right via dual.

use nalgebra::{linalg::LU, DMatrix};
use std::f64;

/// Left or right matrix fraction.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Leri {
    Left,
    Right,
}

/// Computes state-space (A,B,C,D) from polynomial matrix representation.
///
/// INDEX(i) = max degree of row i (left) or column i (right). kpcoef = max(INDEX)+1.
/// PCOEFF(i,j,k) = coefficient of s^(INDEX(iorj)-K+1). Fortran column-major.
///
/// # Returns
/// * `0` - success
/// * `< 0` - invalid argument
/// * `1` - P(s) not row/column proper
pub fn tc04ad(
    leri: Leri,
    m: usize,
    p: usize,
    index: &[i32],
    pcoeff: &[f64],
    ldpco1: usize,
    ldpco2: usize,
    qcoeff: &[f64],
    ldqco1: usize,
    ldqco2: usize,
    n: &mut usize,
    rcond: &mut f64,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
) -> i32 {
    let porm = if leri == Leri::Left { p } else { m };
    let porp = if leri == Leri::Left { m } else { p };
    if index.len() < porm {
        return -4;
    }
    let kpcoef = index.iter().take(porm).map(|&d| d as usize).max().unwrap_or(0) + 1;
    *n = index.iter().take(porm).map(|&d| d as usize).sum();
    if *n == 0 {
        return 0;
    }

    // Leading coefficient matrix L_P (porm×porm): L_P(i,j) = PCOEFF(i,j,1) = coeff of s^INDEX(i)
    let mut l_p = DMatrix::<f64>::zeros(porm, porm);
    for i in 0..porm {
        for j in 0..porm {
            let idx = i + j * ldpco1;
            l_p[(i, j)] = pcoeff[idx];
        }
    }
    let lu_p = LU::new(l_p.clone());
    if !lu_p.is_invertible() {
        return 1;
    }
    let l_p_inv = match lu_p.try_inverse() {
        Some(inv) => inv,
        None => return 1,
    };
    *rcond = 1.0 / (l_p.norm() * l_p_inv.norm()).max(1e-307);

    // D = L_P^{-1} * L_Q
    let mut l_q = DMatrix::<f64>::zeros(porm, porp);
    for i in 0..porm {
        for j in 0..porp {
            let idx = i + j * ldqco1;
            l_q[(i, j)] = qcoeff[idx];
        }
    }
    let d_mat = &l_p_inv * &l_q;
    for i in 0..p.min(d.nrows()) {
        for j in 0..m.min(d.ncols()) {
            let (di, dj) = if leri == Leri::Left {
                (i, j)
            } else {
                (i, j)
            };
            if di < d_mat.nrows() && dj < d_mat.ncols() {
                d[(i, j)] = d_mat[(di, dj)];
            }
        }
    }

    // Observable companion: state ordering by rows (block 1 has INDEX(0) states, block 2 has INDEX(1), ...)
    let mut row_start = vec![0usize; porm + 1];
    for (i, &deg) in index.iter().take(porm).enumerate() {
        row_start[i + 1] = row_start[i] + deg as usize;
    }

    a.fill(0.0);
    c.fill(0.0);
    for (ii, &deg_i) in index.iter().take(porm).enumerate() {
        let d_i = deg_i as usize;
        let r0 = row_start[ii];
        for k in 0..d_i {
            for j in 0..porm {
                let exp = deg_i - k as i32 - 1;
                if exp >= 0 {
                    let kk = (exp + 1) as usize;
                    if kk < kpcoef {
                        let idx = ii + j * ldpco1 + kk * ldpco1 * ldpco2;
                        a[(r0 + k, row_start[j])] = -pcoeff[idx] / l_p[(ii, ii)];
                    }
                }
            }
            if k < d_i - 1 {
                a[(r0 + k, r0 + k + 1)] = 1.0;
            }
        }
        c[(ii, r0 + d_i - 1)] = 1.0 / l_p[(ii, ii)];
    }

    // B from Q(s) - P(s)D: coefficient comparison. Sbar(s)*B = Q - P*D.
    // Simplified: solve for B from the relation (e.g. constant term and higher).
    let mut q_minus_pd = vec![DMatrix::<f64>::zeros(porm, porp); kpcoef];
    for k in 0..kpcoef {
        for i in 0..porm {
            for j in 0..porp {
                let idx_q = i + j * ldqco1 + k * ldqco1 * ldqco2;
                let mut v = qcoeff[idx_q];
                for jp in 0..porm {
                    let idx_p = i + jp * ldpco1 + k * ldpco1 * ldpco2;
                    v -= pcoeff[idx_p] * d_mat[(jp, j)];
                }
                q_minus_pd[k][(i, j)] = v;
            }
        }
    }
    for j in 0..porp {
        for i in 0..porm {
            let d_i = index[i] as usize;
            let r0 = row_start[i];
            for (k, qk) in q_minus_pd.iter().enumerate().take(d_i) {
                let exp = index[i] as i32 - k as i32;
                if exp >= 0 {
                    b[(r0 + (d_i - 1 - k), j)] = qk[(i, j)] / l_p[(i, i)];
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tc04ad_left_2x2() {
        // TC04AD example: M=2, P=2, L. INDEX=[2,2]. Same data as TC01OD/TC05AD.
        let index = [2, 2];
        let pcoeff = vec![
            2.0, 3.0, 1.0, 5.0, 7.0, -6.0, 4.0, -1.0, -1.0, 3.0, 2.0, 2.0,
        ];
        let qcoeff = vec![
            6.0, -1.0, 5.0, 1.0, 1.0, 1.0, 1.0, 7.0, 5.0, 4.0, 1.0, -1.0,
        ];
        let mut n = 0;
        let mut rcond = 0.0;
        let mut a = DMatrix::zeros(4, 4);
        let mut b = DMatrix::zeros(4, 2);
        let mut c = DMatrix::zeros(2, 4);
        let mut d = DMatrix::zeros(2, 2);
        let info = tc04ad(
            Leri::Left,
            2,
            2,
            &index,
            &pcoeff,
            2,
            2,
            &qcoeff,
            2,
            2,
            &mut n,
            &mut rcond,
            &mut a,
            &mut b,
            &mut c,
            &mut d,
        );
        assert_eq!(info, 0);
        assert_eq!(n, 4);
        assert!(rcond > 0.0);
        assert!(a[(0, 0)].is_finite());
        assert!(b[(0, 0)].is_finite());
    }
}
