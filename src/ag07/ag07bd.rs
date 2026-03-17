//! AG07BD — Descriptor system inverse (SLICOT AG07BD).
//!
//! Full SLICOT-equivalent API. Computes inverse (Ai-lambda*Ei,Bi,Ci,Di) of (A-lambda*E,B,C,D).
//! Main path: block assembly per SLICOT formulas (Ai=[A B;C D], Ei, Bi=[0;-I], Ci=[0 I], Di=0).

use nalgebra::DMatrix;

/// Whether E is general or identity.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobE {
    /// E is general.
    General,
    /// E is identity.
    Identity,
}

/// AG07BD: compute inverse of descriptor system.
///
/// # Returns
/// * 0: success (quick return when N=0 or M=0); 1: algorithm not implemented; <0: invalid argument.
pub fn ag07bd(
    jobe: JobE,
    n: usize,
    m: usize,
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    ai: &mut DMatrix<f64>,
    ei: &mut DMatrix<f64>,
    bi: &mut DMatrix<f64>,
    ci: &mut DMatrix<f64>,
    di: &mut DMatrix<f64>,
) -> i32 {
    if n > 0 && (a.nrows() != n || a.ncols() != n) {
        return -4;
    }
    if jobe == JobE::General && n > 0 && (e.nrows() != n || e.ncols() != n) {
        return -6;
    }
    if n > 0 && m > 0 && (b.nrows() != n || b.ncols() != m) {
        return -8;
    }
    if m > 0 && n > 0 && (c.nrows() != m || c.ncols() != n) {
        return -10;
    }
    if m > 0 && (d.nrows() != m || d.ncols() != m) {
        return -12;
    }
    let nm = n + m;
    if ai.nrows() != nm || ai.ncols() != nm {
        return -14;
    }
    if ei.nrows() != nm || ei.ncols() != nm {
        return -16;
    }
    if bi.nrows() != nm || bi.ncols() != m {
        return -18;
    }
    if ci.nrows() != m || ci.ncols() != nm {
        return -20;
    }
    if di.nrows() != m || di.ncols() != m {
        return -22;
    }
    if n == 0 && m == 0 {
        return 0;
    }
    if m == 0 {
        // N>0, M=0: no inputs/outputs; degenerate
        return 0;
    }

    // Main path: Ai = [A B; C D], Ei (block diag E or I, then zeros), Bi = [0; -I], Ci = [0 I], Di = 0.
    if n == 0 {
        // Static system: Ai = D, Ei = I_m, Bi = -I, Ci = I, Di = 0.
        ai.view_mut((0, 0), (m, m)).copy_from(d);
        ei.fill(0.0);
        for i in 0..m {
            ei[(i, i)] = 1.0;
        }
        for i in 0..m {
            for j in 0..m {
                bi[(i, j)] = if i == j { -1.0 } else { 0.0 };
            }
        }
        ci.fill(0.0);
        for i in 0..m {
            ci[(i, i)] = 1.0;
        }
        di.fill(0.0);
        return 0;
    }

    // Ai: top-left N×N = A, top-right N×M = B, bottom-left M×N = C, bottom-right M×M = D
    ai.view_mut((0, 0), (n, n)).copy_from(a);
    ai.view_mut((0, n), (n, m)).copy_from(b);
    ai.view_mut((n, 0), (m, n)).copy_from(c);
    ai.view_mut((n, n), (m, m)).copy_from(d);

    // Ei: top-left N×N = E or I; rest zero (full matrix zeroed then block set)
    ei.fill(0.0);
    if jobe == JobE::Identity {
        for i in 0..n {
            ei[(i, i)] = 1.0;
        }
    } else {
        ei.view_mut((0, 0), (n, n)).copy_from(e);
    }

    // Bi: top N×M = 0, bottom M×M = -I
    bi.fill(0.0);
    for i in 0..m {
        bi[(n + i, i)] = -1.0;
    }

    // Ci: left M×N = 0, right M×M = I
    ci.fill(0.0);
    for i in 0..m {
        ci[(i, n + i)] = 1.0;
    }

    // Di = 0
    di.fill(0.0);

    0
}

/// Compatibility wrapper for benchmarking: (n, m) -> INFO.
#[inline]
pub fn ag07bd_nm(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let nm = n + m;
    let a = DMatrix::zeros(n, n);
    let e = DMatrix::zeros(n, n);
    let b = DMatrix::zeros(n, m);
    let c = DMatrix::zeros(m, n);
    let d = DMatrix::zeros(m, m);
    let mut ai = DMatrix::zeros(nm, nm);
    let mut ei = DMatrix::zeros(nm, nm);
    let mut bi = DMatrix::zeros(nm, m);
    let mut ci = DMatrix::zeros(m, nm);
    let mut di = DMatrix::zeros(m, m);
    ag07bd(
        JobE::Identity,
        n,
        m,
        &a,
        &e,
        &b,
        &c,
        &d,
        &mut ai,
        &mut ei,
        &mut bi,
        &mut ci,
        &mut di,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ag07bd_trivial() {
        let a = DMatrix::zeros(0, 0);
        let e = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 0);
        let c = DMatrix::zeros(0, 0);
        let d = DMatrix::zeros(0, 0);
        let mut ai = DMatrix::zeros(0, 0);
        let mut ei = DMatrix::zeros(0, 0);
        let mut bi = DMatrix::zeros(0, 0);
        let mut ci = DMatrix::zeros(0, 0);
        let mut di = DMatrix::zeros(0, 0);
        assert_eq!(
            ag07bd(JobE::Identity, 0, 0, &a, &e, &b, &c, &d, &mut ai, &mut ei, &mut bi, &mut ci, &mut di),
            0
        );
    }

    #[test]
    fn test_ag07bd_n1_m1_identity_e() {
        let a = DMatrix::from_row_slice(1, 1, &[2.0]);
        let e = DMatrix::zeros(1, 1);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut ai = DMatrix::zeros(2, 2);
        let mut ei = DMatrix::zeros(2, 2);
        let mut bi = DMatrix::zeros(2, 1);
        let mut ci = DMatrix::zeros(1, 2);
        let mut di = DMatrix::zeros(1, 1);
        assert_eq!(
            ag07bd(JobE::Identity, 1, 1, &a, &e, &b, &c, &d, &mut ai, &mut ei, &mut bi, &mut ci, &mut di),
            0
        );
        // Ai = [A B; C D] = [2 1; 1 0]
        assert_eq!(ai[(0, 0)], 2.0);
        assert_eq!(ai[(0, 1)], 1.0);
        assert_eq!(ai[(1, 0)], 1.0);
        assert_eq!(ai[(1, 1)], 0.0);
        // Ei = I in top 1x1, rest 0
        assert_eq!(ei[(0, 0)], 1.0);
        assert_eq!(ei[(1, 0)], 0.0);
        assert_eq!(ei[(0, 1)], 0.0);
        assert_eq!(ei[(1, 1)], 0.0);
        // Bi = [0; -1]
        assert_eq!(bi[(0, 0)], 0.0);
        assert_eq!(bi[(1, 0)], -1.0);
        // Ci = [0 1]
        assert_eq!(ci[(0, 0)], 0.0);
        assert_eq!(ci[(0, 1)], 1.0);
        assert_eq!(di[(0, 0)], 0.0);
    }

    #[test]
    fn test_ag07bd_static_n0_m1() {
        let a = DMatrix::zeros(0, 0);
        let e = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 1);
        let c = DMatrix::zeros(1, 0);
        let d = DMatrix::from_row_slice(1, 1, &[3.0]);
        let mut ai = DMatrix::zeros(1, 1);
        let mut ei = DMatrix::zeros(1, 1);
        let mut bi = DMatrix::zeros(1, 1);
        let mut ci = DMatrix::zeros(1, 1);
        let mut di = DMatrix::zeros(1, 1);
        assert_eq!(
            ag07bd(JobE::Identity, 0, 1, &a, &e, &b, &c, &d, &mut ai, &mut ei, &mut bi, &mut ci, &mut di),
            0
        );
        assert_eq!(ai[(0, 0)], 3.0);
        assert_eq!(ei[(0, 0)], 1.0);
        assert_eq!(bi[(0, 0)], -1.0);
        assert_eq!(ci[(0, 0)], 1.0);
        assert_eq!(di[(0, 0)], 0.0);
    }
}
