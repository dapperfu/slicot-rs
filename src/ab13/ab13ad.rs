//! AB13AD — Hankel norm of the alpha-stable projection of (A,B,C) (SLICOT).
//!
//! Simplified: if all eigenvalues are alpha-stable, returns the Hankel norm of (A,B,C);
//! otherwise returns error (full additive decomposition not implemented).

use nalgebra::DMatrix;

use crate::ab13::ab13ax::{ab13ax, Dico};

/// Continuous or discrete time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum DicoAd {
    Continuous,
    Discrete,
}

/// Equilibration: not implemented (no scaling).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    No,
    Scale,
}

/// Computes the Hankel norm of the alpha-stable projection. alpha: stability boundary
/// (continuous: real part < alpha; discrete: modulus < alpha).
/// Returns 0 on success; hnorm = Hankel norm; info=1 Schur failed, 2 separation failed (unstable present).
pub fn ab13ad(
    dico: DicoAd,
    _equil: Equil,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    alpha: f64,
    ns: &mut usize,
    hsv: &mut [f64],
    hnorm: &mut f64,
    dwork: &mut [f64],
) -> i32 {
    *hnorm = 0.0;
    *ns = 0;
    let n = a.nrows();
    if n == 0 {
        return 0;
    }
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -1;
    }
    let schur = match a.clone().try_schur(1e-10, 100) {
        Some(s) => s,
        None => return 1,
    };
    let eig = schur.complex_eigenvalues();
    let stable_count = match dico {
        DicoAd::Continuous => eig.iter().filter(|c| c.re < alpha).count(),
        DicoAd::Discrete => eig.iter().filter(|c| c.norm_sqr().sqrt() < alpha).count(),
    };
    if stable_count < n {
        return 2;
    }
    *ns = n;
    let dico_ax = match dico {
        DicoAd::Continuous => Dico::Continuous,
        DicoAd::Discrete => Dico::Discrete,
    };
    let info = ab13ax(dico_ax, a, b, c, hnorm, dwork);
    if info != 0 {
        return if info == 2 { 4 } else { info };
    }
    if hsv.len() >= n {
        hsv[0] = *hnorm;
        for i in 1..n {
            hsv[i] = 0.0;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13ad_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let b = DMatrix::<f64>::zeros(0, 0);
        let c = DMatrix::<f64>::zeros(0, 0);
        let mut ns = 99;
        let mut hsv = vec![0.0; 1];
        let mut hnorm = -1.0;
        let mut dwork = vec![0.0; 1];
        assert_eq!(
            ab13ad(
                DicoAd::Continuous,
                Equil::No,
                &a,
                &b,
                &c,
                0.0,
                &mut ns,
                &mut hsv,
                &mut hnorm,
                &mut dwork,
            ),
            0
        );
        assert_eq!(ns, 0);
        assert_eq!(hnorm, 0.0);
    }

    #[test]
    fn test_ab13ad_all_stable() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut ns = 0;
        let mut hsv = vec![0.0; 2];
        let mut hnorm = 0.0;
        let mut dwork = vec![0.0; 4];
        assert_eq!(
            ab13ad(
                DicoAd::Continuous,
                Equil::No,
                &a,
                &b,
                &c,
                0.0,
                &mut ns,
                &mut hsv,
                &mut hnorm,
                &mut dwork,
            ),
            0
        );
        assert_eq!(ns, 1);
        assert!(hnorm > 0.0);
    }
}
