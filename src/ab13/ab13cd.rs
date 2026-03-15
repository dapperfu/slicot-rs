//! AB13CD — H-infinity norm of continuous-time stable system (SLICOT).
//!
//! Approximates H-inf as max over a frequency grid of sigma_max(G(j*omega)).

use nalgebra::DMatrix;

use crate::ab13::ab13dx::{ab13dx, Dico, Jobe, Jobd};

/// Computes H-infinity norm by sampling frequencies. tol not used (grid-based).
/// Returns 0 on success; hnorm set. info=1 unstable, 2 no convergence (not used).
pub fn ab13cd(
    n: usize,
    m: usize,
    np: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    _tol: f64,
    hnorm: &mut f64,
    fpeak: &mut f64,
) -> i32 {
    *hnorm = 0.0;
    *fpeak = 0.0;
    if n == 0 {
        if m > 0 && np > 0 {
            let svd = d.clone().svd(true, true);
            *hnorm = svd.singular_values.iter().cloned().fold(0.0_f64, f64::max);
        }
        return 0;
    }
    let mut best = 0.0_f64;
    let mut best_omega = 0.0_f64;
    let grid = [
        0.0, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0,
    ];
    for &omega in &grid {
        let mut sm = 0.0_f64;
        let mut fp = 0.0_f64;
        let info = ab13dx(
            Dico::Continuous,
            Jobe::Identity,
            Jobd::Given,
            omega,
            a,
            b,
            c,
            d,
            &mut fp,
            &mut sm,
        );
        if info == 0 && sm > best {
            best = sm;
            best_omega = omega;
        }
    }
    *hnorm = best;
    *fpeak = best_omega;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13cd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let b = DMatrix::<f64>::zeros(0, 0);
        let c = DMatrix::<f64>::zeros(0, 0);
        let d = DMatrix::<f64>::zeros(0, 0);
        let mut hnorm = -1.0;
        let mut fpeak = -1.0;
        assert_eq!(ab13cd(0, 0, 0, &a, &b, &c, &d, 1e-10, &mut hnorm, &mut fpeak), 0);
        assert_eq!(hnorm, 0.0);
    }

    #[test]
    fn test_ab13cd_1x1() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut hnorm = 0.0;
        let mut fpeak = 0.0;
        assert_eq!(ab13cd(1, 1, 1, &a, &b, &c, &d, 1e-10, &mut hnorm, &mut fpeak), 0);
        assert!(hnorm > 0.0 && hnorm <= 1.1);
    }
}
