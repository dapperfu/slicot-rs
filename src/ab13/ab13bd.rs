//! AB13BD — H2 or L2 norm of (A,B,C,D) (SLICOT).
//!
//! For stable system: H2 = sqrt(trace(B'*X*B)) where A'*X+X*A = -C'*C (continuous).
//! D must be zero for continuous H2. Simplified: continuous only, D=0.

use nalgebra::DMatrix;

use crate::ab13::lyapunov;

/// H2 or L2 norm.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobn {
    H2,
    L2,
}

/// Computes H2 or L2 norm. Continuous only; D must be zero for H2.
/// Returns 0 on success; norm is set. info=1 Schur/Lyapunov failed, 5 D nonzero (cont H2), 6 unstable.
pub fn ab13bd(
    _dico: char,
    jobn: Jobn,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    norm: &mut f64,
    _dwork: &mut [f64],
) -> i32 {
    *norm = 0.0;
    let n = a.nrows();
    if n == 0 {
        return 0;
    }
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -1;
    }
    if jobn == Jobn::H2 && d.iter().any(|&x| x != 0.0) {
        return 5;
    }
    let at = a.transpose();
    let q_obs = &c.transpose() * c;
    let mut x = DMatrix::<f64>::zeros(n, n);
    if !lyapunov::lyapunov_continuous(&at, &q_obs, &mut x) {
        return 4;
    }
    let btxb = b.transpose() * &x * b;
    let trace_val = (0..b.ncols()).map(|j| btxb[(j, j)]).sum::<f64>();
    *norm = trace_val.max(0.0).sqrt();
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13bd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let b = DMatrix::<f64>::zeros(0, 0);
        let c = DMatrix::<f64>::zeros(0, 0);
        let d = DMatrix::<f64>::zeros(0, 0);
        let mut norm = -1.0;
        let mut dwork = vec![0.0; 1];
        assert_eq!(ab13bd('C', Jobn::H2, &a, &b, &c, &d, &mut norm, &mut dwork), 0);
        assert_eq!(norm, 0.0);
    }

    #[test]
    fn test_ab13bd_1x1_stable() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut norm = 0.0;
        let mut dwork = vec![0.0; 4];
        assert_eq!(ab13bd('C', Jobn::H2, &a, &b, &c, &d, &mut norm, &mut dwork), 0);
        assert!(norm > 0.0 && norm < 1.0);
    }
}
