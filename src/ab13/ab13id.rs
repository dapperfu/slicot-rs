//! AB13ID — Test if descriptor transfer function G(lambda) = C*(lambda*E-A)^{-1}*B is proper (SLICOT).
//!
//! Returns true if proper (no infinite poles), false otherwise. Simplified: N=0 => true;
//! if E is full column rank (rank = N) then proper; if E is singular (rank < N) return false.

use nalgebra::DMatrix;

/// Default tolerance for rank determination (relative to largest singular value).
pub const DEFAULT_TOL: f64 = 1e-10;

/// Tests whether the descriptor system G(lambda) = C*(lambda*E - A)^{-1}*B is proper.
///
/// A descriptor system is proper if E has full rank N (no infinite poles). Uses SVD of E
/// to compute rank. Only E is required for this simplified test.
///
/// # Arguments
/// * `e` - N×N descriptor matrix E
/// * `tol` - tolerance for rank: singular values > tol * sigma_max are counted (use DEFAULT_TOL if unsure)
///
/// # Returns
/// true if proper (N=0 or rank(E) = N), false if improper (rank(E) < N).
pub fn ab13id(e: &DMatrix<f64>, tol: f64) -> bool {
    let n = e.nrows();
    if n == 0 {
        return true;
    }
    if e.ncols() != n {
        return false;
    }
    let svd = e.clone().svd(true, false);
    let sigmas = &svd.singular_values;
    let sigma_max = sigmas.iter().cloned().fold(0.0_f64, f64::max);
    let threshold = tol * sigma_max.max(1.0) + tol; // avoid 0*eps
    let rank = sigmas.iter().filter(|&&s| s > threshold).count();
    rank >= n
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab13id_n0_returns_true() {
        let e = DMatrix::<f64>::zeros(0, 0);
        assert!(ab13id(&e, DEFAULT_TOL));
    }

    #[test]
    fn test_ab13id_e_identity_a_zero_proper() {
        // E = I, A = 0 => proper (full rank E).
        let e = DMatrix::<f64>::identity(2, 2);
        assert!(ab13id(&e, DEFAULT_TOL));
    }

    #[test]
    fn test_ab13id_singular_e_improper() {
        // E = 0 (singular) => improper.
        let e = DMatrix::<f64>::zeros(2, 2);
        assert!(!ab13id(&e, DEFAULT_TOL));
    }
}
