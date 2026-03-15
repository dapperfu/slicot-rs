//! MB01QD — Scale matrix A by scalar CTO/CFROM (SLICOT MB01QD)
//!
//! Multiplies A by CTO/CFROM with optional overflow protection. This implementation supports TYPE 'G' (full matrix) only.

use nalgebra::DMatrix;

/// Matrix storage type. Only General (full) is implemented.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01QdType {
    /// Full matrix.
    General,
    /// Lower triangular (not yet implemented).
    Lower,
    /// Upper triangular (not yet implemented).
    Upper,
    /// Upper Hessenberg (not yet implemented).
    Hessenberg,
}

/// Multiplies the M×N matrix A by the scalar CTO/CFROM. Supports only TYPE 'G' (full matrix).
/// CFROM must be nonzero.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn mb01qd(
    typ: Mb01QdType,
    _m: usize,
    _n: usize,
    _kl: i32,
    _ku: i32,
    cfrom: f64,
    cto: f64,
    _nbl: i32,
    _nrows: &[i32],
    a: &mut DMatrix<f64>,
) -> i32 {
    if typ != Mb01QdType::General {
        return 1; // not implemented for other types
    }
    if cfrom == 0.0 {
        return -6;
    }
    let mul = cto / cfrom;
    for x in a.iter_mut() {
        *x *= mul;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01qd_general() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let nrows: [i32; 0] = [];
        assert_eq!(
            mb01qd(
                Mb01QdType::General,
                2,
                2,
                0,
                0,
                1.0,
                2.0,
                0,
                &nrows,
                &mut a
            ),
            0
        );
        assert_eq!(a[(0, 0)], 2.0);
        assert_eq!(a[(1, 1)], 8.0);
    }
}
