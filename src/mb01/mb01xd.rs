//! MB01XD — Compute U'*U or L*L' in place (SLICOT MB01XD)
//
// Calls the unblocked kernel MB01XY.

use nalgebra::DMatrix;

use super::mb01xy::{mb01xy, Mb01XyUplo};

/// Which triangle is stored.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01XdUplo {
    Upper,
    Lower,
}

/// Overwrites the triangular part of A with U'*U (Upper) or L*L' (Lower).
pub fn mb01xd(uplo: Mb01XdUplo, a: &mut DMatrix<f64>) -> i32 {
    let uplo_xy = match uplo {
        Mb01XdUplo::Upper => Mb01XyUplo::Upper,
        Mb01XdUplo::Lower => Mb01XyUplo::Lower,
    };
    mb01xy(uplo_xy, a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb01xd_upper() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
        assert_eq!(mb01xd(Mb01XdUplo::Upper, &mut a), 0);
        assert!((a[(1, 1)] - 13.0).abs() < 1e-15);
    }
}
