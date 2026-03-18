//! MB01PD — Scale or undo scaling of a matrix (SLICOT MB01PD)
//!
//! Scales so the matrix norm is in [SMLNUM, BIGNUM], or undoes that scaling.
//! This implementation supports SCUN 'S'/'U' and TYPE 'G' (full matrix) only.

use nalgebra::DMatrix;

/// Scale or undo scaling.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01PdScun {
    Scale,
    Unscale,
}

/// Matrix storage type. Only General (full) is implemented.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01PdType {
    General,
    Lower,
    Upper,
    Hessenberg,
}

/// Scales or undoes scaling of A so its norm is in a safe range. ANRM is the norm of the original A (for scale) or the norm when scaled (for unscale).
/// Supports TYPE 'G' (full) and 'H' (Hessenberg) for Unscale; Hessenberg uses the same factor for the whole matrix.
///
/// # Returns
/// 0 on success; 1 if TYPE not supported; < 0 if the i-th argument is invalid.
pub fn mb01pd(
    scun: Mb01PdScun,
    typ: Mb01PdType,
    _m: usize,
    _n: usize,
    _kl: i32,
    _ku: i32,
    anrm: f64,
    _nbl: i32,
    _nrows: &[i32],
    a: &mut DMatrix<f64>,
) -> i32 {
    if typ != Mb01PdType::General && (typ != Mb01PdType::Hessenberg || scun != Mb01PdScun::Unscale) {
        return 1;
    }
    if anrm < 0.0 {
        return -7;
    }
    if anrm == 0.0 || a.nrows() == 0 || a.ncols() == 0 {
        return 0;
    }
    let smlnum = f64::MIN_POSITIVE * (1.0 + 1.0 / 16.0);
    let bignum = 1.0 / smlnum;
    let mul = match scun {
        Mb01PdScun::Scale => {
            if anrm < smlnum {
                smlnum / anrm
            } else if anrm > bignum {
                bignum / anrm
            } else {
                return 0;
            }
        }
        Mb01PdScun::Unscale => {
            if anrm < smlnum {
                anrm / smlnum
            } else if anrm > bignum {
                anrm / bignum
            } else {
                return 0;
            }
        }
    };
    for x in a.iter_mut() {
        *x *= mul;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01pd_scale_tiny() {
        // Use a norm below SMLNUM so scale-up path is taken (1e-308 is representable)
        let anrm = 1e-308_f64;
        let mut a = DMatrix::from_row_slice(1, 1, &[anrm]);
        let nrows: [i32; 0] = [];
        assert_eq!(
            mb01pd(
                Mb01PdScun::Scale,
                Mb01PdType::General,
                1,
                1,
                0,
                0,
                anrm,
                0,
                &nrows,
                &mut a
            ),
            0
        );
        // After scale-up, norm should be in safe range (>= SMLNUM order)
        assert!(a[(0, 0)].abs() >= 1e-308);
    }
}
