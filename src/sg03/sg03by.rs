//! SG03BY — Complex plane rotation in real arithmetic (SLICOT SG03BY)
//!
//! Computes parameters (CR, CI, SR, SI, Z) for the complex Givens rotation
//! such that the second component of the result is zero; Z is non-negative real.

/// Computes parameters for the complex Givens rotation
///
/// ```text
/// ( CR-CI*I   SR-SI*I )   ( XR+XI*I )   ( Z )
/// (                    ) * (         ) = (   )
/// ( -SR-SI*I  CR+CI*I )   ( YR+YI*I )   ( 0 )
/// ```
/// Z is a non-negative real number.
///
/// # Arguments
/// * `xr`, `xi` — real and imaginary parts of first component
/// * `yr`, `yi` — real and imaginary parts of second component
///
/// # Outputs (via mutable references)
/// * `cr`, `ci` — first column of rotation (complex)
/// * `sr`, `si` — second column of rotation (complex)
/// * `z` — non-negative real (norm of result)
pub fn sg03by(
    xr: f64,
    xi: f64,
    yr: f64,
    yi: f64,
    cr: &mut f64,
    ci: &mut f64,
    sr: &mut f64,
    si: &mut f64,
    z: &mut f64,
) {
    const ONE: f64 = 1.0;
    const ZERO: f64 = 0.0;

    *z = xr.abs().max(xi.abs()).max(yr.abs()).max(yi.abs());

    if *z == ZERO {
        *cr = ONE;
        *ci = ZERO;
        *sr = ZERO;
        *si = ZERO;
    } else {
        let xr_z = xr / *z;
        let xi_z = xi / *z;
        let yr_z = yr / *z;
        let yi_z = yi / *z;
        *z *= (xr_z * xr_z + xi_z * xi_z + yr_z * yr_z + yi_z * yi_z).sqrt();
        *cr = xr / *z;
        *ci = xi / *z;
        *sr = yr / *z;
        *si = yi / *z;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sg03by_zeros() {
        let mut cr = 0.0;
        let mut ci = 0.0;
        let mut sr = 0.0;
        let mut si = 0.0;
        let mut z = 0.0;
        sg03by(0.0, 0.0, 0.0, 0.0, &mut cr, &mut ci, &mut sr, &mut si, &mut z);
        assert_eq!(cr, 1.0);
        assert_eq!(ci, 0.0);
        assert_eq!(sr, 0.0);
        assert_eq!(si, 0.0);
        assert_eq!(z, 0.0);
    }

    #[test]
    fn sg03by_unit() {
        let (xr, xi) = (1.0, 0.0);
        let (yr, yi) = (0.0, 0.0);
        let mut cr = 0.0;
        let mut ci = 0.0;
        let mut sr = 0.0;
        let mut si = 0.0;
        let mut z = 0.0;
        sg03by(xr, xi, yr, yi, &mut cr, &mut ci, &mut sr, &mut si, &mut z);
        assert!((cr - 1.0).abs() < 1e-10);
        assert!(ci.abs() < 1e-10);
        assert!(sr.abs() < 1e-10);
        assert!(si.abs() < 1e-10);
        assert!((z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sg03by_norm() {
        let (xr, xi) = (3.0, 4.0);
        let (yr, yi) = (1.0, 0.0);
        let mut cr = 0.0;
        let mut ci = 0.0;
        let mut sr = 0.0;
        let mut si = 0.0;
        let mut z = 0.0;
        sg03by(xr, xi, yr, yi, &mut cr, &mut ci, &mut sr, &mut si, &mut z);
        let norm_sq = xr * xr + xi * xi + yr * yr + yi * yi;
        assert!((z * z - norm_sq).abs() < 1e-10);
    }
}
