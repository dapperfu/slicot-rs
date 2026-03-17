//! SB03OV — Complex plane rotation (SLICOT SB03OV).
//!
//! Constructs a complex plane rotation so that for complex a and real b,
//! ( conjg(c)  s; -s c ) * (a; b) = (d; 0) with d real. On return a(1)=d, a(2)=0.

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// DLAPY3(x,y,z) = sqrt(x^2 + y^2 + z^2) avoiding overflow. Public for SB03OY.
#[inline]
pub(crate) fn dlapy3(x: f64, y: f64, z: f64) -> f64 {
    let x_abs = x.abs();
    let y_abs = y.abs();
    let z_abs = z.abs();
    let max_abs = x_abs.max(y_abs).max(z_abs);
    if max_abs == ZERO {
        return ZERO;
    }
    let scale = max_abs;
    let x_s = x / scale;
    let y_s = y / scale;
    let z_s = z / scale;
    scale * (x_s * x_s + y_s * y_s + z_s * z_s).sqrt()
}

/// SB03OV: complex plane rotation. A = [real(a), imag(a)] on entry; on exit A = [d, 0].
/// C = [real(c), imag(c)], S = sine. If norm([a;b]) < small, rotation is identity.
pub fn sb03ov(a: &mut [f64; 2], b: f64, small: f64, c: &mut [f64; 2], s: &mut f64) {
    let d = dlapy3(a[0], a[1], b);
    if d < small {
        c[0] = ONE;
        c[1] = ZERO;
        *s = ZERO;
        if d > ZERO {
            a[0] = d;
            a[1] = ZERO;
        }
    } else {
        c[0] = a[0] / d;
        c[1] = a[1] / d;
        *s = b / d;
        a[0] = d;
        a[1] = ZERO;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dlapy3() {
        assert!((dlapy3(3.0, 4.0, 0.0) - 5.0).abs() < 1e-10);
        assert!((dlapy3(0.0, 0.0, 0.0) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_sb03ov_zero() {
        let mut a = [1e-20, 0.0];
        let mut c = [0.0; 2];
        let mut s = 0.0;
        sb03ov(&mut a, 0.0, 1e-10, &mut c, &mut s);
        assert_eq!(c[0], 1.0);
        assert_eq!(c[1], 0.0);
        assert_eq!(s, 0.0);
    }

    #[test]
    fn test_sb03ov_real() {
        let mut a = [3.0, 0.0];
        let mut c = [0.0; 2];
        let mut s = 0.0;
        sb03ov(&mut a, 4.0, 1e-10, &mut c, &mut s);
        assert!((a[0] - 5.0).abs() < 1e-10);
        assert_eq!(a[1], 0.0);
        assert!((c[0] - 3.0 / 5.0).abs() < 1e-10);
        assert_eq!(c[1], 0.0);
        assert!((s - 4.0 / 5.0).abs() < 1e-10);
    }
}
