//! SG03BR — Complex Givens rotation in real arithmetic (SLICOT SG03BR)
//!
//! Computes parameters (C, SR, SI, ZR, ZI) for the complex Givens rotation
//! that zeros the second component of (XR+XI*I, YR+YI*I)^T.
//! Adaptation for real data of LAPACK ZLARTG; avoids overflow.

/// Computes parameters for the complex Givens rotation
///
/// ```text
/// (   C     SR+SI*I )   ( XR+XI*I )   ( ZR+ZI*I )
/// (                 ) * (         ) = (         )
/// ( -SR+SI*I   C    )   ( YR+YI*I )   (    0    )
/// ```
/// with C^2 + |SR+SI*I|^2 = 1.
///
/// # Arguments
/// * `xr`, `xi` — real and imaginary parts of first component
/// * `yr`, `yi` — real and imaginary parts of second component
///
/// # Outputs (via mutable references)
/// * `c`, `sr`, `si` — rotation parameters
/// * `zr`, `zi` — real and imaginary parts of the first component of the result
#[allow(clippy::too_many_arguments)]
pub fn sg03br(
    xr: f64,
    xi: f64,
    yr: f64,
    yi: f64,
    c: &mut f64,
    sr: &mut f64,
    si: &mut f64,
    zr: &mut f64,
    zi: &mut f64,
) {
    const ONE: f64 = 1.0;
    const TWO: f64 = 2.0;
    const ZERO: f64 = 0.0;

    let safmin = f64::MIN_POSITIVE;
    let eps = f64::EPSILON;
    let base = 2.0_f64;
    let exp = (safmin / eps).ln() / base.ln() / TWO;
    let safmn2 = base.powf(exp.floor());
    let safmx2 = ONE / safmn2;

    fn dlapy2(x: f64, y: f64) -> f64 {
        x.hypot(y)
    }

    let mut scale = xr.abs().max(xi.abs()).max(yr.abs()).max(yi.abs());
    let mut xrs = xr;
    let mut xis = xi;
    let mut yrs = yr;
    let mut yis = yi;
    let mut count = 0_i32;

    if scale >= safmx2 {
        loop {
            count += 1;
            xrs *= safmn2;
            xis *= safmn2;
            yrs *= safmn2;
            yis *= safmn2;
            scale *= safmn2;
            if scale < safmx2 {
                break;
            }
        }
    } else if scale <= safmn2 {
        if yr == ZERO && yi == ZERO {
            *c = ONE;
            *sr = ZERO;
            *si = ZERO;
            *zr = xr;
            *zi = xi;
            return;
        }
        loop {
            count -= 1;
            xrs *= safmx2;
            xis *= safmx2;
            yrs *= safmx2;
            yis *= safmx2;
            scale *= safmx2;
            if scale > safmn2 {
                break;
            }
        }
    }

    let x2 = xrs * xrs + xis * xis;
    let y2 = yrs * yrs + yis * yis;

    if x2 <= y2.max(ONE) * safmin {
        if xr == ZERO && xi == ZERO {
            *c = ZERO;
            *zr = dlapy2(yr, yi);
            *zi = ZERO;
            let d = dlapy2(yrs, yis);
            *sr = yrs / d;
            *si = -yis / d;
            return;
        }
        let x2s = dlapy2(xrs, xis);
        let y2s = y2.sqrt();
        *c = x2s / y2s;
        let (tr, ti) = if xr.abs().max(xi.abs()) > ONE {
            let d = dlapy2(xr, xi);
            (xr / d, xi / d)
        } else {
            let dr = safmx2 * xr;
            let di = safmx2 * xi;
            let d = dlapy2(dr, di);
            (dr / d, di / d)
        };
        *sr = tr * (yrs / y2s) + ti * (yis / y2s);
        *si = ti * (yrs / y2s) - tr * (yis / y2s);
        *zr = *c * xr + *sr * yr - *si * yi;
        *zi = *c * xi + *si * yr + *sr * yi;
        return;
    }

    let x2s = (ONE + y2 / x2).sqrt();
    *zr = x2s * xrs;
    *zi = x2s * xis;
    *c = ONE / x2s;
    let d = x2 + y2;
    let mut sr_val = *zr / d;
    let mut si_val = *zi / d;
    let dr = sr_val * yrs + si_val * yis;
    si_val = si_val * yrs - sr_val * yis;
    sr_val = dr;
    *sr = sr_val;
    *si = si_val;

    if count != 0 {
        if count > 0 {
            for _ in 0..count {
                *zr *= safmx2;
                *zi *= safmx2;
            }
        } else {
            for _ in 0..(-count) {
                *zr *= safmn2;
                *zi *= safmn2;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sg03br_zeros_y() {
        let (xr, xi) = (3.0, 4.0);
        let (yr, yi) = (0.0, 0.0);
        let mut c = 0.0;
        let mut sr = 0.0;
        let mut si = 0.0;
        let mut zr = 0.0;
        let mut zi = 0.0;
        sg03br(xr, xi, yr, yi, &mut c, &mut sr, &mut si, &mut zr, &mut zi);
        assert_eq!(c, 1.0);
        assert_eq!(sr, 0.0);
        assert_eq!(si, 0.0);
        assert!((zr - 3.0).abs() < 1e-10);
        assert!((zi - 4.0).abs() < 1e-10);
    }

    #[test]
    fn sg03br_rotation_norm() {
        let (xr, xi) = (1.0, 0.0);
        let (yr, yi) = (0.0, 1.0);
        let mut c = 0.0;
        let mut sr = 0.0;
        let mut si = 0.0;
        let mut zr = 0.0;
        let mut zi = 0.0;
        sg03br(xr, xi, yr, yi, &mut c, &mut sr, &mut si, &mut zr, &mut zi);
        assert!((c * c + sr * sr + si * si - 1.0).abs() < 1e-10);
        assert!((zr * yr + zi * yi - (-sr * xr + si * xi)).abs() < 1e-10);
    }
}
