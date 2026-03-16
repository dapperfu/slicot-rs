//! TD05AD — Evaluation of transfer function G(jW) at a specified frequency (SLICOT TD05AD)
//!
//! G(jW) = (B(1)+B(2)*(jW)+...)/(A(1)+A(2)*(jW)+...). Output Cartesian (real, imag) or polar (dB, degrees).

/// Frequency unit.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum UnitF {
    /// Radians per second.
    Radians,
    /// Hertz (W = 2*pi*f).
    Hertz,
}

/// Output format.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Output {
    /// Cartesian: (real part, imaginary part).
    Cartesian,
    /// Polar: (magnitude in dB, phase in degrees).
    Polar,
}

/// Evaluates G(jW) for denominator A and numerator B (coefficients in ascending powers of jW).
///
/// # Arguments
/// * `unit_f` - Radians or Hertz (W = 2*pi*f if Hertz)
/// * `output` - Cartesian (VALR=Re(G), VALI=Im(G)) or Polar (VALR=|G| dB, VALI=phase deg)
/// * `a` - Denominator coefficients: A(i) = coefficient of (jW)^(i-1), length NP1 = N+1
/// * `b` - Numerator coefficients: B(i) = coefficient of (jW)^(i-1), length MP1 = M+1
/// * `w` - Frequency value (rad/s or Hz per unit_f)
/// * `valr` - Output: real part or magnitude (dB)
/// * `vali` - Output: imaginary part or phase (degrees)
///
/// # Returns
/// * `0` - success
/// * `< 0` - invalid argument
/// * `1` - W is a pole (denominator zero) or all A coefficients zero
pub fn td05ad(
    unit_f: UnitF,
    output: Output,
    a: &[f64],
    b: &[f64],
    w: f64,
    valr: &mut f64,
    vali: &mut f64,
) -> i32 {
    let np1 = a.len();
    let mp1 = b.len();
    if np1 < 1 || mp1 < 1 {
        return -3;
    }

    let omega = match unit_f {
        UnitF::Radians => w,
        UnitF::Hertz => w * 2.0 * std::f64::consts::PI,
    };

    // (jW)^k: k=0->1, k=1->j*W, k=2->-W^2, k=3->-j*W^3, k=4->W^4, ...
    let mut den_re = 0.0;
    let mut den_im = 0.0;
    let mut w_pow = 1.0;
    for k in 0..np1 {
        let (re, im) = match k % 4 {
            0 => (1.0, 0.0),
            1 => (0.0, 1.0),
            2 => (-1.0, 0.0),
            _ => (0.0, -1.0),
        };
        den_re += a[k] * re * w_pow;
        den_im += a[k] * im * w_pow;
        w_pow *= omega;
    }

    let den_norm_sq = den_re * den_re + den_im * den_im;
    if den_norm_sq <= 0.0 {
        return 1;
    }
    // Check for near-zero leading coefficient (all A zero would make den = 0)
    let a_leading = a[np1 - 1];
    if a_leading.abs() < 1e-100 && np1 == 1 {
        return 1;
    }

    let mut num_re = 0.0;
    let mut num_im = 0.0;
    let mut w_pow = 1.0;
    for k in 0..mp1 {
        let (re, im) = match k % 4 {
            0 => (1.0, 0.0),
            1 => (0.0, 1.0),
            2 => (-1.0, 0.0),
            _ => (0.0, -1.0),
        };
        num_re += b[k] * re * w_pow;
        num_im += b[k] * im * w_pow;
        w_pow *= omega;
    }

    // G = num / den
    let g_re = (num_re * den_re + num_im * den_im) / den_norm_sq;
    let g_im = (num_im * den_re - num_re * den_im) / den_norm_sq;

    match output {
        Output::Cartesian => {
            *valr = g_re;
            *vali = g_im;
        }
        Output::Polar => {
            let mag_sq = g_re * g_re + g_im * g_im;
            *valr = if mag_sq > 0.0 {
                10.0 * mag_sq.log10()
            } else {
                -1e30 // -inf in practice
            };
            *vali = g_im.atan2(g_re).to_degrees();
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_td05ad_cartesian_radians() {
        // SLICOT example: NP1=6, MP1=4, W=1, R, C
        // A = 1,1,0,0,2,1  B = 6,2,3,1
        let a = vec![1.0, 1.0, 0.0, 0.0, 2.0, 1.0];
        let b = vec![6.0, 2.0, 3.0, 1.0];
        let mut valr = 0.0;
        let mut vali = 0.0;
        assert_eq!(
            td05ad(UnitF::Radians, Output::Cartesian, &a, &b, 1.0, &mut valr, &mut vali),
            0
        );
        // Expected: 0.8462 -0.2308*j
        assert!((valr - 0.8462).abs() < 0.001);
        assert!((vali - (-0.2308)).abs() < 0.001);
    }

    #[test]
    fn test_td05ad_polar() {
        let a = vec![1.0, 0.0, 1.0]; // 1 + 0*j*w - w^2 = 1 - w^2 at w=0
        let b = vec![1.0, 0.0];      // 1
        let mut valr = 0.0;
        let mut vali = 0.0;
        assert_eq!(
            td05ad(UnitF::Radians, Output::Polar, &a, &b, 0.0, &mut valr, &mut vali),
            0
        );
        // G(0) = 1/1 = 1 -> 0 dB, 0 deg
        assert!((valr - 0.0).abs() < 1e-6);
        assert!((vali - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_td05ad_hertz() {
        // W in Hz -> omega = 2*pi*W. At W=1/(2*pi), omega=1.
        let a = vec![1.0, 1.0];
        let b = vec![1.0];
        let mut valr = 0.0;
        let mut vali = 0.0;
        let f = 1.0 / (2.0 * std::f64::consts::PI);
        assert_eq!(
            td05ad(UnitF::Hertz, Output::Cartesian, &a, &b, f, &mut valr, &mut vali),
            0
        );
        // G(j*1) = 1/(1+j) = (1-j)/2
        assert!((valr - 0.5).abs() < 1e-10);
        assert!((vali - (-0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_td05ad_pole() {
        let a = vec![0.0, 0.0]; // zero denominator
        let b = vec![1.0];
        let mut valr = 0.0;
        let mut vali = 0.0;
        assert_eq!(
            td05ad(UnitF::Radians, Output::Cartesian, &a, &b, 1.0, &mut valr, &mut vali),
            1
        );
    }
}
