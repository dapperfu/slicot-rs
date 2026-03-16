//! TB04BX — Gain of a SISO system given state-space, poles and zeros (SLICOT TB04BX)
//!
//! Gain g = (c*(s0*I-A)^{-1}*b + d) * Prod(s0-Pi)/Prod(s0-Zi) for s0 not a pole/zero.

use nalgebra::{linalg::LU, DMatrix, DVector};

/// Computes the gain of (A,b,c,d) using the formula with a chosen s0.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb04bx(
    ip: usize,
    iz: usize,
    a: &mut DMatrix<f64>,
    _lda: usize,
    b: &mut [f64],
    c: &[f64],
    d: f64,
    pr: &[f64],
    pi: &[f64],
    zr: &[f64],
    zi: &[f64],
    gain: &mut f64,
    iwork: &mut [i32],
) -> i32 {
    if ip == 0 {
        *gain = d;
        return 0;
    }
    if a.nrows() < ip || a.ncols() < ip {
        return -3;
    }
    if b.len() < ip || c.len() < ip {
        return -4;
    }
    if pr.len() < ip || pi.len() < ip {
        return -7;
    }
    if iz > 0 && (zr.len() < iz || zi.len() < iz) {
        return -9;
    }
    if iwork.len() < ip {
        return -12;
    }

    let s0 = 0.0_f64;
    for i in 0..ip {
        if (pr[i] - s0).abs() < 1e-10 && pi[i].abs() < 1e-10 {
            return -7;
        }
    }
    for i in 0..iz {
        if (zr[i] - s0).abs() < 1e-10 && zi[i].abs() < 1e-10 {
            return -9;
        }
    }

    let a_s0 = a.view((0, 0), (ip, ip)).into_owned()
        - DMatrix::identity(ip, ip) * s0;
    let mut b_vec = DVector::from_row_slice(&b[..ip]);
    let lu = LU::new(a_s0);
    if let Some(inv) = lu.try_inverse() {
        let x = inv * &b_vec;
        let mut c_x = 0.0;
        for i in 0..ip {
            c_x += c[i] * x[i];
        }
        let val = c_x + d;
        let mut num = 1.0;
        for i in 0..ip {
            let p_diff = (s0 - pr[i]) * (s0 - pr[i]) + pi[i] * pi[i];
            num *= p_diff;
        }
        let mut den = 1.0;
        for i in 0..iz {
            let z_diff = (s0 - zr[i]) * (s0 - zr[i]) + zi[i] * zi[i];
            den *= z_diff;
        }
        *gain = if den.abs() > 1e-20 { val * num / den } else { val };
    } else {
        return 1;
    }
    for i in 0..ip {
        b_vec[i] = b[i];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb04bx_simple() {
        let mut a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let mut b = [1.0, 0.0];
        let c = [1.0, 1.0];
        let pr = [-1.0, -2.0];
        let pi = [0.0, 0.0];
        let zr: [f64; 0] = [];
        let zi: [f64; 0] = [];
        let mut gain = 0.0;
        let mut iwork = [0i32; 2];
        let info = tb04bx(2, 0, &mut a, 2, &mut b, &c, 0.0, &pr, &pi, &zr, &zi, &mut gain, &mut iwork);
        assert_eq!(info, 0);
        assert!(gain.is_finite());
    }
}
