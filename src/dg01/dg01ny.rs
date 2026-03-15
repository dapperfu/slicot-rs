//! DG01NY — Auxiliary for real FFT: combines/splits real signal form with complex FFT (SLICOT DG01NY).
//!
//! Called by DG01ND. No parameter checks; N must be power of 2, N >= 2. XR, XI length >= N+1.

use std::f64::consts::PI;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dg01NyIndi {
    Direct,
    Inverse,
}

/// In-place conversion between real signal representation and complex FFT form.
pub fn dg01ny(indi: Dg01NyIndi, n: usize, xr: &mut [f64], xi: &mut [f64]) {
    if n < 2 || xr.len() < n + 1 || xi.len() < n + 1 {
        return;
    }
    let lindi = indi == Dg01NyIndi::Direct;
    let pi2 = if lindi { -2.0 * PI } else { 2.0 * PI };
    let whelp = pi2 / (2 * n) as f64;
    let wstpi = whelp.sin();
    let wstpr = -2.0 * (whelp / 2.0).sin().powi(2);
    let mut wi = 0.0_f64;
    let mut wr = if lindi { 1.0 } else { -1.0 };
    if lindi {
        xr[n] = xr[0];
        xi[n] = xi[0];
    }
    let n2 = n / 2 + 1;
    for i in 1..=n2 {
        let j = n + 2 - i;
        if j < i {
            break;
        }
        let i0 = i - 1;
        let j0 = j - 1;
        let mut ar = xr[i0] + xr[j0];
        let mut ai = xi[i0] - xi[j0];
        let mut br = xi[i0] + xi[j0];
        let mut bi = xr[j0] - xr[i0];
        if lindi {
            ar *= 0.5;
            ai *= 0.5;
            br *= 0.5;
            bi *= 0.5;
        }
        let helpr = wr * br - wi * bi;
        let helpi = wr * bi + wi * br;
        xr[i0] = ar + helpr;
        xi[i0] = ai + helpi;
        xr[j0] = ar - helpr;
        xi[j0] = helpi - ai;
        let whelp_ = wr;
        wr = wr + wr * wstpr - wi * wstpi;
        wi = wi + wi * wstpr + whelp_ * wstpi;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dg01ny_direct_n2() {
        let mut xr = vec![1.0, 2.0, 1.0];
        let mut xi = vec![0.0, 0.0, 0.0];
        dg01ny(Dg01NyIndi::Direct, 2, &mut xr, &mut xi);
        assert_eq!(xr.len(), 3);
        assert_eq!(xi.len(), 3);
    }

    #[test]
    fn test_dg01ny_inverse_n2() {
        let mut xr = vec![1.0, 2.0, 1.0];
        let mut xi = vec![0.0, 0.0, 0.0];
        dg01ny(Dg01NyIndi::Inverse, 2, &mut xr, &mut xi);
        assert_eq!(xr.len(), 3);
        assert_eq!(xi.len(), 3);
    }

    #[test]
    fn test_dg01ny_n_too_small() {
        let mut xr = vec![1.0];
        let mut xi = vec![0.0];
        dg01ny(Dg01NyIndi::Direct, 1, &mut xr, &mut xi);
        assert_eq!(xr[0], 1.0);
    }
}
