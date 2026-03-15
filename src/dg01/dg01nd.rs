//! DG01ND — Discrete Fourier transform of a real signal (SLICOT DG01ND).
//!
//! Uses DG01MD and DG01NY. N is half the number of real samples (so 2*N samples). XR, XI length >= N+1.

use crate::dg01::dg01md::{dg01md, Dg01MdIndi};
use crate::dg01::dg01ny::{dg01ny, Dg01NyIndi};

fn is_power_of_two(n: usize) -> bool {
    n >= 2 && (n & (n - 1)) == 0
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dg01NdIndi {
    Direct,
    Inverse,
}

/// Real FFT or inverse. XR(0..N) and XI(0..N) on entry: for Direct, odd/even parts of real signal;
/// on exit: real and imaginary parts of FFT. Lengths must be >= N+1.
pub fn dg01nd(indi: Dg01NdIndi, n: usize, xr: &mut [f64], xi: &mut [f64]) -> i32 {
    if !matches!(indi, Dg01NdIndi::Direct | Dg01NdIndi::Inverse) {
        return -1;
    }
    if !is_power_of_two(n) {
        return -2;
    }
    if xr.len() < n + 1 || xi.len() < n + 1 {
        return -3;
    }
    let lindi = indi == Dg01NdIndi::Direct;
    if !lindi {
        dg01ny(Dg01NyIndi::Inverse, n, xr, xi);
    }
    let info = dg01md(
        if lindi { Dg01MdIndi::Direct } else { Dg01MdIndi::Inverse },
        n,
        &mut xr[..n],
        &mut xi[..n],
    );
    if info != 0 {
        return info;
    }
    if lindi {
        dg01ny(Dg01NyIndi::Direct, n, xr, xi);
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dg01nd_n2() {
        let mut xr = [1.0, 0.0, 0.0];
        let mut xi = [0.0, 0.0, 0.0];
        assert_eq!(dg01nd(Dg01NdIndi::Direct, 2, &mut xr, &mut xi), 0);
        assert_eq!(dg01nd(Dg01NdIndi::Inverse, 2, &mut xr, &mut xi), 0);
    }
}
