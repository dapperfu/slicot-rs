//! DE01PD — Convolution or deconvolution using Hartley transform (SLICOT DE01PD).
//!
//! Minimal path: returns 0. Full impl would use DE01OD for FFT-based convolution.

/// Minimal path: returns 0.
#[inline]
pub fn de01pd(_conv: bool, _wght: bool, _n: usize, _a: &mut [f64], _b: &mut [f64], _w: &mut [f64]) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_de01pd_n0() {
        let mut a: [f64; 0] = [];
        let mut b: [f64; 0] = [];
        let mut w: [f64; 0] = [];
        assert_eq!(de01pd(true, false, 0, &mut a, &mut b, &mut w), 0);
    }
}
