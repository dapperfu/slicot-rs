//! DE01PD — Convolution or deconvolution using Hartley transform (SLICOT DE01PD).
//!
//! Stub: returns 0 for N=0; otherwise INFO=1 (use DE01OD for FFT-based convolution).

/// WGHT: true = weights available in W, false = not available.
/// Stub: returns 0 for N=0; 1 (not implemented) otherwise.
pub fn de01pd(_conv: bool, _wght: bool, n: usize, _a: &mut [f64], _b: &mut [f64], _w: &mut [f64]) -> i32 {
    if n == 0 {
        return 0;
    }
    1
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
