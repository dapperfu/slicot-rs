//! DG01OD — Discrete Hartley transform (SLICOT DG01OD).
//!
//! Minimal path: returns 0. Full transform not yet implemented.

/// Minimal path: returns 0.
#[inline]
pub fn dg01od(_scr: u8, _wght: u8, _n: usize, _a: &mut [f64], _w: &mut [f64]) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dg01od_n0() {
        let mut a: [f64; 0] = [];
        let mut w: [f64; 0] = [];
        assert_eq!(dg01od(0, 0, 0, &mut a, &mut w), 0);
    }
}
