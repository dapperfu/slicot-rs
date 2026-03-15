//! DG01OD — Discrete Hartley transform (SLICOT DG01OD).
//!
//! Stub: returns 0 for N=0 or N=1; INFO=1 otherwise.

/// Stub: returns 0 for N=0 or N=1; 1 (not implemented) otherwise.
pub fn dg01od(_scr: u8, _wght: u8, n: usize, _a: &mut [f64], _w: &mut [f64]) -> i32 {
    if n == 0 || n == 1 {
        return 0;
    }
    1
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
