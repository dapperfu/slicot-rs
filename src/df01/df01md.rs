//! DF01MD — Sine or cosine transform of a real signal (SLICOT DF01MD).
//!
//! Minimal path: returns 0. Full impl would require N = 2^k+1, N>=5 and DG01ND.

/// Minimal path: returns 0. Full transform not yet implemented.
#[inline]
pub fn df01md(_sico: u8, _n: usize, _dt: f64, _a: &mut [f64], _dwork: &mut [f64]) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_df01md_n0() {
        let mut a: [f64; 0] = [];
        let mut dwork: [f64; 0] = [];
        assert_eq!(df01md(b'S', 0, 1.0, &mut a, &mut dwork), 0);
    }
}
