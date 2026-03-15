//! DF01MD — Sine or cosine transform of a real signal (SLICOT DF01MD).
//!
//! Uses DG01ND. Stub: returns 0 when N=0; INFO=1 otherwise (full impl would require N = 2^k+1, N>=5).

/// Stub: returns 0 when N=0; 1 (not implemented) otherwise.
pub fn df01md(_sico: u8, n: usize, _dt: f64, _a: &mut [f64], _dwork: &mut [f64]) -> i32 {
    if n == 0 {
        return 0;
    }
    1
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
