//! MC01SW — Scale real polynomial (SLICOT MC01SW)
//
// P(x) := alpha * P(x). Coefficients in increasing powers.

/// In-place scale by alpha.
pub fn mc01sw(dp: i32, alpha: f64, p: &mut [f64]) -> i32 {
    let dp = dp as usize;
    if dp + 1 > p.len() {
        return -3;
    }
    for i in 0..=dp {
        p[i] *= alpha;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01sw() {
        let mut p = [1.0, 2.0, 3.0];
        assert_eq!(mc01sw(2, 2.0, &mut p), 0);
        assert!((p[0] - 2.0).abs() < 1e-10);
        assert!((p[1] - 4.0).abs() < 1e-10);
    }
}
