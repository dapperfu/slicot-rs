//! MC01SX — Real polynomial evaluation at real point (SLICOT MC01SX)
//
// Evaluates P(x) at x = x0 using Horner.

/// Returns P(x0). P in increasing powers.
pub fn mc01sx(dp: i32, x0: f64, p: &[f64], val: &mut f64) -> i32 {
    let dp = dp as usize;
    if dp + 1 > p.len() {
        return -3;
    }
    let mut v = p[dp];
    for j in (0..dp).rev() {
        v = v * x0 + p[j];
    }
    *val = v;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01sx() {
        let p = [1.0, 2.0, 1.0];
        let mut val = 0.0;
        assert_eq!(mc01sx(2, 1.0, &p, &mut val), 0);
        assert!((val - 4.0).abs() < 1e-10);
    }
}
