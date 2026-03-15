//! MC01SD — Derivative of real polynomial (SLICOT MC01SD)
//
// B(x) = d/dx A(x). Coefficients in increasing powers.

/// DP degree of A; A input, B output (length DP for derivative).
pub fn mc01sd(dp: i32, a: &[f64], b: &mut [f64]) -> i32 {
    let dp = dp as usize;
    if dp + 1 > a.len() || dp > b.len() {
        return -2;
    }
    for i in 0..dp {
        b[i] = (i + 1) as f64 * a[i + 1];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01sd() {
        let a = [0.0, 1.0, 2.0];
        let mut b = [0.0; 2];
        assert_eq!(mc01sd(2, &a, &mut b), 0);
        assert!((b[0] - 1.0).abs() < 1e-10);
        assert!((b[1] - 4.0).abs() < 1e-10);
    }
}
