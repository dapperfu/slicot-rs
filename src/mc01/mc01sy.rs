//! MC01SY — Sum of two real polynomials (SLICOT MC01SY)
//
// C(x) = alpha*A(x) + beta*B(x). Coefficients in increasing powers.

/// DP common degree; C has length DP+1.
pub fn mc01sy(
    dp: i32,
    alpha: f64,
    beta: f64,
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
) -> i32 {
    let dp = dp as usize;
    if dp + 1 > a.len() || dp + 1 > b.len() || dp + 1 > c.len() {
        return -5;
    }
    for i in 0..=dp {
        c[i] = alpha * a[i] + beta * b[i];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01sy() {
        let a = [1.0, 1.0];
        let b = [0.0, 1.0];
        let mut c = [0.0; 2];
        assert_eq!(mc01sy(1, 1.0, 1.0, &a, &b, &mut c), 0);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
    }
}
