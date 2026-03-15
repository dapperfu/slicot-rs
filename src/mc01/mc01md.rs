//! MC01MD — Leading coefficients of shifted polynomial (SLICOT MC01MD)
//
// For real polynomial P(x) and scalar alpha, computes leading K coefficients
// of P(x) = q(1) + q(2)*(x-alpha) + ... + q(K)*(x-alpha)^(K-1) + ... using Horner.

/// Leading K coefficients of shifted polynomial. P and Q are coefficients in increasing powers.
/// Q must have length >= DP+1; first K elements filled, rest used as workspace.
pub fn mc01md(dp: i32, alpha: f64, k: i32, p: &[f64], q: &mut [f64]) -> i32 {
    let dp = dp as usize;
    let k = k as usize;
    if dp.saturating_add(1) > p.len() {
        return -4;
    }
    if dp + 1 > q.len() {
        return -5;
    }
    if k < 1 || k > dp + 1 {
        return -3;
    }
    // Horner at alpha gives q(1)=P(alpha). Then synthetic division by (x-alpha) gives next poly.
    let mut work: Vec<f64> = p[..=dp].to_vec();
    for i in 0..k {
        let deg = dp - i;
        // Evaluate current polynomial at alpha -> q[i]
        let mut val = work[deg];
        for j in (0..deg).rev() {
            val = val * alpha + work[j];
        }
        q[i] = val;
        if i + 1 < k && deg > 0 {
            // Quotient D(x) = (current - val)/(x-alpha): D[deg-1]=work[deg], D[j]=work[j+1]+alpha*D[j+1]
            let mut saved = work[deg - 1];
            work[deg - 1] = work[deg];
            for j in (0..deg - 1).rev() {
                let newd = saved + alpha * work[j + 1];
                saved = work[j];
                work[j] = newd;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01md_example() {
        // DP=5, ALPHA=2, K=6, P = 6,5,4,3,2,1 (increasing powers)
        let p = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let mut q = [0.0; 6];
        assert_eq!(mc01md(5, 2.0, 6, &p, &mut q), 0);
        assert!((q[0] - 120.0).abs() < 1e-6);
        assert!((q[1] - 201.0).abs() < 1e-6);
        assert!((q[2] - 150.0).abs() < 1e-6);
        assert!((q[3] - 59.0).abs() < 1e-6);
        assert!((q[4] - 12.0).abs() < 1e-6);
        assert!((q[5] - 1.0).abs() < 1e-6);
    }
}
