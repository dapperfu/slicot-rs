//! MC01ND — Value of real polynomial at complex point (SLICOT MC01ND)
//
// Computes P(x0) for real polynomial P and complex x0 using Horner.

/// P has coefficients in increasing powers; returns P(x0) in (vr, vi).
pub fn mc01nd(dp: i32, xr: f64, xi: f64, p: &[f64], vr: &mut f64, vi: &mut f64) -> i32 {
    let dp = dp as usize;
    if dp.saturating_add(1) > p.len() {
        return -4;
    }
    let mut qr = p[dp];
    let mut qi = 0.0_f64;
    for j in (0..dp).rev() {
        // q = x0*q + p[j]
        let (new_qr, new_qi) = (xr * qr - xi * qi + p[j], xr * qi + xi * qr);
        qr = new_qr;
        qi = new_qi;
    }
    *vr = qr;
    *vi = qi;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01nd_example() {
        // DP=4, x0=-1.56+0.29j, P = 5,3,-1,2,1
        let p = [5.0, 3.0, -1.0, 2.0, 1.0];
        let mut vr = 0.0;
        let mut vi = 0.0;
        assert_eq!(mc01nd(4, -1.56, 0.29, &p, &mut vr, &mut vi), 0);
        assert!((vr - (-4.1337)).abs() < 1e-3);
        assert!((vi - 1.7088).abs() < 1e-3);
    }
}
