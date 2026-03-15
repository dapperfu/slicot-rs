//! MC01VD — Complex polynomial evaluation at complex point (SLICOT MC01VD)
//
// P(x0) for complex P and complex x0. Horner.

/// P stored as (PR, PI) in increasing powers. Result in (VR, VI).
pub fn mc01vd(
    dp: i32,
    xr: f64,
    xi: f64,
    pr: &[f64],
    pi: &[f64],
    vr: &mut f64,
    vi: &mut f64,
) -> i32 {
    let dp = dp as usize;
    if dp + 1 > pr.len() || dp + 1 > pi.len() {
        return -4;
    }
    let mut qr = pr[dp];
    let mut qi = pi[dp];
    for j in (0..dp).rev() {
        let (new_r, new_i) = (
            xr * qr - xi * qi + pr[j],
            xr * qi + xi * qr + pi[j],
        );
        qr = new_r;
        qi = new_i;
    }
    *vr = qr;
    *vi = qi;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01vd() {
        let pr = [1.0, 0.0, 1.0];
        let pi = [0.0; 3];
        let mut vr = 0.0;
        let mut vi = 0.0;
        assert_eq!(mc01vd(2, 1.0, 0.0, &pr, &pi, &mut vr, &mut vi), 0);
        assert!((vr - 2.0).abs() < 1e-10);
        assert!((vi - 0.0).abs() < 1e-10);
    }
}
