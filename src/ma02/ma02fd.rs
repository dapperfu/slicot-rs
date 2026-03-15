//! MA02FD — Modified hyperbolic plane rotation (SLICOT MA02FD)
//
// Computes c, s (c^2 + s^2 = 1) so that y1 = sqrt(x1^2 - x2^2), y2 = 0. Requires |x2| < |x1| or x1=x2=0.

/// Computes c, s for modified hyperbolic rotation; overwrites x1 with y1 = sqrt(x1^2 - x2^2).
/// Returns 0 on success; 1 if |x2| >= |x1| and (x1,x2) != (0,0).
pub fn ma02fd(x1: &mut f64, x2: f64, c: &mut f64, s: &mut f64) -> i32 {
    if (*x1 != 0.0 || x2 != 0.0) && x2.abs() >= x1.abs() {
        return 1;
    }
    *c = 1.0;
    *s = 0.0;
    if *x1 != 0.0 {
        *s = x2 / *x1;
        let t = (1.0 - *s) * (1.0 + *s);
        *c = if t >= 0.0 { t.sqrt() } else { 0.0 };
        if *x1 < 0.0 {
            *c = -(*c);
        }
        *x1 = *c * *x1;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02fd_zero() {
        let mut x1 = 0.0;
        let mut c = 0.0;
        let mut s = 0.0;
        assert_eq!(ma02fd(&mut x1, 0.0, &mut c, &mut s), 0);
        assert_eq!(c, 1.0);
        assert_eq!(s, 0.0);
    }
}
