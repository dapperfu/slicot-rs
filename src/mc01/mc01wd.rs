//! MC01WD — Complex polynomial from zeros (SLICOT MC01WD)
//
// Same as MC01OD: P(x) = prod(x - z_i). Output in increasing powers.

/// Delegates to MC01OD.
pub fn mc01wd(
    k: i32,
    rez: &[f64],
    imz: &[f64],
    rep: &mut [f64],
    imp: &mut [f64],
    dwork: &mut [f64],
) -> i32 {
    super::mc01od::mc01od(k, rez, imz, rep, imp, dwork)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01wd() {
        let mut rep = [0.0; 2];
        let mut imp = [0.0; 2];
        let mut dwork = [0.0; 4];
        assert_eq!(mc01wd(1, &[1.0], &[0.0], &mut rep, &mut imp, &mut dwork), 0);
        assert!((rep[0] - (-1.0)).abs() < 1e-10);
    }
}
