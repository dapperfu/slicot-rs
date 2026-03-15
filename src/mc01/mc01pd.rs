//! MC01PD — Real polynomial from zeros, increasing powers (SLICOT MC01PD)
//
// Same as MC01PY but output in increasing powers of x.

/// Coefficients in increasing powers: P[0]=constant, P[K]=leading.
pub fn mc01pd(k: i32, rez: &[f64], imz: &[f64], p: &mut [f64], dwork: &mut [f64]) -> i32 {
    let k = k as usize;
    if k + 1 > p.len() {
        return -4;
    }
    let info = super::mc01py::mc01py(k as i32, rez, imz, p, dwork);
    if info != 0 {
        return info;
    }
    p[0..=k].reverse();
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01pd() {
        let mut p = [0.0; 2];
        let mut dwork = [0.0; 1];
        assert_eq!(mc01pd(1, &[1.0], &[0.0], &mut p, &mut dwork), 0);
        assert!((p[0] - (-1.0)).abs() < 1e-10);
        assert!((p[1] - 1.0).abs() < 1e-10);
    }
}
