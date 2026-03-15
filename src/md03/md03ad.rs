//! MD03AD — Stub (SLICOT MD03AD) — nonlinear least squares driver

/// Stub: returns 0.
pub fn md03ad(
    _m: i32,
    _n: i32,
    _x: &mut [f64],
    _ldx: usize,
    _y: &[f64],
    _tol: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
    _ldwork: i32,
) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_md03ad_stub() {
        let mut x = [0.0; 2];
        let mut iwork = [0; 10];
        let mut dwork = [0.0; 100];
        assert_eq!(md03ad(1, 1, &mut x, 1, &[], 1e-10, &mut iwork, &mut dwork, 100), 0);
    }
}
