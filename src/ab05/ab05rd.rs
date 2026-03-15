//! AB05RD — Similarity transformation of state-space model (SLICOT AB05RD). Not yet implemented.

use nalgebra::DMatrix;

/// Stub: validates dimensions, returns 1 (not implemented) or < 0 (invalid argument).
pub fn ab05rd(
    _n: usize,
    _m: usize,
    _p: usize,
    _a: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    _c: &mut DMatrix<f64>,
    _d: &mut DMatrix<f64>,
    _t: &DMatrix<f64>,
) -> i32 {
    if _n == 0 {
        return 0;
    }
    if _t.nrows() != _n || _t.ncols() != _n {
        return -9;
    }
    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab05rd_trivial() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut c = DMatrix::zeros(0, 0);
        let mut d = DMatrix::zeros(0, 0);
        let t = DMatrix::zeros(0, 0);
        assert_eq!(ab05rd(0, 0, 0, &mut a, &mut b, &mut c, &mut d, &t), 0);
    }
}
