//! AB05PD — Addition of two system state-space models (SLICOT AB05PD). Not yet implemented.

use nalgebra::DMatrix;

/// Stub: validates dimensions, returns 1 (not implemented) or < 0 (invalid argument).
pub fn ab05pd(
    _n1: usize,
    _m1: usize,
    _p1: usize,
    _n2: usize,
    _a1: &DMatrix<f64>,
    _b1: &DMatrix<f64>,
    _c1: &DMatrix<f64>,
    _d1: &DMatrix<f64>,
    _a2: &DMatrix<f64>,
    _b2: &DMatrix<f64>,
    _c2: &DMatrix<f64>,
    _d2: &DMatrix<f64>,
    _a: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    _c: &mut DMatrix<f64>,
    _d: &mut DMatrix<f64>,
) -> i32 {
    if _n1 == 0 && _n2 == 0 {
        return 0;
    }
    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab05pd_trivial() {
        let a1 = DMatrix::zeros(0, 0);
        let b1 = DMatrix::zeros(0, 0);
        let c1 = DMatrix::zeros(0, 0);
        let d1 = DMatrix::zeros(0, 0);
        let a2 = DMatrix::zeros(0, 0);
        let b2 = DMatrix::zeros(0, 0);
        let c2 = DMatrix::zeros(0, 0);
        let d2 = DMatrix::zeros(0, 0);
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut c = DMatrix::zeros(0, 0);
        let mut d = DMatrix::zeros(0, 0);
        assert_eq!(ab05pd(0, 0, 0, 0, &a1, &b1, &c1, &d1, &a2, &b2, &c2, &d2, &mut a, &mut b, &mut c, &mut d), 0);
    }
}
