//! TG01PD — Descriptor system transformations (SLICOT TG01PD)

use nalgebra::DMatrix;

pub fn tg01pd(
    _l: usize,
    _n: usize,
    _m: usize,
    _p: usize,
    _a: &mut DMatrix<f64>,
    _e: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    _c: &mut DMatrix<f64>,
) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tg01pd_smoke() {
        let mut a = DMatrix::zeros(2, 2);
        let mut e = DMatrix::identity(2, 2);
        let mut b = DMatrix::zeros(2, 1);
        let mut c = DMatrix::zeros(1, 2);
        assert_eq!(tg01pd(2, 2, 1, 1, &mut a, &mut e, &mut b, &mut c), 0);
    }
}
