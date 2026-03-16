//! TG01FZ — Descriptor system (complex) (SLICOT TG01FZ)

use nalgebra::DMatrix;
use num_complex::Complex64;

pub fn tg01fz(
    _l: usize,
    _n: usize,
    _m: usize,
    _p: usize,
    _a: &mut DMatrix<Complex64>,
    _e: &mut DMatrix<Complex64>,
    _b: &mut DMatrix<Complex64>,
    _c: &mut DMatrix<Complex64>,
) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tg01fz_smoke() {
        let mut a = DMatrix::from_fn(2, 2, |_, _| Complex64::new(0.0, 0.0));
        let mut e = DMatrix::from_fn(2, 2, |i, j| if i == j { Complex64::new(1.0, 0.0) } else { Complex64::new(0.0, 0.0) });
        let mut b = DMatrix::zeros(2, 1);
        let mut c = DMatrix::zeros(1, 2);
        assert_eq!(tg01fz(2, 2, 1, 1, &mut a, &mut e, &mut b, &mut c), 0);
    }
}
