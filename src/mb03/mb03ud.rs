//! MB03UD — Matrix operation (SLICOT MB03). Stub.

use nalgebra::DMatrix;

/// Stub. Returns 0.
pub fn mb03ud(_n: usize, _a: &DMatrix<f64>, _b: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb03ud_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb03ud(0, &a, &mut b), 0);
    }
}
