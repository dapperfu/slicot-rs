//! MB04OD — Matrix operation (SLICOT MB04). Stub.

use nalgebra::DMatrix;

/// Stub. Returns 0.
pub fn mb04od(_n: usize, _a: &DMatrix<f64>, _b: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb04od_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb04od(0, &a, &mut b), 0);
    }
}
