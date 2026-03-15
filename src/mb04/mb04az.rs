//! MB04AZ — Matrix operation (SLICOT MB04). Stub.

use nalgebra::DMatrix;

/// Stub. Returns 0.
pub fn mb04az(_n: usize, _a: &DMatrix<f64>, _b: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb04az_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb04az(0, &a, &mut b), 0);
    }
}
