//! MB05OD — Matrix operation (SLICOT MB05). Stub.

use nalgebra::DMatrix;

/// Stub. Returns 0.
pub fn mb05od(_n: usize, _a: &DMatrix<f64>, _b: &mut DMatrix<f64>) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb05od_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb05od(0, &a, &mut b), 0);
    }
}
