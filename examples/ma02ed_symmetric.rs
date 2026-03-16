//! Example: MA02ED — complete a symmetric matrix from one triangle.
//!
//! Build a 2×2 symmetric matrix from the upper triangle only; MA02ED fills the lower triangle.

use nalgebra::DMatrix;
use slicot_rs::ma02::ma02ed::{ma02ed, Ma02EdUplo};

fn main() {
    // Upper triangle: [1, 2; 0, 3] (row-major). Lower (1,0) is 0, will be filled with 2.
    let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
    let info = ma02ed(Ma02EdUplo::Upper, &mut a);
    assert_eq!(info, 0);
    assert!((a[(1, 0)] - 2.0_f64).abs() < 1e-15);
    println!("Symmetric matrix (upper → lower filled):\n{a}");
}
