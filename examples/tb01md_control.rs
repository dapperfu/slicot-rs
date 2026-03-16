//! Example: TB01MD — reduce (B, A) to controller Hessenberg form.
//!
//! Small LTI system: state dimension 3, inputs 2. Transforms (B, A) to upper controller Hessenberg.

use nalgebra::DMatrix;
use slicot_rs::tb01::tb01md::{tb01md, JobU, Uplo};

fn main() {
    let mut a = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
    let mut b = DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
    let mut u = None;
    let info = tb01md(JobU::No, Uplo::Upper, &mut a, &mut b, &mut u);
    assert_eq!(info, 0);
    println!("A (controller Hessenberg form):\n{a}");
    println!("B (transformed):\n{b}");
}
