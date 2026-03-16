//! # MA02ED — Complete a symmetric matrix from one triangle
//!
//! ## What the SLICOT routine does
//!
//! **MA02ED** takes a matrix that stores only the **upper** or **lower** triangle of a symmetric
//! matrix and fills in the other triangle so that `A(i,j) = A(j,i)`. This is essential when
//! downstream routines (e.g. eigensolvers, Cholesky) expect a full symmetric matrix but your
//! data is stored in packed form to save memory or to match the output of routines that only
//! write one triangle.
//!
//! ## Why it exists
//!
//! In control and linear algebra libraries, symmetric matrices (e.g. Gramians, covariance
//! matrices, Hessians) are often computed or stored in triangular form. MA02ED provides a
//! standard way to "expand" that representation to full storage without changing the
//! mathematical object. It is used in SLICOT in model reduction (Gramians), Lyapunov/Sylvester
//! solvers, and anywhere symmetric matrices are passed to LAPACK-style routines that expect
//! full matrices.
//!
//! ## When to use it
//!
//! - You have the upper (or lower) triangle of a symmetric matrix and need the full matrix.
//! - You are interfacing with code that expects `A(i,j)` and `A(j,i)` both defined.
//!
//! This example builds a small symmetric matrix from the upper triangle and prints the result.

use nalgebra::DMatrix;
use slicot_rs::ma02::ma02ed::{ma02ed, Ma02EdUplo};

fn main() {
    // Upper triangle: [1, 2; 0, 3] (row-major). Lower (1,0) is 0; MA02ED fills it with 2.
    let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
    let info = ma02ed(Ma02EdUplo::Upper, &mut a);
    assert_eq!(info, 0);
    assert!((a[(1, 0)] - 2.0_f64).abs() < 1e-15);
    println!("Symmetric matrix (upper → lower filled):\n{a}");
}
