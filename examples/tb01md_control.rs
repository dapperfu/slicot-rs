//! # TB01MD — Reduce (A, B) to controller Hessenberg form
//!
//! ## What the SLICOT routine does
//!
//! **TB01MD** reduces the state-space pair \((A, B)\) to **upper controller Hessenberg form** by
//! orthogonal similarity. That is, it finds an orthogonal matrix \(U\) such that
//! \((\tilde{A}, \tilde{B}) = (U^\top A U, U^\top B)\) has \(B\) with a staircase structure and
//! \(\tilde{A}\) in upper Hessenberg form. This structure is the first step in many
//! controllability and stabilizability algorithms (e.g. staircase forms, pole placement).
//!
//! ## Why it exists
//!
//! Controller Hessenberg form is a standard condensed form in control theory. It reveals the
//! controllable subspace and reduces the problem size for subsequent steps (e.g. AB01ND-style
//! staircase, or feedback design). TB01MD is the multi-input counterpart of the single-input
//! controller Hessenberg reduction and is used throughout SLICOT for analysis and synthesis.
//!
//! ## When to use it
//!
//! - Before computing controllability indices or staircase forms.
//! - As a preprocessing step for pole assignment or LQR-type design.
//! - When you need a condensed \((A,B)\) pair for numerical stability.
//!
//! This example reduces a small LTI system (3 states, 2 inputs) to controller Hessenberg form.

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
