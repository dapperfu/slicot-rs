//! # AB01ND — Orthogonal controllability staircase form (multi-input)
//!
//! ## What the SLICOT routine does
//!
//! **AB01ND** reduces the state-space pair \((A, B)\) to **orthogonal controllability staircase
//! form**. It computes an orthogonal matrix \(Z\) such that \((Z^\top A Z, Z^\top B)\) is in
//! block upper Hessenberg (staircase) form, and returns the **dimension of the controllable
//! subspace** `ncont`. The uncontrollable part of the system is isolated in the trailing
//! block of the transformed \(A\); the first `ncont` states correspond to the controllable
//! subsystem.
//!
//! ## Why it exists
//!
//! In control theory, knowing whether a system is controllable (and how many states are
//! controllable) is fundamental for stabilizability, pole placement, and model reduction.
//! The staircase form is numerically stable (orthogonal transformations) and reveals the
//! controllability indices. AB01ND is the multi-input version and is used by many SLICOT
//! routines that need a minimal or controllable realization.
//!
//! ## When to use it
//!
//! - You need the controllable subspace dimension or a controllable/uncontrollable split.
//! - You are preparing for minimal realization, model reduction, or feedback design.
//! - You want to check controllability without computing the controllability matrix (which
//!   is ill-conditioned for large \(n\)).
//!
//! This example runs AB01ND on a small (4×4, 2 inputs) pair and prints `ncont` and the
//! transformed matrices.

use nalgebra::DMatrix;
use slicot_rs::ab01::ab01nd::{ab01nd, JobZ};

fn main() {
    let n = 4_usize;
    let m = 2_usize;
    let mut a = DMatrix::from_row_slice(n, n, &[
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 0.0,
    ]);
    let mut b = DMatrix::from_row_slice(n, m, &[
        0.0, 0.0,
        1.0, 0.0,
        0.0, 0.0,
        0.0, 1.0,
    ]);
    let mut ncont = 0_usize;
    let mut indcon = 0_usize;
    let mut nblk = vec![0i32; n];
    let info = ab01nd(
        JobZ::No,
        n,
        m,
        &mut a,
        &mut b,
        &mut ncont,
        &mut indcon,
        &mut nblk,
        None,
        1e-10_f64,
    );
    if info == 0 {
        println!("Controllable subspace dimension: ncont = {}", ncont);
        println!("A (staircase form):\n{}", a);
        println!("B (transformed):\n{}", b);
    } else {
        println!("AB01ND returned INFO = {} (0 = success, 1 = not yet implemented)", info);
        println!("When implemented, ncont and the transformed (A,B) will be computed.");
    }
}
