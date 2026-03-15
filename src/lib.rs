//! # slicot-rs
//!
//! Pure Rust 1:1 mapping of SLICOT (Subroutine Library in Control Theory) routines.
//! Uses [nalgebra](https://crates.io/crates/nalgebra) for linear algebra; no FFI to LAPACK/BLAS.

pub mod ab01;
pub mod ab04;
pub mod ab05;
pub mod ab07;
pub mod ab08;
pub mod ab09;
pub mod ab13;
pub mod ab8n;
pub mod tb01;
