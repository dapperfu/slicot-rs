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
pub mod ag07;
pub mod ag08;
pub mod ag8b;
pub mod bb01;
pub mod bb02;
pub mod bb03;
pub mod bb04;
pub mod bd01;
pub mod bd02;
pub mod de01;
pub mod df01;
pub mod dg01;
pub mod dgeg;
pub mod dk01;
pub mod dlac;
pub mod dlat;
pub mod fb01;
pub mod fd01;
pub mod ib01;
pub mod ib03;
pub mod ma01;
pub mod ma02;
pub mod mb01;
pub mod mb02;
pub mod mb03;
pub mod mb04;
pub mod mb05;
pub mod mb3j;
pub mod mb3l;
pub mod mb3o;
pub mod mb3p;
pub mod mb4d;
pub mod mc01;
pub mod mc03;
pub mod md03;
pub mod nf01;
pub mod sb01;
pub mod sb02;
pub mod sb03;
pub mod sb04;
pub mod sb06;
pub mod sb08;
pub mod sb09;
pub mod sg02;
pub mod sg03;
// pub mod sb10;  // disabled: missing submodules
// pub mod sb16;  // disabled: compile errors in sb16bd/sb16cd
pub mod tb01;
pub mod tb03;
pub mod tb04;
pub mod tb05;
pub mod tc01;
pub mod tc04;
pub mod tc05;
pub mod td03;
pub mod td04;
pub mod td05;
pub mod tf01;
// pub mod tg01;  // disabled: missing submodules (tg01kd, tg01kz, ...)
pub mod ud01;
pub mod ue01;
pub mod zgeg;
pub mod zlat;

/// SLICOT .dat/.res I/O for fuzzer and tests (pilot: AB01ND).
pub mod slicot_io;
