//! sb02mu — Stub (SLICOT sb02mu)
pub fn sb02mu(_n: usize, _a: &mut [f64], _lda: usize) -> i32 { -1 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert!(sb02mu(1, &mut [], 1) != 0); } }
