//! sb02mr — Stub (SLICOT sb02mr)
pub fn sb02mr(_n: usize, _a: &mut [f64], _lda: usize) -> i32 { -1 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert!(sb02mr(1, &mut [], 1) != 0); } }
