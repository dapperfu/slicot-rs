//! sb02mt — Stub (SLICOT sb02mt)
pub fn sb02mt(_n: usize, _a: &mut [f64], _lda: usize) -> i32 { -1 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert!(sb02mt(1, &mut [], 1) != 0); } }
