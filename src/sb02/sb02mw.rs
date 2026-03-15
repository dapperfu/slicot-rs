//! sb02mw — Stub (SLICOT sb02mw)
pub fn sb02mw(_n: usize, _a: &mut [f64], _lda: usize) -> i32 { -1 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert!(sb02mw(1, &mut [], 1) != 0); } }
