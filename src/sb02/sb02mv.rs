//! sb02mv — Stub (SLICOT sb02mv)
pub fn sb02mv(_n: usize, _a: &mut [f64], _lda: usize) -> i32 { -1 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert!(sb02mv(1, &mut [], 1) != 0); } }
