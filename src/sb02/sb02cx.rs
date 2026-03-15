//! sb02cx — Stub (SLICOT sb02cx)
pub fn sb02cx(_n: usize, _a: &mut [f64], _lda: usize) -> i32 { -1 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert!(sb02cx(1, &mut [], 1) != 0); } }
