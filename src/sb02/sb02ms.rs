//! sb02ms — Stub (SLICOT sb02ms)
pub fn sb02ms(_n: usize, _a: &mut [f64], _lda: usize) -> i32 { -1 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert!(sb02ms(1, &mut [], 1) != 0); } }
