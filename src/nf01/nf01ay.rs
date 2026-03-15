//! nf01ay — Stub (SLICOT nf01ay)
pub fn nf01ay(_n: i32, _x: &[f64], _dwork: &mut [f64]) -> i32 { 0 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert_eq!(nf01ay(1, &[], &mut []), 0); } }
