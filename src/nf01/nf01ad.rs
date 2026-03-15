//! nf01ad — Stub (SLICOT nf01ad)
pub fn nf01ad(_n: i32, _x: &[f64], _dwork: &mut [f64]) -> i32 { 0 }
#[cfg(test)] mod tests { use super::*; #[test] fn stub() { assert_eq!(nf01ad(1, &[], &mut []), 0); } }
