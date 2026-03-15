//! SB08ND — SLICOT SB08ND. Stub.
use nalgebra::DMatrix;
pub fn sb08nd(_n: usize, _a: &DMatrix<f64>, _x: &mut DMatrix<f64>) -> i32 { 0 }
#[cfg(test)]
mod tests { use super::*; #[test] fn test_sb08nd() { let a = DMatrix::zeros(1,1); let mut x = DMatrix::zeros(1,1); assert_eq!(sb08nd(1,&a,&mut x),0); } }
