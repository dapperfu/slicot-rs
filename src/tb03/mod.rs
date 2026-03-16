//! TB03 — Polynomial matrix representation from state-space (SLICOT TB03* routines)
//!
//! Left/right polynomial matrix representation (inv(P(s))*Q(s) or Q(s)*inv(P(s)))
//! with the same transfer matrix as a given state-space representation.

pub mod tb03ad;
pub mod tb03ay;
