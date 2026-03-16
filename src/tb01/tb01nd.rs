//! TB01ND — Reduce (A,C) to observer Hessenberg form (SLICOT TB01ND)
//!
//! Dual of controller Hessenberg: reduce pair (A,C) to upper or lower observer Hessenberg using unitary U.

use nalgebra::DMatrix;
use crate::tb01::tb01md::{JobU, Uplo};

/// Reduces (A,C) to observer Hessenberg form. Observer form is the dual of controller form:
/// apply the same algorithm to (A', C') and then transpose back, i.e. reduce (A', C') to controller Hessenberg
/// => (U'*A'*U, C'*U) so (U'*A*U, C*U) in terms of (A,C) is (U'*A*U, C*U). So we can call controller Hessenberg
/// on (A', C') with B replaced by C', then get A' -> U'*A'*U, "B" -> U'*C'. So A -> U'*A*U, and we need C*U.
/// Controller Hessenberg on (A', C') gives: A' <- U'*A'*U, C' <- U'*C' (as the "B" part). So A = (U'*A'*U)' = U'*A*U, and C' becomes U'*C' so C = (U'*C')' = C*U. So we transpose A and C, call TB01MD-style on (A', C') treating C' as B, then transpose A back.
pub fn tb01nd(
    jobu: JobU,
    uplo: Uplo,
    a: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    u: &mut Option<&mut DMatrix<f64>>,
) -> i32 {
    let n = a.nrows();
    let p = c.nrows();
    if a.ncols() != n || c.ncols() != n {
        return -5;
    }
    if (jobu == JobU::Init || jobu == JobU::Update) && u.is_none() {
        return -9;
    }
    if let Some(ref uu) = *u {
        if uu.nrows() != n || uu.ncols() != n {
            return -9;
        }
    }
    if n == 0 || p == 0 {
        return 0;
    }
    if uplo != Uplo::Upper {
        return -2;
    }
    // Dual: (A,C) observer Hessenberg <=> (A', C') controller Hessenberg. So form At = A', Bt = C' (n×p).
    let mut at = a.transpose();
    let mut bt = c.transpose();
    if jobu == JobU::Init {
        if let Some(ref mut uu) = *u {
            uu.fill(0.0);
            for i in 0..n {
                uu[(i, i)] = 1.0;
            }
        }
    }
    let info = crate::tb01::tb01md::tb01md(jobu, uplo, &mut at, &mut bt, u);
    if info != 0 {
        return info;
    }
    *a = at.transpose();
    *c = bt.transpose();
    0
}
