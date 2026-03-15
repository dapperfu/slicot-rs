//! AB04MD — Bilinear transformation of state-space system (SLICOT AB04MD)
//!
//! Discrete <-> continuous via bilinear transform of (A,B,C,D).

use nalgebra::{DMatrix, linalg::LU};

/// Type of transformation: discrete-to-continuous or continuous-to-discrete.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum BilinearType {
    /// Discrete-time -> continuous-time.
    D2C,
    /// Continuous-time -> discrete-time.
    C2D,
}

/// Performs bilinear transformation on (A,B,C,D). Overwrites A, B, C, D in place.
///
/// # Arguments
/// * `typ` - D2C or C2D
/// * `alpha` - Parameter (must be nonzero)
/// * `beta` - Parameter (must be nonzero)
/// * `a` - N×N state matrix (overwritten)
/// * `b` - N×M input matrix (overwritten)
/// * `c` - P×N output matrix (overwritten)
/// * `d` - P×M feedthrough (overwritten)
///
/// # Returns
/// 0 on success; < 0 invalid argument; 1 = (alpha*I+A) singular (D2C); 2 = (beta*I-A) singular (C2D).
pub fn ab04md(
    typ: BilinearType,
    alpha: f64,
    beta: f64,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n || b.nrows() != n || c.ncols() != n || d.nrows() != p || d.ncols() != m {
        return -7;
    }
    if alpha == 0.0 {
        return -5;
    }
    if beta == 0.0 {
        return -6;
    }
    if n == 0 && m == 0 && p == 0 {
        return 0;
    }
    if n == 0 {
        return 0;
    }

    let (palpha, pbeta) = match typ {
        BilinearType::D2C => (alpha, beta),
        BilinearType::C2D => (-beta, -alpha),
    };
    let ab2 = palpha * pbeta * 2.0_f64;
    let sqrab2_abs = ab2.abs().sqrt();
    let sqrab2 = if palpha >= 0.0 {
        sqrab2_abs
    } else {
        -sqrab2_abs
    };

    let mut a_work = a.clone();
    for i in 0..n {
        a_work[(i, i)] += palpha;
    }
    let lu = LU::new(a_work.clone());
    if !lu.is_invertible() {
        return if typ == BilinearType::D2C { 1 } else { 2 };
    }
    let a_inv = match lu.try_inverse() {
        Some(inv) => inv,
        None => return if typ == BilinearType::D2C { 1 } else { 2 },
    };

    // B := (alpha*I+A)^{-1} * B
    let sol_b = &a_inv * b.clone();
    b.copy_from(&sol_b);

    // D := D - C * (alpha*I+A)^{-1} * B  = D - C * B (B already overwritten)
    *d -= &*c * b.clone();

    // B := sqrt(2*alpha*beta) * B
    *b *= sqrab2;

    // C := sqrt(2*alpha*beta) * C * (alpha*I+A)^{-1}
    let new_c = (&*c * &a_inv) * sqrab2;
    c.copy_from(&new_c);

    // A := beta*I - 2*alpha*beta*(alpha*I+A)^{-1}
    let mut new_a = &a_inv * (-ab2);
    for i in 0..n {
        new_a[(i, i)] += pbeta;
    }
    a.copy_from(&new_a);

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab04md_n0() {
        let mut a = DMatrix::from_row_slice(0, 0, &[]);
        let mut b = DMatrix::from_row_slice(0, 1, &[]);
        let mut c = DMatrix::from_row_slice(1, 0, &[]);
        let mut d = DMatrix::from_row_slice(1, 1, &[1.0]);
        let info = ab04md(BilinearType::D2C, 1.0, 1.0, &mut a, &mut b, &mut c, &mut d);
        assert_eq!(info, 0);
    }

    #[test]
    fn test_ab04md_singular() {
        let mut a = DMatrix::from_row_slice(1, 1, &[-1.0]); // alpha*I+A = 0 when alpha=1
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let info = ab04md(BilinearType::D2C, 1.0, 1.0, &mut a, &mut b, &mut c, &mut d);
        assert!(info == 1 || info == 2);
    }
}
