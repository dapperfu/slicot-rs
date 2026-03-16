//! MA02JD — Residual || Q^T Q - I ||_F for orthogonal symplectic Q (SLICOT MA02JD)
//
// Q = [ op(Q1) op(Q2); -op(Q2) op(Q1) ]. Computes Frobenius norm of Q^T*Q - I.

use nalgebra::DMatrix;

/// Computes the residual || Q^T Q - I ||_F. Q is 2n×2n built from Q1, Q2 as above.
/// res is workspace n×n (used internally).
pub fn ma02jd(
    tran1: bool,
    tran2: bool,
    q1: &DMatrix<f64>,
    q2: &DMatrix<f64>,
    res: &mut DMatrix<f64>,
) -> f64 {
    let n = q1.nrows();
    if n == 0 || q1.ncols() != n || q2.nrows() != n || q2.ncols() != n || res.nrows() != n || res.ncols() != n {
        return 0.0;
    }
    let op1 = if tran1 { q1.transpose() } else { q1.clone() };
    let op2 = if tran2 { q2.transpose() } else { q2.clone() };
    // Q = [ op1 op2; -op2 op1 ], 2n×2n
    // Q^T*Q - I: block (1,1) = op1'*op1 + op2'*op2 - I
    // block (1,2) = op1'*op2 - op2'*op1
    // block (2,1) = -op2'*op1 + op1'*op2 = -(op2'*op1 - op1'*op2)
    // block (2,2) = op2'*op2 + op1'*op1 - I
    let p11 = &op1.transpose() * &op1 + &op2.transpose() * &op2;
    let p12 = &op1.transpose() * &op2 - &op2.transpose() * &op1;
    let mut sum_sq = 0.0;
    for i in 0..n {
        for j in 0..n {
            let v11 = p11[(i, j)] - if i == j { 1.0 } else { 0.0 };
            sum_sq += v11 * v11;
        }
    }
    for i in 0..n {
        for j in 0..n {
            sum_sq += p12[(i, j)] * p12[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..n {
            let v21 = -p12[(j, i)];
            sum_sq += v21 * v21;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let v22 = p11[(i, j)] - if i == j { 1.0 } else { 0.0 };
            sum_sq += v22 * v22;
        }
    }
    sum_sq.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02jd_identity() {
        let q1 = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let q2 = DMatrix::from_row_slice(2, 2, &[0.0, 0.0, 0.0, 0.0]);
        let mut res = DMatrix::zeros(2, 2);
        let r = ma02jd(false, false, &q1, &q2, &mut res);
        assert!(r < 1e-10, "identity Q should give residual ~0, got {}", r);
    }

    #[test]
    fn test_ma02jd_zero_dim() {
        let q1 = DMatrix::<f64>::zeros(0, 0);
        let q2 = DMatrix::<f64>::zeros(0, 0);
        let mut res = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(ma02jd(false, false, &q1, &q2, &mut res), 0.0);
    }
}
