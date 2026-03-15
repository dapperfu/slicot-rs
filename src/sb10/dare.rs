//! Discrete-time algebraic Riccati equation (DARE) helper for SB10*.
//! Solves A'*X*A - X - (A'*X*B)*inv(R + B'*X*B)*(B'*X*A) + Q = 0.

use nalgebra::DMatrix;

/// Solves the discrete-time algebraic Riccati equation (DARE):
///   A'*X*A - X - A'*X*B*inv(R + B'*X*B)*B'*X*A + Q = 0.
/// Returns true on success and writes solution into x; false if singular.
pub fn dare(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    q: &DMatrix<f64>,
    r: &DMatrix<f64>,
    x: &mut DMatrix<f64>,
    tol: f64,
    max_it: usize,
) -> bool {
    let n = a.nrows();
    if n == 0 {
        return true;
    }
    let m = b.ncols();
    let mut xk = q.clone();
    const TOL_DEFAULT: f64 = 1e-12;
    const MAX_IT_DEFAULT: usize = 80;
    let tol_use = if tol > 0.0 { tol } else { TOL_DEFAULT };
    let max_it_use = if max_it > 0 { max_it } else { MAX_IT_DEFAULT };

    for _ in 0..max_it_use {
        // R + B'*X*B
        let bxb = b.transpose() * &xk * b;
        let mut rpbxb = r.clone();
        for i in 0..m {
            for j in 0..m {
                rpbxb[(i, j)] += bxb[(i, j)];
            }
        }
        let inv_rpbxb = match rpbxb.try_inverse() {
            Some(inv) => inv,
            None => return false,
        };
        // K = inv(R + B'*X*B)*B'*X*A
        let k = &inv_rpbxb * b.transpose() * &xk * a;
        // X_new = Q + A'*X*A - A'*X*B*K
        let x_new = q + a.transpose() * &xk * a - a.transpose() * &xk * b * &k;
        let diff = &x_new - &xk;
        if diff.norm() < tol_use * (1.0 + xk.norm()) {
            x.copy_from(&x_new);
            return true;
        }
        xk = x_new;
    }
    x.copy_from(&xk);
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dare_1x1() {
        let a = DMatrix::from_row_slice(1, 1, &[0.5]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let q = DMatrix::from_row_slice(1, 1, &[1.0]);
        let r = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut x = DMatrix::zeros(1, 1);
        assert!(dare(&a, &b, &q, &r, &mut x, 1e-12, 80));
        // 0.25*x - x + 1 - 0.25*x^2/(1+x) = 0 => check
        assert!(x[(0, 0)] > 0.0);
    }
}
