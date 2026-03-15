//! Lyapunov equation solver: A'*X + X*A = -Q (continuous).
//! Uses Kronecker product: (A'⊗I + I⊗A')*vec(X) = vec(-Q).

use nalgebra::DMatrix;

/// Solves the continuous-time Lyapunov equation A'*X + X*A = -Q.
/// Returns true on success, false if the Kronecker matrix is singular.
/// Works for any n (memory O(n^4), time O(n^6)); suitable for small/medium n.
pub fn lyapunov_continuous(a: &DMatrix<f64>, q: &DMatrix<f64>, x: &mut DMatrix<f64>) -> bool {
    let n = a.nrows();
    if n == 0 {
        return true;
    }
    if a.ncols() != n || q.nrows() != n || q.ncols() != n || x.nrows() != n || x.ncols() != n {
        return false;
    }
    let at = a.transpose();
    // K = A'⊗I + I⊗A' is n^2 × n^2. vec(X) = n^2, vec(-Q) = n^2.
    let n2 = n * n;
    let mut k = DMatrix::<f64>::zeros(n2, n2);
    for i in 0..n {
        for j in 0..n {
            for p in 0..n {
                for qq in 0..n {
                    let row = i * n + j;
                    let col = p * n + qq;
                    let v = if j == qq { at[(i, p)] } else { 0.0 } + if i == p { at[(j, qq)] } else { 0.0 };
                    k[(row, col)] = v;
                }
            }
        }
    }
    let mut rhs = DMatrix::<f64>::zeros(n2, 1);
    for i in 0..n {
        for j in 0..n {
            rhs[(i * n + j, 0)] = -q[(i, j)];
        }
    }
    let lu = nalgebra::linalg::LU::new(k.clone());
    let vec_x = match lu.solve(&rhs) {
        Some(v) => v,
        None => return false,
    };
    for i in 0..n {
        for j in 0..n {
            x[(i, j)] = vec_x[(i * n + j, 0)];
        }
    }
    true
}

/// Solves the discrete-time Lyapunov (Stein) equation A*X*A' - X = -Q.
/// Returns true on success, false if the Kronecker matrix is singular.
pub fn lyapunov_discrete_stein(a: &DMatrix<f64>, q: &DMatrix<f64>, x: &mut DMatrix<f64>) -> bool {
    let n = a.nrows();
    if n == 0 {
        return true;
    }
    if a.ncols() != n || q.nrows() != n || q.ncols() != n || x.nrows() != n || x.ncols() != n {
        return false;
    }
    let n2 = n * n;
    let mut k = DMatrix::<f64>::zeros(n2, n2);
    for i in 0..n {
        for j in 0..n {
            for p in 0..n {
                for qq in 0..n {
                    let row = i * n + j;
                    let col = p * n + qq;
                    k[(row, col)] = a[(i, p)] * a[(j, qq)] - if (i, j) == (p, qq) { 1.0 } else { 0.0 };
                }
            }
        }
    }
    let mut rhs = DMatrix::<f64>::zeros(n2, 1);
    for i in 0..n {
        for j in 0..n {
            rhs[(i * n + j, 0)] = -q[(i, j)];
        }
    }
    let lu = nalgebra::linalg::LU::new(k);
    if let Some(vec_x) = lu.solve(&rhs) {
        for i in 0..n {
            for j in 0..n {
                x[(i, j)] = vec_x[(i * n + j, 0)];
            }
        }
        true
    } else {
        false
    }
}

/// Solves A*X + X*A' = -Q (continuous-time controllability Lyapunov).
/// Returns true on success.
pub fn lyapunov_continuous_dual(a: &DMatrix<f64>, q: &DMatrix<f64>, x: &mut DMatrix<f64>) -> bool {
    let at = a.transpose();
    lyapunov_continuous(&at, q, x)
}

/// Solves A'*X*A - X = -Q (discrete-time observability Stein).
/// Returns true on success.
pub fn lyapunov_discrete_stein_dual(a: &DMatrix<f64>, q: &DMatrix<f64>, x: &mut DMatrix<f64>) -> bool {
    let n = a.nrows();
    if n == 0 {
        return true;
    }
    if a.ncols() != n || q.nrows() != n || q.ncols() != n || x.nrows() != n || x.ncols() != n {
        return false;
    }
    let at = a.transpose();
    let n2 = n * n;
    let mut k = DMatrix::<f64>::zeros(n2, n2);
    for i in 0..n {
        for j in 0..n {
            for p in 0..n {
                for qq in 0..n {
                    let row = i * n + j;
                    let col = p * n + qq;
                    k[(row, col)] = at[(i, p)] * at[(j, qq)] - if (i, j) == (p, qq) { 1.0 } else { 0.0 };
                }
            }
        }
    }
    let mut rhs = DMatrix::<f64>::zeros(n2, 1);
    for i in 0..n {
        for j in 0..n {
            rhs[(i * n + j, 0)] = -q[(i, j)];
        }
    }
    let lu = nalgebra::linalg::LU::new(k);
    if let Some(vec_x) = lu.solve(&rhs) {
        for i in 0..n {
            for j in 0..n {
                x[(i, j)] = vec_x[(i * n + j, 0)];
            }
        }
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lyapunov_1x1() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let q = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut x = DMatrix::<f64>::zeros(1, 1);
        assert!(lyapunov_continuous(&a, &q, &mut x));
        assert!((x[(0, 0)] - 0.5).abs() < 1e-10);
    }
}
