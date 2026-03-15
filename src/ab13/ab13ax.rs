//! AB13AX — Hankel norm of a stable system with A in real Schur form (SLICOT).
//!
//! Hankel norm = max Hankel singular value = sqrt(max eigenvalue of W*X)
//! where X = observability Gramian (A'*X + X*A = -C'*C), W = controllability (A*W + W*A' = -B*B').

use nalgebra::DMatrix;

use crate::ab13::lyapunov;

/// Power iteration to get the largest eigenvalue (in magnitude) of a matrix.
fn power_iteration_max_eig(m: &DMatrix<f64>, max_iter: usize) -> f64 {
    let n = m.nrows();
    if n == 0 {
        return 0.0;
    }
    let mut v = DMatrix::<f64>::from_fn(n, 1, |_, _| 1.0 / (n as f64).sqrt());
    for _ in 0..max_iter {
        let u = m * &v;
        let nrm = u.norm();
        if nrm < 1e-15 {
            break;
        }
        v = u / nrm;
    }
    let rq = (v.transpose() * m * &v)[(0, 0)];
    rq.abs()
}

/// Continuous ('C') or discrete ('D') time. Discrete uses A*W*A' - W = -B*B' and A'*X*A - X = -C'*C.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// Computes the Hankel norm of (A,B,C). A must be in real Schur form and stable.
/// Returns 0 on success; hnorm is set to the Hankel norm. info=1 if A not stable, 2 if Lyapunov failed.
pub fn ab13ax(
    dico: Dico,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    hnorm: &mut f64,
    dwork: &mut [f64],
) -> i32 {
    *hnorm = 0.0;
    let n = a.nrows();
    if n == 0 {
        return 0;
    }
    if a.ncols() != n || b.nrows() != n || c.ncols() != n {
        return -1;
    }
    let at = a.transpose();

    let mut x = DMatrix::<f64>::zeros(n, n);
    let mut w = DMatrix::<f64>::zeros(n, n);
    let q_obs = &c.transpose() * c;
    let q_ctr = b * b.transpose();

    let ok_x = match dico {
        Dico::Continuous => lyapunov::lyapunov_continuous(&at, &q_obs, &mut x),
        Dico::Discrete => {
            // A'*X*A - X = -C'*C => solve discrete Lyapunov (we use vec form: (A'⊗A' - I)*vec(X) = vec(-C'*C))
            let n2 = n * n;
            if dwork.len() < n2 * n2 + n2 {
                return -2;
            }
            let mut k = DMatrix::<f64>::zeros(n2, n2);
            for i in 0..n {
                for j in 0..n {
                    for p in 0..n {
                        for q in 0..n {
                            let row = i * n + j;
                            let col = p * n + q;
                            k[(row, col)] = at[(i, p)] * at[(j, q)] - if (i, j) == (p, q) { 1.0 } else { 0.0 };
                        }
                    }
                }
            }
            let mut rhs = DMatrix::<f64>::zeros(n2, 1);
            for i in 0..n {
                for j in 0..n {
                    rhs[(i * n + j, 0)] = -q_obs[(i, j)];
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
    };
    if !ok_x {
        return 2;
    }

    let ok_w = match dico {
        Dico::Continuous => {
            lyapunov::lyapunov_continuous(a, &q_ctr, &mut w)
        }
        Dico::Discrete => {
            let n2 = n * n;
            let mut k = DMatrix::<f64>::zeros(n2, n2);
            for i in 0..n {
                for j in 0..n {
                    for p in 0..n {
                        for q in 0..n {
                            let row = i * n + j;
                            let col = p * n + q;
                            k[(row, col)] = a[(i, p)] * a[(j, q)] - if (i, j) == (p, q) { 1.0 } else { 0.0 };
                        }
                    }
                }
            }
            let mut rhs = DMatrix::<f64>::zeros(n2, 1);
            for i in 0..n {
                for j in 0..n {
                    rhs[(i * n + j, 0)] = -q_ctr[(i, j)];
                }
            }
            let lu = nalgebra::linalg::LU::new(k);
            if let Some(vec_w) = lu.solve(&rhs) {
                for i in 0..n {
                    for j in 0..n {
                        w[(i, j)] = vec_w[(i * n + j, 0)];
                    }
                }
                true
            } else {
                false
            }
        }
    };
    if !ok_w {
        return 2;
    }

    let wx = &w * &x;
    let max_eig = power_iteration_max_eig(&wx, 50);
    *hnorm = max_eig.max(0.0).sqrt();
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13ax_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let b = DMatrix::<f64>::zeros(0, 0);
        let c = DMatrix::<f64>::zeros(0, 0);
        let mut hnorm = -1.0;
        let mut dwork = vec![0.0; 1];
        assert_eq!(ab13ax(Dico::Continuous, &a, &b, &c, &mut hnorm, &mut dwork), 0);
        assert_eq!(hnorm, 0.0);
    }

    #[test]
    fn test_ab13ax_1x1_stable() {
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut hnorm = 0.0;
        let mut dwork = vec![0.0; 4];
        assert_eq!(ab13ax(Dico::Continuous, &a, &b, &c, &mut hnorm, &mut dwork), 0);
        assert!(hnorm > 0.0 && hnorm < 1.0);
    }
}
