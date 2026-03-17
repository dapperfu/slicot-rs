//! SB03OU — Solve Lyapunov op(A)'*X*op(A) - X = -scale^2*op(B)'*op(B) (or continuous).
//! QR or RQ of B, then SB03OT, then make diagonal of U non-negative.

use nalgebra::DMatrix;
use std::cmp::min;

use crate::sb03::sb03ot::sb03ot;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Solve for Cholesky U: op(A)'*U'*U + U'*U*op(A) = -scale^2*op(B)'*op(B) (continuous)
/// or op(A)'*U'*U*op(A) - U'*U = -scale^2*op(B)'*op(B) (discrete).
///
/// A is N×N real Schur (input); B is M×N (if !ltrans) or N×M (if ltrans), may be overwritten;
/// U is N×N output Cholesky factor (column-major LDU); TAU has length min(N,M); DWORK length >= 4*N.
/// Returns INFO: 0 = ok, 1 = singular warning, 2 = A not stable/convergent, 3 = block > 2×2, 4 = real 2×2 block, <0 = invalid argument.
pub fn sb03ou(
    discr: bool,
    ltrans: bool,
    n: usize,
    m: usize,
    a: &[f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    tau: &mut [f64],
    u: &mut [f64],
    ldu: usize,
    scale: &mut f64,
    dwork: &mut [f64],
) -> i32 {
    if lda < n.max(1) {
        return -6;
    }
    if (!ltrans && ldb < m.max(1)) || (ltrans && ldb < n.max(1)) {
        return -8;
    }
    if ldu < n.max(1) {
        return -11;
    }
    if dwork.len() < 4 * n {
        return -14;
    }
    if n == 0 {
        *scale = ONE;
        if !dwork.is_empty() {
            dwork[0] = ONE;
        }
        return 0;
    }
    if m == 0 {
        *scale = ONE;
        for j in 0..n {
            for i in 0..=j {
                u[i + j * ldu] = ZERO;
            }
        }
        if !dwork.is_empty() {
            dwork[0] = (4 * n) as f64;
        }
        return 0;
    }

    let mn = min(n, m);

    if ltrans {
        // B is N×M. RQ factorization: B = R*Q. We get R and build F in U.
        // RQ via QR of B': B' = Q*R => B = R'*Q'. First N rows of R (when M>=N) are N×N upper.
        let bt = DMatrix::from_fn(m, n, |i, j| b[j + i * ldb]);
        let qr = bt.qr();
        let r = qr.r();

        if m >= n {
            // R is M×N; first N rows are N×N upper. Copy to U.
            for j in 0..n {
                for i in 0..=j {
                    u[i + j * ldu] = r[(i, j)];
                }
                for i in (j + 1)..n {
                    u[i + j * ldu] = ZERO;
                }
            }
        } else {
            // R is M×N; M×M upper block at R(0..M, 0..M). Put in U(N-M..N, N-M..N); zero rest.
            for j in 0..n {
                for i in 0..n {
                    u[i + j * ldu] = ZERO;
                }
            }
            for j in 0..m {
                for i in 0..=j {
                    u[(n - m + i) + (n - m + j) * ldu] = r[(i, j)];
                }
            }
        }
        // Overwrite B with RQ factor details (tau etc.) for API compatibility - we don't need to preserve B
        // DGERQF overwrites B; we already built U from our QR-of-B' approach, so B is left as-is for now.
    } else {
        // B is M×N. QR factorization: B = Q*R. Copy upper MN×N of R to U.
        let b_mat = DMatrix::from_fn(m, n, |i, j| b[i + j * ldb]);
        let qr = b_mat.qr();
        let r = qr.r();

        for j in 0..n {
            for i in 0..=min(j, mn - 1) {
                u[i + j * ldu] = r[(i, j)];
            }
            for i in mn..n {
                u[i + j * ldu] = ZERO;
            }
        }
        if m < n {
            // Zero trailing (N-M)×(N-M) at U(M+1,M+1)
            for j in m..n {
                for i in m..n {
                    u[i + j * ldu] = ZERO;
                }
            }
        }
    }

    // Solve canonical Lyapunov for U
    let info = sb03ot(discr, ltrans, n, a, lda, u, ldu, scale, dwork);
    if info != 0 && info != 1 {
        return info;
    }

    // Make diagonal of U non-negative
    for j in 0..n {
        if u[j + j * ldu] < ZERO {
            for i in 0..=j {
                u[i + j * ldu] = -u[i + j * ldu];
            }
        }
    }

    if !dwork.is_empty() {
        dwork[0] = (4 * n) as f64;
    }
    info
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03ou_1x1_cont() {
        let n = 1usize;
        let m = 1usize;
        let a = [-1.0];
        let mut b = [1.0];
        let mut tau = [0.0];
        let mut u = [0.0];
        let mut scale = 0.0;
        let mut dwork = vec![0.0; 4 * n];
        let info = sb03ou(
            false,
            false,
            n,
            m,
            &a,
            1,
            &mut b,
            1,
            &mut tau,
            &mut u,
            1,
            &mut scale,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert!((scale - 1.0).abs() < 1e-10);
        assert!(u[0].abs() > 0.0);
    }
}
