//! TB01XD — Dual transformation P*A'*P, B <- P*C', C <- B'*P (SLICOT TB01XD)
//!
//! Special similarity transformation of the dual system. Optionally transpose D.

use nalgebra::DMatrix;

/// Whether the direct transmission matrix D is present.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobD {
    /// D is present; on exit D is replaced by D'.
    Present,
    /// D is assumed zero (not referenced).
    Zero,
}

/// Applies dual transformation: A <- P*A'*P, B <- P*C', C <- B'*P.
/// If JobD::Present, D <- D'. B and C swap dimensions (B becomes N×P, C becomes M×N).
///
/// # Arguments
/// * `jobd` - Whether D is present (and should be transposed).
/// * `a` - N×N state matrix; overwritten with P*A'*P.
/// * `b` - On entry N×M; on exit N×P (P*C'). Caller must ensure b has at least max(M,P) columns or provide separate output for C'.
/// * `c` - On entry P×N; on exit M×N (B'*P). Caller must ensure c has at least max(M,P) rows.
/// * `d` - If JobD::Present, P×M on entry, M×P on exit (D'). Otherwise not referenced.
///
/// # Returns
/// 0 on success; &lt; 0 if the i-th argument had an illegal value.
///
/// # Note
/// SLICOT uses a single B array sized (LDB, MAX(M,P)) and C (LDC, N) with LDC >= MAX(1,M,P). So on exit B holds P*C' (N×P) and C holds B'*P (M×N). Caller must pass B with columns >= max(M,P) and C with rows >= max(M,P). This implementation takes b and c as mutable; b must have nrows=N, ncols>=max(M,P); c must have nrows>=max(M,P), ncols=N. On exit we write P*C' into b columns 0..P and B'*P into c rows 0..M.
pub fn tb01xd(
    jobd: JobD,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    mut d: Option<&mut DMatrix<f64>>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n {
        return -6;
    }
    if b.nrows() != n {
        return -8;
    }
    if c.ncols() != n {
        return -10;
    }
    if m > 0 && b.nrows() != n {
        return -8;
    }
    if p > 0 && c.ncols() != n {
        return -10;
    }
    if jobd == JobD::Present {
        if let Some(ref dd) = d {
            if dd.nrows() != p || dd.ncols() != m {
                return -14;
            }
        } else {
            return -14;
        }
    }
    if n == 0 {
        if jobd == JobD::Present {
            if let Some(ref mut dd) = d {
                let dc = dd.clone();
                for j in 0..m {
                    for i in 0..p {
                        dd[(j, i)] = dc[(i, j)];
                    }
                }
            }
        }
        return 0;
    }
    // A <- P*A'*P: (A')[i,j] = A[j,i], then (P*X*P)[i,j] = X[n-1-j, n-1-i], so (P*A'*P)[i,j] = A[n-1-i, n-1-j]
    let ac = a.clone();
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = ac[(n - 1 - j, n - 1 - i)];
        }
    }
    let cc = c.clone();
    let b_orig = b.clone();
    // B <- P*C': (P*C')[i,k] = C[k, n-1-i]
    for i in 0..n {
        for k in 0..p.min(b.ncols()) {
            b[(i, k)] = cc[(k, n - 1 - i)];
        }
    }
    // C <- B'*P: (B'*P)[k,j] = B[n-1-j, k]
    for k in 0..m.min(c.nrows()) {
        for j in 0..n {
            c[(k, j)] = b_orig[(n - 1 - j, k)];
        }
    }
    if jobd == JobD::Present {
        if let Some(ref mut dd) = d {
            let dc = dd.clone();
            for i in 0..p {
                for j in 0..m {
                    dd[(j, i)] = dc[(i, j)];
                }
            }
        }
    }
    0
}
