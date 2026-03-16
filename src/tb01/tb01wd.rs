//! TB01WD — Reduce A to real Schur form, apply to B and C (SLICOT TB01WD)
//!
//! A <- U'*A*U (real Schur), B <- U'*B, C <- C*U; outputs eigenvalues in WR, WI.

use nalgebra::DMatrix;

/// Reduces A to upper real Schur form by orthogonal similarity U'*A*U and applies
/// the same transformation to B and C. Fills WR, WI with real and imaginary parts of eigenvalues.
///
/// # Returns
/// 0 on success; < 0 invalid argument; > 0 QR failed to converge (INFO = first unconverged index).
pub fn tb01wd(
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    u: &mut DMatrix<f64>,
    wr: &mut [f64],
    wi: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n {
        return -5;
    }
    if b.nrows() != n || c.ncols() != n {
        return -7;
    }
    if u.nrows() != n || u.ncols() != n {
        return -10;
    }
    if wr.len() < n || wi.len() < n {
        return -11;
    }
    if n == 0 {
        return 0;
    }
    let schur = match a.clone().try_schur(1e-14, 100) {
        Some(s) => s,
        None => return 1,
    };
    let eigs = schur.complex_eigenvalues();
    let (q, r) = schur.unpack();
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = r[(i, j)];
        }
    }
    for i in 0..n {
        for j in 0..n {
            u[(i, j)] = q[(i, j)];
        }
    }
    let qt = q.transpose();
    let b_new = &qt * b.clone();
    let c_new = c.clone() * &q;
    for i in 0..n {
        for j in 0..b.ncols() {
            b[(i, j)] = b_new[(i, j)];
        }
    }
    for i in 0..c.nrows() {
        for j in 0..n {
            c[(i, j)] = c_new[(i, j)];
        }
    }
    for (i, z) in eigs.iter().enumerate().take(n) {
        wr[i] = z.re;
        wi[i] = z.im;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb01wd_small() {
        let n = 2;
        let m = 1;
        let p = 1;
        let mut a = DMatrix::from_row_slice(n, n, &[1.0, 1.0, 0.0, 2.0]);
        let mut b = DMatrix::from_row_slice(n, m, &[1.0, 0.0]);
        let mut c = DMatrix::from_row_slice(p, n, &[1.0, 0.0]);
        let mut u = DMatrix::zeros(n, n);
        let mut wr = vec![0.0; n];
        let mut wi = vec![0.0; n];
        let info = tb01wd(&mut a, &mut b, &mut c, &mut u, &mut wr, &mut wi);
        assert_eq!(info, 0);
        assert!((a[(1, 0)]).abs() < 1e-10);
    }
}
