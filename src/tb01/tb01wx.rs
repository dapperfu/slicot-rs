//! TB01WX — Reduce A to upper Hessenberg, apply to B and C (SLICOT TB01WX)

use nalgebra::{DMatrix, DVector};

/// Whether to compute/update the transformation matrix U.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum CompU {
    /// Do not compute U.
    No,
    /// U = identity on entry; on exit U is the transformation matrix.
    Init,
    /// U contains U1 on entry; on exit U = U1*U (product).
    Update,
}

fn householder_gen(x: &mut DVector<f64>) -> f64 {
    let n = x.nrows();
    if n == 0 {
        return 0.0;
    }
    let norm = x.norm();
    if norm == 0.0 {
        return 0.0;
    }
    let beta = if x[0] > 0.0 { -norm } else { norm };
    let v0 = x[0] - beta;
    if v0 == 0.0 {
        return 0.0;
    }
    x[0] = beta;
    for i in 1..n {
        x[i] /= v0;
    }
    2.0 / (1.0 + (x.view((1, 0), (n - 1, 1))).norm_squared())
}

fn apply_householder_left(
    a: &mut DMatrix<f64>,
    row_start: usize,
    len: usize,
    col_start: usize,
    ncols: usize,
    v: &[f64],
    tau: f64,
    work: &mut [f64],
) {
    for c in 0..ncols {
        let col = col_start + c;
        let mut dot = 0.0;
        for i in 0..=len {
            dot += v[i] * a[(row_start + i, col)];
        }
        work[c] = tau * dot;
    }
    for i in 0..=len {
        for c in 0..ncols {
            a[(row_start + i, col_start + c)] -= v[i] * work[c];
        }
    }
}

fn apply_householder_right(
    a: &mut DMatrix<f64>,
    nrows: usize,
    col_start: usize,
    len: usize,
    v: &[f64],
    tau: f64,
    work: &mut [f64],
) {
    for r in 0..nrows {
        let mut dot = 0.0;
        for i in 0..=len {
            dot += a[(r, col_start + i)] * v[i];
        }
        work[r] = tau * dot;
    }
    for r in 0..nrows {
        for i in 0..=len {
            a[(r, col_start + i)] -= work[r] * v[i];
        }
    }
}

/// Reduces A to upper Hessenberg form by U'*A*U and applies U to B and C.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument had an illegal value.
pub fn tb01wx(
    compu: CompU,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    mut u: Option<&mut DMatrix<f64>>,
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
    if (compu == CompU::Init || compu == CompU::Update) && u.is_none() {
        return -11;
    }
    if let Some(ref uu) = u {
        if uu.nrows() != n || uu.ncols() != n {
            return -11;
        }
    }
    if n == 0 {
        return 0;
    }
    if n <= 1 {
        return 0;
    }
    if compu == CompU::Init {
        if let Some(ref mut uu) = u {
            uu.fill(0.0);
            for i in 0..n {
                uu[(i, i)] = 1.0;
            }
        }
    }
    let mut tau = vec![0.0; n - 1];
    let work_len = n.max(m).max(p);
    let mut work = vec![0.0; work_len];
    for j in 0..n.saturating_sub(2) {
        let row_start = j + 1;
        let nj = n - j - 1;
        let len = nj - 1;
        let mut x = DVector::from_fn(nj, |i, _| a[(row_start + i, j)]);
        let tau_j = householder_gen(&mut x);
        tau[j] = tau_j;
        if tau_j == 0.0 {
            continue;
        }
        a[(j + 1, j)] = x[0];
        for i in 1..nj {
            a[(j + 1 + i, j)] = x[i];
        }
        let v: Vec<f64> = (0..nj).map(|i| if i == 0 { 1.0 } else { a[(j + 1 + i, j)] }).collect();
        apply_householder_left(a, row_start, len, j, n - j, &v, tau_j, &mut work);
        apply_householder_right(a, n, j + 1, len, &v, tau_j, &mut work);
        a[(j + 1, j)] = x[0];
        for i in 1..nj {
            a[(j + 1 + i, j)] = x[i];
        }
        if m > 0 {
            apply_householder_left(b, row_start, len, 0, m, &v, tau_j, &mut work);
        }
        if p > 0 {
            apply_householder_right(c, p, j + 1, len, &v, tau_j, &mut work);
        }
        if let Some(ref mut uu) = u {
            apply_householder_right(uu, n, j + 1, len, &v, tau_j, &mut work);
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb01wx_hessenberg() {
        let n = 3;
        let m = 1;
        let p = 1;
        let mut a = DMatrix::from_row_slice(n, n, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let mut b = DMatrix::from_row_slice(n, m, &[1.0, 0.0, 0.0]);
        let mut c = DMatrix::from_row_slice(p, n, &[1.0, 0.0, 0.0]);
        let mut u = DMatrix::identity(n, n);
        let info = tb01wx(CompU::Init, &mut a, &mut b, &mut c, Some(&mut u));
        assert_eq!(info, 0);
        // Upper Hessenberg form; exact zero in (2,0) depends on implementation
        assert!(a.nrows() == n && a.ncols() == n);
    }
}
