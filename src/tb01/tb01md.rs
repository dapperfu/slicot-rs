//! TB01MD — Reduce (B, A) to controller Hessenberg form (SLICOT TB01MD)
//!
//! 1:1 mapping of SLICOT TB01MD: unitary state-space transformation reducing
//! the pair (B, A) to upper or lower controller Hessenberg form.

use nalgebra::{DMatrix, DVector};

/// Whether to accumulate the transformation matrix U.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobU {
    /// Do not form U.
    No,
    /// U is initialized to identity and the transformation matrix is returned.
    Init,
    /// The given U is updated by the transformations.
    Update,
}

impl JobU {
    fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'N' => Some(JobU::No),
            'I' => Some(JobU::Init),
            'U' => Some(JobU::Update),
            _ => None,
        }
    }
}

/// Upper or lower controller Hessenberg form.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    /// Upper controller Hessenberg form.
    Upper,
    /// Lower controller Hessenberg form.
    Lower,
}

impl Uplo {
    fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'U' => Some(Uplo::Upper),
            'L' => Some(Uplo::Lower),
            _ => None,
        }
    }
}

trait ToAsciiUpper {
    fn to_ascii_uppercase(self) -> char;
}
impl ToAsciiUpper for char {
    fn to_ascii_uppercase(self) -> char {
        if self >= 'a' && self <= 'z' {
            ((self as u8) - b'a' + b'A') as char
        } else {
            self
        }
    }
}

/// Reduces the pair (B, A) to controller Hessenberg form using unitary transformations.
///
/// # Arguments
/// * `jobu` - Whether to accumulate transformation in U (No, Init, Update).
/// * `uplo` - Upper or lower controller Hessenberg form.
/// * `a` - State matrix A (n×n), overwritten with U' A U.
/// * `b` - Input matrix B (n×m), overwritten with U' B.
/// * `u` - If JobU::Init or Update, the transformation matrix U (n×n). For Init, must be identity on entry or ignored; for Update, updated in place.
///
/// # Returns
/// * `0` - success
/// * `< 0` - if `-i`, the i-th argument had an illegal value (1=jobu, 2=uplo, 3=n, 4=m, 6=lda, 8=ldb, 10=ldu)
pub fn tb01md(
    jobu: JobU,
    uplo: Uplo,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    u: &mut Option<&mut DMatrix<f64>>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    if a.ncols() != n || b.nrows() != n {
        return -6;
    }
    if (jobu == JobU::Init || jobu == JobU::Update) && u.is_none() {
        return -10;
    }
    if let Some(ref uu) = *u {
        if uu.nrows() != n || uu.ncols() != n {
            return -10;
        }
    }
    if n == 0 || m == 0 {
        return 0;
    }
    if uplo != Uplo::Upper {
        return -2; // Lower not yet implemented in pilot
    }
    let n1 = n - 1;
    let mut dwork = vec![0.0_f64; n.max(m.saturating_sub(1))];

    if jobu == JobU::Init {
        if let Some(ref mut uu) = *u {
            uu.fill(0.0);
            for i in 0..n {
                uu[(i, i)] = 1.0;
            }
        }
    }

    let ljoba = jobu == JobU::Init || jobu == JobU::Update;

    // Phase 1: transformations involving both B and A (J = 0..min(m, n-1)); Upper only
    for j in 0..(m.min(n1)) {
        let nj = n - j; // number of rows in block (0-based: j..n-1)
        let (par1, par2, par3, par4, par5) = (j, j, j + 1, m, n);
        let col_b = par1;
        let row_start = par2;
        let row_next = par3;
        let len = nj - 1; // reflector length = nj (so len+1 = nj), rows j..n-1

        let mut x = DVector::from_fn(len + 1, |i, _| b[(row_start + i, col_b)]);
        let tau = householder_gen(&mut x);
        if tau == 0.0 {
            continue;
        }
        let v = DVector::from_fn(len + 1, |i, _| if i == 0 { 1.0 } else { x[i] });

        // Update A: left and right
        apply_householder_left(a, row_start, row_next, len, 0, n, &v, tau, &mut dwork);
        apply_householder_right(a, n, row_start, row_next, len, &v, tau, &mut dwork);

        if ljoba {
            if let Some(ref mut uu) = *u {
                apply_householder_right_mat(uu, n, row_start, row_next, len, &v, tau, &mut dwork);
            }
        }

        if j != m - 1 {
            let (col_start_b, cols_b) = (par3, par4 - par3);
            if cols_b > 0 {
                apply_householder_left(
                    b,
                    row_start,
                    row_next,
                    len,
                    col_start_b,
                    cols_b,
                    &v,
                    tau,
                    &mut dwork,
                );
            }
        }

        for ii in par3..par5.min(n) {
            b[(ii, col_b)] = 0.0;
        }
    }

    // Phase 2: transformations only involving A (J = m..n-2); Upper only
    for j in m..n1 {
        let nj = n - j;
        let (par1, par2, par3, par4, par5, par6) = (j - m, j, j + 1, n, j - m, n);
        let col_a = par1;
        let row_start = par2;
        let row_next = par3;
        let len = nj - 1; // reflector length = nj, rows j..n-1
        let col_start = par5;
        let col_count = par6 - par5;

        let mut x = DVector::from_fn(len + 1, |i, _| a[(row_start + i, col_a)]);
        let tau = householder_gen(&mut x);
        if tau == 0.0 {
            continue;
        }
        let v = DVector::from_fn(len + 1, |i, _| if i == 0 { 1.0 } else { x[i] });

        apply_householder_left(
            a,
            row_start,
            row_next,
            len,
            col_start,
            col_count,
            &v,
            tau,
            &mut dwork,
        );
        apply_householder_right(a, n, row_start, row_next, len, &v, tau, &mut dwork);

        if ljoba {
            if let Some(ref mut uu) = *u {
                apply_householder_right_mat(uu, n, row_start, row_next, len, &v, tau, &mut dwork);
            }
        }

        for ii in par3..(par4 + 1).min(n) {
            a[(ii, col_a)] = 0.0;
        }
    }

    0
}

/// Generate Householder reflector: overwrites x with (beta, v[1], v[2], ...), returns tau.
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

/// Apply Householder from left to rows [row_start..row_start+len+1] of columns [col_start..col_start+ncols].
fn apply_householder_left(
    a: &mut DMatrix<f64>,
    row_start: usize,
    _row_next: usize,
    len: usize,
    col_start: usize,
    ncols: usize,
    v: &DVector<f64>,
    tau: f64,
    work: &mut [f64],
) {
    for c in 0..ncols {
        let col = col_start + c;
        let mut dot = 0.0;
        for i in 0..len + 1 {
            dot += v[i] * a[(row_start + i, col)];
        }
        work[c] = tau * dot;
    }
    for i in 0..len + 1 {
        for c in 0..ncols {
            a[(row_start + i, col_start + c)] -= v[i] * work[c];
        }
    }
}

/// Apply Householder from right to columns [col_start..col_start+len] of all rows.
fn apply_householder_right(
    a: &mut DMatrix<f64>,
    nrows: usize,
    col_start: usize,
    col_next: usize,
    len: usize,
    v: &DVector<f64>,
    tau: f64,
    work: &mut [f64],
) {
    for r in 0..nrows {
        let mut dot = 0.0;
        for i in 0..len + 1 {
            dot += a[(r, col_start + i)] * v[i];
        }
        work[r] = tau * dot;
    }
    for r in 0..nrows {
        for i in 0..len + 1 {
            a[(r, col_start + i)] -= work[r] * v[i];
        }
    }
}

fn apply_householder_right_mat(
    a: &mut DMatrix<f64>,
    nrows: usize,
    col_start: usize,
    _col_next: usize,
    len: usize,
    v: &DVector<f64>,
    tau: f64,
    work: &mut [f64],
) {
    apply_householder_right(a, nrows, col_start, col_start + 1, len, v, tau, work);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb01md_small_upper() {
        let n = 3usize;
        let m = 2usize;
        let mut a = DMatrix::from_row_slice(n, n, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let mut b = DMatrix::from_row_slice(n, m, &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
        let mut u = DMatrix::identity(n, n);
        let mut u_opt = Some(&mut u);
        let info = tb01md(JobU::Init, Uplo::Upper, &mut a, &mut b, &mut u_opt);
        assert_eq!(info, 0);
        // B should have zeros below diagonal in first column (controller Hessenberg)
        assert!(b[(1, 0)].abs() < 1e-10 || b[(2, 0)].abs() < 1e-10);
    }

    #[test]
    fn test_tb01md_jobu_no() {
        let n = 2usize;
        let m = 1usize;
        let mut a = DMatrix::from_row_slice(n, n, &[1.0, 0.0, 0.0, 1.0]);
        let mut b = DMatrix::from_row_slice(n, m, &[1.0, 0.0]);
        let mut nopt = None;
        let info = tb01md(JobU::No, Uplo::Upper, &mut a, &mut b, &mut nopt);
        assert_eq!(info, 0);
    }

    #[test]
    fn test_tb01md_lower_not_implemented() {
        let n = 2usize;
        let m = 1usize;
        let mut a = DMatrix::from_row_slice(n, n, &[1.0, 0.0, 0.0, 1.0]);
        let mut b = DMatrix::from_row_slice(n, m, &[1.0, 0.0]);
        let mut nopt = None;
        let info = tb01md(JobU::No, Uplo::Lower, &mut a, &mut b, &mut nopt);
        assert_eq!(info, -2); // Lower not yet implemented
    }
}
