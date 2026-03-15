//! AB01MD — Find controllable realization for single-input system (SLICOT AB01MD)
//!
//! Reduces (A, B) to orthogonal canonical form using Householder transformations:
//! B to one non-zero element, A to upper Hessenberg; returns NCONT (controllable order).

use nalgebra::{DMatrix, DVector};

/// Whether to accumulate the orthogonal transformation matrix Z.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobZ {
    /// Do not form Z.
    No,
    /// Store transformations in factored form (Z and TAU); not yet fully supported.
    Factored,
    /// Z is set to identity on entry; on exit contains the transformation matrix Z.
    Init,
}

impl JobZ {
    fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'N' => Some(JobZ::No),
            'F' => Some(JobZ::Factored),
            'I' => Some(JobZ::Init),
            _ => None,
        }
    }
}

fn to_ascii_upper(c: char) -> char {
    if c >= 'a' && c <= 'z' {
        ((c as u8) - b'a' + b'A') as char
    } else {
        c
    }
}

/// Max (infinity) norm of a matrix (max absolute row sum for 1-norm of columns).
fn norm_inf_matrix(a: &DMatrix<f64>) -> f64 {
    let (n, m) = (a.nrows(), a.ncols());
    if n == 0 || m == 0 {
        return 0.0;
    }
    let mut max = 0.0_f64;
    for j in 0..m {
        let mut s = 0.0;
        for i in 0..n {
            s += a[(i, j)].abs();
        }
        if s > max {
            max = s;
        }
    }
    max
}

/// Frobenius norm.
fn norm_frobenius(a: &DMatrix<f64>) -> f64 {
    a.norm()
}

/// 1-norm of a vector (sum of absolute values).
fn norm_1_vec(b: &DVector<f64>) -> f64 {
    b.iter().map(|x| x.abs()).sum()
}

/// Scale leading (m,n) block of A by 1/scale (for undo: multiply by scale).
fn scale_block(a: &mut DMatrix<f64>, m: usize, n: usize, scale: f64) {
    if scale == 0.0 {
        return;
    }
    let s = 1.0 / scale;
    for j in 0..n {
        for i in 0..m {
            a[(i, j)] *= s;
        }
    }
}

/// Undo scale: multiply block by scale.
fn unscale_block(a: &mut DMatrix<f64>, m: usize, n: usize, scale: f64) {
    for j in 0..n {
        for i in 0..m {
            a[(i, j)] *= scale;
        }
    }
}

/// Generate Householder reflector: overwrites x with [beta, v[1], v[2], ...], returns tau.
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
    let n2 = (x.view((1, 0), (n - 1, 1))).norm_squared();
    2.0 / (1.0 + n2)
}

/// Apply Householder from left to rows [row_start..row_start+len+1] of columns [col_start..].
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

/// Apply Householder from right to columns [col_start..col_start+len+1] of all rows.
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

/// Reduce A to upper Hessenberg in place (DGEHRD-style). Stores reflectors in lower part of A;
/// tau[0..n-1] receives the scalar factors. ilo=0, ihi=n-1.
fn hessenberg_reduce(a: &mut DMatrix<f64>, tau: &mut [f64], work: &mut [f64]) {
    let n = a.nrows();
    if n <= 2 {
        return;
    }
    for j in 0..n - 2 {
        let row_start = j + 1; // used in apply_householder_left/right
        let len = n - j - 2; // reflector length in elements below diagonal = (n-1)-(j+1) = n-j-2, so len+1 = n-j-1
        let nj = n - j - 1;   // number of rows from j to n-1
        let mut x = DVector::from_fn(nj, |i, _| a[(row_start + i, j)]);
        let tau_j = householder_gen(&mut x);
        tau[j] = tau_j;
        if tau_j == 0.0 {
            continue;
        }
        // Store reflector in A(j+2:n, j) — x[0] is beta, x[1..] is v(2..)
        a[(j + 1, j)] = x[0];
        for i in 1..nj {
            a[(j + 1 + i, j)] = x[i];
        }
        let v: Vec<f64> = (0..nj).map(|i| if i == 0 { 1.0 } else { a[(j + 1 + i, j)] }).collect();
        // Apply from left to A(j+1:n, j:n)
        apply_householder_left(a, row_start, len, j, n - j, &v, tau_j, work);
        // Apply from right to A(0:n, j+1:n)
        apply_householder_right(a, n, j + 1, len, &v, tau_j, work);
        // Restore stored reflector (left mul may have overwritten)
        a[(j + 1, j)] = x[0];
        for i in 1..nj {
            a[(j + 1 + i, j)] = x[i];
        }
    }
}

/// Form orthogonal Q from Hessenberg reflectors stored in A (lower part) and tau, into Z (Z := Q).
/// Assumes Z is n×n; we apply reflectors from the right to identity.
fn form_q_hessenberg(a: &DMatrix<f64>, tau: &[f64], z: &mut DMatrix<f64>, work: &mut [f64]) {
    let n = a.nrows();
    z.fill(0.0);
    for i in 0..n {
        z[(i, i)] = 1.0;
    }
    if n <= 2 {
        return;
    }
    for j in (0..n - 2).rev() {
        let _row_start = j + 1;
        let nj = n - j - 1;
        let len = nj - 1;
        let tau_j = tau[j];
        let mut v = vec![0.0; nj];
        v[0] = 1.0;
        for i in 1..nj {
            v[i] = a[(j + 1 + i, j)];
        }
        // Apply (I - tau*v*v') from the right to Z(:, j+1:n)
        apply_householder_right(z, n, j + 1, len, &v, tau_j, work);
    }
}

/// Find controllable realization for single-input system (A, B). Overwrites A and B.
///
/// # Arguments
/// * `jobz` - No: do not form Z; Init: return transformation matrix in Z; Factored: store in factored form (Z/TAU).
/// * `a` - N×N state matrix (overwritten with upper Hessenberg canonical form).
/// * `b` - N×1 input vector (overwritten; only B(1) non-zero on exit).
/// * `z` - If JobZ::Init (or Factored), N×N matrix; for Init, set to identity on entry to get Z on exit.
/// * `tol` - Tolerance for controllability; if <= 0, default is used (N*eps*max(||A||_F, ||B||_1)).
/// * `ncont` - On exit, order of controllable realization (0..=N).
///
/// # Returns
/// 0 on success; < 0 if -i indicates the i-th argument was invalid.
pub fn ab01md(
    jobz: JobZ,
    a: &mut DMatrix<f64>,
    b: &mut DVector<f64>,
    mut z: Option<&mut DMatrix<f64>>,
    tol: f64,
    ncont: &mut usize,
) -> i32 {
    let n = a.nrows();
    if a.ncols() != n {
        return -3;
    }
    if b.nrows() != n {
        return -5;
    }
    let ljobi = jobz == JobZ::Init;
    let ljobf = jobz == JobZ::Factored;
    let ljobz = ljobi || ljobf;
    if !ljobz && jobz != JobZ::No {
        return -1;
    }
    if ljobz && (z.is_none() || z.as_ref().map(|z| z.nrows() != n || z.ncols() != n).unwrap_or(true)) {
        return -7;
    }

    *ncont = 0;
    if n == 0 {
        return 0;
    }

    let eps = f64::EPSILON;
    let mut dwork = vec![0.0_f64; n.max(1)];

    let anorm = norm_inf_matrix(a);
    let bnorm = norm_inf_matrix(&b.view((0, 0), (n, 1)).into_owned());
    if bnorm == 0.0 {
        if ljobi {
            if let Some(ref mut zz) = z {
                zz.fill(0.0);
                for i in 0..n {
                    zz[(i, i)] = 1.0;
                }
            }
        }
        return 0;
    }

    scale_block(a, n, n, anorm);
    for i in 0..n {
        b[i] /= bnorm;
    }

    let fanorm = norm_frobenius(a);
    let fbnorm = norm_1_vec(b);
    let mut toldef = tol;
    if toldef <= 0.0 {
        toldef = (n as f64) * eps * fanorm.max(fbnorm);
    }

    if fbnorm <= toldef {
        unscale_block(a, n, n, anorm);
        for i in 0..n {
            b[i] *= bnorm;
        }
        if ljobi {
            if let Some(ref mut zz) = z {
                zz.fill(0.0);
                for i in 0..n {
                    zz[(i, i)] = 1.0;
                }
            }
        }
        return 0;
    }

    let mut tau = vec![0.0_f64; n];
    let b1 = if n > 1 {
        let mut x = DVector::from_fn(n, |i, _| b[i]);
        let h = householder_gen(&mut x);
        let b1_val = x[0];
        b[0] = 1.0;
        for i in 1..n {
            b[i] = x[i];
        }
        let v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { b[i] }).collect();
        apply_householder_left(a, 0, n - 2, 0, n, &v, h, &mut dwork);
        apply_householder_right(a, n, 0, n - 2, &v, h, &mut dwork);
        b[0] = b1_val;
        tau[0] = h;
        b1_val
    } else {
        b[0]
    };

    hessenberg_reduce(a, &mut tau[1..], &mut dwork);

    if ljobz {
        if let Some(ref mut zz) = z {
            if n > 1 {
                zz[(0, 0)] = 1.0;
                for i in 1..n {
                    zz[(i, 0)] = b[i];
                }
            }
            if n > 2 {
                for j in 1..n - 1 {
                    for i in j + 2..n {
                        zz[(i, j)] = a[(i, j - 1)];
                    }
                }
            }
            if ljobi {
                form_q_hessenberg(a, &tau[1..], zz, &mut dwork);
                // Z should be Z1 * Z2. We have Z2 in zz now. Z1 = I - tau[0]*v*v' with v = [1, b[1], ..., b[n-1]].
                // So Z = Z1 * Z2: zz <- Z1 * zz. Apply Z1 from the left: (I - tau[0]*v*v') * zz.
                if n > 1 {
                    let v: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { b[i] }).collect();
                    let h = tau[0];
                    for j in 0..n {
                        let mut dot = 0.0;
                        for i in 0..n {
                            dot += v[i] * zz[(i, j)];
                        }
                        let s = h * dot;
                        for i in 0..n {
                            zz[(i, j)] -= v[i] * s;
                        }
                    }
                }
            }
        }
    }

    // Zero below first subdiagonal and B(2:n)
    if n > 2 {
        for j in 0..n - 2 {
            for i in j + 3..n {
                a[(i, j)] = 0.0;
            }
        }
    }
    for i in 1..n {
        b[i] = 0.0;
    }

    if tol <= 0.0 {
        toldef = (n as f64) * eps * fanorm.max(b1.abs());
    }
    let mut j = 0_usize;
    while j < n - 1 {
        if a[(j + 1, j)].abs() <= toldef {
            break;
        }
        j += 1;
    }
    *ncont = j + 1;
    if j + 1 < n {
        a[(j + 1, j)] = 0.0;
    }

    unscale_block(a, *ncont, *ncont, anorm);
    b[0] *= bnorm;
    if *ncont < n {
        for c in *ncont..n {
            for r in 0..n {
                a[(r, c)] *= anorm;
            }
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab01md_n1() {
        let mut a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DVector::from_row_slice(&[1.0]);
        let mut ncont = 0_usize;
        let info = ab01md(JobZ::No, &mut a, &mut b, None, 0.0, &mut ncont);
        assert_eq!(info, 0);
        assert_eq!(ncont, 1);
    }

    #[test]
    fn test_ab01md_n2_controllable() {
        // (A,B) that reduces to Hessenberg with non-zero subdiagonal -> NCONT=2
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let mut b = DVector::from_row_slice(&[1.0, 0.0]);
        let mut ncont = 0_usize;
        let info = ab01md(JobZ::No, &mut a, &mut b, None, 0.0, &mut ncont);
        assert_eq!(info, 0);
        assert_eq!(ncont, 2);
        assert!(b[1].abs() < 1e-10);
    }

    #[test]
    fn test_ab01md_zero_b() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut b = DVector::from_row_slice(&[0.0, 0.0]);
        let mut ncont = 99_usize;
        let info = ab01md(JobZ::No, &mut a, &mut b, None, 0.0, &mut ncont);
        assert_eq!(info, 0);
        assert_eq!(ncont, 0);
    }

    #[test]
    fn test_ab01md_jobz_init() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let mut b = DVector::from_row_slice(&[1.0, 0.0]);
        let mut z = DMatrix::identity(2, 2);
        let mut ncont = 0_usize;
        let info = ab01md(JobZ::Init, &mut a, &mut b, Some(&mut z), 0.0, &mut ncont);
        assert_eq!(info, 0);
        assert_eq!(ncont, 2);
        assert!(b[1].abs() < 1e-10);
    }
}
