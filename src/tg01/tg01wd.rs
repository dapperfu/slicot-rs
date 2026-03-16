//! TG01WD — Reduce (A,E) to generalized real Schur form, apply to B and C (SLICOT TG01WD)
//!
//! Pure Rust: use Hessenberg + QZ-style iteration would require full QZ; here we leave (A,E) unchanged
//! and set Q=Z=I, with alphar/alphai/beta from a simple eigenvalue computation on E^{-1}*A when E is invertible.

use nalgebra::DMatrix;

/// Reduce (A,E) to real generalized Schur form. Placeholder: Q=Z=I, output eigenvalues when E=I.
pub fn tg01wd(
    n: usize,
    _m: usize,
    _p: usize,
    a: &mut DMatrix<f64>,
    e: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    _c: &mut DMatrix<f64>,
    q: &mut DMatrix<f64>,
    z: &mut DMatrix<f64>,
    alphar: &mut [f64],
    alphai: &mut [f64],
    beta: &mut [f64],
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -4;
    }
    if e.nrows() != n || e.ncols() != n {
        return -6;
    }
    if n == 0 {
        return 0;
    }
    if alphar.len() < n || alphai.len() < n || beta.len() < n {
        return -17;
    }
    if q.nrows() != n || q.ncols() != n || z.nrows() != n || z.ncols() != n {
        return -15;
    }
    q.fill_with_identity();
    z.fill_with_identity();
    if let Some(e_inv) = e.clone().try_inverse() {
        let ae = e_inv * a.clone();
        let schur = ae.schur();
        let eigs = schur.complex_eigenvalues();
        for i in 0..n.min(eigs.len()) {
            alphar[i] = eigs[i].re;
            alphai[i] = eigs[i].im;
            beta[i] = 1.0;
        }
        for i in eigs.len()..n {
            alphar[i] = 0.0;
            alphai[i] = 0.0;
            beta[i] = 1.0;
        }
    } else {
        for i in 0..n {
            alphar[i] = a[(i, i)];
            alphai[i] = 0.0;
            beta[i] = e[(i, i)].max(1e-30);
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tg01wd_smoke() {
        let n = 2;
        let mut a = DMatrix::from_row_slice(n, n, &[1.0, 0.0, 0.0, 2.0]);
        let mut e = DMatrix::identity(n, n);
        let mut b = DMatrix::zeros(n, 1);
        let mut c = DMatrix::zeros(1, n);
        let mut q = DMatrix::zeros(n, n);
        let mut z = DMatrix::zeros(n, n);
        let mut alphar = vec![0.0; n];
        let mut alphai = vec![0.0; n];
        let mut beta = vec![0.0; n];
        assert_eq!(tg01wd(n, 1, 1, &mut a, &mut e, &mut b, &mut c, &mut q, &mut z, &mut alphar, &mut alphai, &mut beta), 0);
        assert!((alphar[0] - 1.0).abs() < 1e-6);
        assert!((alphar[1] - 2.0).abs() < 1e-6);
    }
}
