//! TB05AD — Frequency response matrix of state-space (A,B,C) (SLICOT TB05AD)
//!
//! G(freq) = C * inv(freq*I - A) * B at complex frequency freq.

use nalgebra::{linalg::LU, DMatrix};
use num_complex::Complex64;

/// Whether to balance A and/or compute eigenvalues and condition estimate.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Baleig {
    /// No balance, no eigenvalues, no RCOND.
    N,
    /// No balance, compute RCOND only.
    C,
    /// Balance and compute eigenvalues (and optionally RCOND with A).
    B,
    /// Balance and compute eigenvalues and RCOND.
    A,
}

/// Whether A is general or already upper Hessenberg.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Inita {
    /// A is general.
    G,
    /// A is already upper Hessenberg.
    H,
}

/// Computes the frequency response matrix G(freq) = C * inv(freq*I - A) * B.
///
/// # Arguments
/// * `baleig` - Balance/eigenvalue/RCOND options (this implementation supports N and C).
/// * `inita` - General or Hessenberg A (only G supported: direct solve).
/// * `a` - State matrix (n×n), may be overwritten if INITA='G' and BALEIG requests transform.
/// * `b` - Input matrix (n×m).
/// * `c` - Output matrix (p×n).
/// * `freq` - Complex frequency at which to evaluate.
/// * `g` - Output: G(freq) (p×m) complex.
/// * `rcond` - If BALEIG='C' or 'A', receives reciprocal condition of (freq*I - A).
/// * `hinvb` - If provided, receives inv(freq*I - A)*B (n×m) complex.
/// * `evre`, `evim` - If BALEIG requests eigenvalues, real/imag parts (not implemented).
///
/// # Returns
/// * `0` - success
/// * `< 0` - invalid argument
/// * `1` - too many iterations (eigenvalue)
/// * `2` - freq too near eigenvalue or RCOND < eps
pub fn tb05ad(
    baleig: Baleig,
    _inita: Inita,
    _a: &mut DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    freq: Complex64,
    g: &mut DMatrix<Complex64>,
    rcond: Option<&mut f64>,
    mut hinvb: Option<&mut DMatrix<Complex64>>,
    _evre: Option<&mut [f64]>,
    _evim: Option<&mut [f64]>,
) -> i32 {
    let n = _a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if _a.ncols() != n {
        return -8;
    }
    if b.nrows() != n || c.ncols() != n {
        return -9;
    }
    if g.nrows() != p || g.ncols() != m {
        return -14;
    }
    if let Some(ref h) = hinvb {
        if h.nrows() != n || h.ncols() != m {
            return -18;
        }
    }
    if n == 0 {
        return 0;
    }

    let a = _a;
    // H = freq*I - A (complex)
    let mut h = DMatrix::<Complex64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            h[(i, j)] = if i == j {
                freq - Complex64::new(a[(i, j)], 0.0)
            } else {
                Complex64::new(-a[(i, j)], 0.0)
            };
        }
    }

    let lu = LU::new(h.clone());
    if !lu.is_invertible() {
        return 2;
    }
    let h_inv = match lu.try_inverse() {
        Some(inv) => inv,
        None => return 2,
    };

    if let Some(rc) = rcond {
        if matches!(baleig, Baleig::C | Baleig::A) {
            let norm_h = h.norm();
            let norm_inv = h_inv.norm();
            *rc = if norm_h > 0.0 && norm_inv > 0.0 {
                1.0 / (norm_h * norm_inv)
            } else {
                0.0
            };
            if *rc < 1e-15 {
                return 2;
            }
        }
    }

    let b_c: DMatrix<Complex64> = b.map(|x| Complex64::new(x, 0.0));
    let x = &h_inv * &b_c;

    if let Some(ref mut hout) = hinvb {
        hout.copy_from(&x);
    }

    let c_c: DMatrix<Complex64> = c.map(|x| Complex64::new(x, 0.0));
    let g_mat = &c_c * &x;
    g.copy_from(&g_mat);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb05ad_simple() {
        // A=0, B=[1;0], C=[1 0] -> G(freq) = 1/freq. At freq=1+0j, G=1.
        let n = 1usize;
        let m = 1usize;
        let p = 1usize;
        let mut a = DMatrix::from_element(1, 1, 0.0);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let freq = Complex64::new(1.0, 0.0);
        let mut g = DMatrix::from_element(1, 1, Complex64::new(0.0, 0.0));
        let mut rcond = 0.0;
        let info = tb05ad(
            Baleig::C,
            Inita::G,
            &mut a,
            &b,
            &c,
            freq,
            &mut g,
            Some(&mut rcond),
            None,
            None,
            None,
        );
        assert_eq!(info, 0);
        assert!((g[(0, 0)].re - 1.0).abs() < 1e-10);
        assert!((g[(0, 0)].im - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_tb05ad_slicot_example() {
        // SLICOT TB05AD example: N=3, M=1, P=2, freq=(0,0.5)
        let mut a = DMatrix::from_row_slice(
            3,
            3,
            &[1.0, 2.0, 0.0, 4.0, -1.0, 0.0, 0.0, 0.0, 1.0],
        );
        let b = DMatrix::from_row_slice(3, 1, &[1.0, 0.0, 1.0]);
        let c = DMatrix::from_row_slice(2, 3, &[1.0, 0.0, -1.0, 0.0, 0.0, 1.0]);
        let freq = Complex64::new(0.0, 0.5);
        let mut g = DMatrix::from_element(2, 1, Complex64::new(0.0, 0.0));
        let mut rcond = 0.0;
        let mut hinvb = DMatrix::from_element(3, 1, Complex64::new(0.0, 0.0));
        let info = tb05ad(
            Baleig::A,
            Inita::G,
            &mut a,
            &b,
            &c,
            freq,
            &mut g,
            Some(&mut rcond),
            Some(&mut hinvb),
            None,
            None,
        );
        assert_eq!(info, 0);
        // Expected G: (0.69, 0.35) and (-0.80, -0.40)
        assert!((g[(0, 0)].re - 0.69).abs() < 0.02);
        assert!((g[(0, 0)].im - 0.35).abs() < 0.02);
        assert!((g[(1, 0)].re - (-0.80)).abs() < 0.02);
        assert!((g[(1, 0)].im - (-0.40)).abs() < 0.02);
    }
}
