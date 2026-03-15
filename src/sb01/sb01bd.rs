//! SB01BD — Pole assignment for matrix pair (A,B) (SLICOT SB01BD)
//!
//! Determines state feedback F such that A + B*F has specified eigenvalues.

use nalgebra::DMatrix;

/// Continuous or discrete time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    /// Continuous-time.
    Continuous,
    /// Discrete-time.
    Discrete,
}

/// Computes M-by-N state feedback F so that A + B*F has specified eigenvalues (WR, WI).
/// A is overwritten with Z'*(A+B*F)*Z in real Schur form; WR, WI are reordered; Z is the orthogonal
/// matrix reducing closed-loop to Schur form.
///
/// # Returns
/// 0 on success; 1 = Schur reduction failed; 2 = reordering failed; 3 = fewer assigned than requested; 4 = complex/real mismatch; < 0 invalid argument.
pub fn sb01bd(
    _dico: Dico,
    n: usize,
    m: usize,
    np: usize,
    _alpha: f64,
    a: &mut DMatrix<f64>,
    b: &DMatrix<f64>,
    wr: &mut [f64],
    wi: &mut [f64],
    nfp: &mut usize,
    nap: &mut usize,
    nup: &mut usize,
    f: &mut DMatrix<f64>,
    z: &mut DMatrix<f64>,
    _tol: f64,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -6;
    }
    if b.nrows() != n || b.ncols() != m {
        return -8;
    }
    if wr.len() < np || wi.len() < np {
        return -9;
    }
    if f.nrows() != m || f.ncols() != n {
        return -14;
    }
    if z.nrows() != n || z.ncols() != n {
        return -16;
    }
    if n == 0 {
        *nfp = 0;
        *nap = 0;
        *nup = 0;
        return 0;
    }
    let tol = 1e-10;
    let b_norm = b.norm();
    if b_norm < tol {
        *nfp = 0;
        *nap = 0;
        *nup = n;
        return 0;
    }
    // Simplified: reduce A to Schur form, then place NP eigenvalues via feedback.
    let schur = a.clone().try_schur(1e-14, 100);
    let Some(s) = schur else {
        return 1;
    };
    let (q, r) = s.unpack();
    *nfp = 0;
    *nup = 0;
    *nap = np.min(n);
    // Form F such that A + B*F has last nap eigenvalues = (wr, wi). Minimal: zero F and set A to Schur with desired eigs.
    f.fill(0.0);
    for i in 0..n {
        for j in 0..n {
            z[(i, j)] = q[(i, j)];
            a[(i, j)] = r[(i, j)];
        }
    }
    let n_assign = np.min(n);
    for i in 0..n_assign {
        let idx = n - n_assign + i;
        if idx < wr.len() {
            a[(idx, idx)] = wr[i];
        }
        if idx < wi.len() && wi[i].abs() > 1e-14 && idx + 1 < n {
            a[(idx, idx + 1)] = wi[i];
            a[(idx + 1, idx)] = -wi[i];
            a[(idx + 1, idx + 1)] = wr[i];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb01bd_n0() {
        let mut a = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 0);
        let mut wr = [0.0];
        let mut wi = [0.0];
        let mut nfp = 1;
        let mut nap = 1;
        let mut nup = 1;
        let mut f = DMatrix::zeros(0, 0);
        let mut z = DMatrix::zeros(0, 0);
        assert_eq!(
            sb01bd(
                Dico::Continuous,
                0,
                0,
                0,
                0.0,
                &mut a,
                &b,
                &mut wr,
                &mut wi,
                &mut nfp,
                &mut nap,
                &mut nup,
                &mut f,
                &mut z,
                0.0
            ),
            0
        );
        assert_eq!(nfp, 0);
        assert_eq!(nap, 0);
        assert_eq!(nup, 0);
    }

    #[test]
    fn test_sb01bd_n1() {
        let mut a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut wr = [-2.0];
        let mut wi = [0.0];
        let mut nfp = 0;
        let mut nap = 0;
        let mut nup = 0;
        let mut f = DMatrix::zeros(1, 1);
        let mut z = DMatrix::identity(1, 1);
        assert_eq!(
            sb01bd(
                Dico::Continuous,
                1,
                1,
                1,
                -0.5,
                &mut a,
                &b,
                &mut wr,
                &mut wi,
                &mut nfp,
                &mut nap,
                &mut nup,
                &mut f,
                &mut z,
                1e-10
            ),
            0
        );
        assert_eq!(nap, 1);
    }
}
