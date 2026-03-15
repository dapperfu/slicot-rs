//! AB13DX — Max singular value of G(lambda) at a given frequency omega (SLICOT).
//!
//! For standard system (E=I): lambda = j*omega (continuous) or exp(j*omega) (discrete).
//! Forms M = lambda*I - A (complex), solves M*X = B, then G = C*X + D and returns
//! the maximum singular value of G.

use nalgebra::{linalg::LU, DMatrix};
use num_complex::Complex;

/// Continuous ('C') or discrete ('D') time.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    /// Continuous: lambda = j*omega.
    Continuous,
    /// Discrete: lambda = exp(j*omega).
    Discrete,
}

/// E matrix: identity only supported.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobe {
    /// E = I (identity).
    Identity,
}

/// D matrix: zero or used.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Jobd {
    /// D = 0 (zero).
    Zero,
    /// D is given.
    Given,
}

/// Returns the maximum singular value of G(lambda) at frequency omega.
///
/// For E=I: forms M = lambda*I - A (complex), solves M*X = B, G = C*X + D,
/// then returns sigma_max(G). Implements JOBE='I', DICO 'C'/'D', JOBD 'Z'/'D'.
///
/// # Returns
/// info = 0 on success; &lt; 0 invalid argument;
/// &gt; 0 if M = lambda*I - A is singular. On success, `sigma_max` is the max singular value; `fpeak` is set to omega.
pub fn ab13dx(
    dico: Dico,
    _jobe: Jobe,
    jobd: Jobd,
    omega: f64,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    fpeak: &mut f64,
    sigma_max: &mut f64,
) -> i32 {
    *fpeak = omega;
    *sigma_max = 0.0;
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n {
        return -8;
    }
    if b.nrows() != n || c.ncols() != n || d.nrows() != p || d.ncols() != m {
        return -9;
    }
    if n == 0 {
        return 0;
    }

    let lambda: Complex<f64> = match dico {
        Dico::Continuous => Complex::new(0.0, omega),
        Dico::Discrete => Complex::new(omega.cos(), omega.sin()),
    };

    // M = lambda*I - A (complex n×n)
    let mut m = DMatrix::<Complex<f64>>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let v = if i == j {
                lambda - Complex::new(a[(i, j)], 0.0)
            } else {
                Complex::new(-a[(i, j)], 0.0)
            };
            m[(i, j)] = v;
        }
    }

    let lu = LU::new(m.clone());
    if !lu.is_invertible() {
        return 1;
    }
    let m_inv = match lu.try_inverse() {
        Some(inv) => inv,
        None => return 1,
    };

    // B as complex n×m
    let b_c: DMatrix<Complex<f64>> = b.map(|x| Complex::new(x, 0.0));
    // X = M^{-1} * B
    let x = &m_inv * &b_c;

    // G = C*X + D (p×m complex). C real p×n, X complex n×m, D real p×m
    let c_c: DMatrix<Complex<f64>> = c.map(|x| Complex::new(x, 0.0));
    let g = &c_c * &x;
    let g = if matches!(jobd, Jobd::Given) {
        let d_c: DMatrix<Complex<f64>> = d.map(|x| Complex::new(x, 0.0));
        g + d_c
    } else {
        g
    };

    // Max singular value of G (p×m complex); singular values are real
    let svd = g.svd(true, true);
    let sigmas = &svd.singular_values;
    *sigma_max = sigmas.iter().cloned().fold(0.0_f64, f64::max);
    *fpeak = omega;
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab13dx_n0_returns_zero() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let b = DMatrix::<f64>::zeros(0, 0);
        let c = DMatrix::<f64>::zeros(0, 0);
        let d = DMatrix::<f64>::zeros(0, 0);
        let mut fpeak = -1.0;
        let mut sigma_max = -1.0;
        assert_eq!(
            ab13dx(
                Dico::Continuous,
                Jobe::Identity,
                Jobd::Zero,
                1.0,
                &a,
                &b,
                &c,
                &d,
                &mut fpeak,
                &mut sigma_max,
            ),
            0
        );
        assert_eq!(sigma_max, 0.0);
        assert_eq!(fpeak, 1.0);
    }

    #[test]
    fn test_ab13dx_1x1_known_sigma() {
        // G(s) = 1/(s+1) => at omega=0, G(0)=1, sigma_max = 1. A=-1, B=1, C=1, D=0.
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut fpeak = 0.0;
        let mut sigma_max = 0.0;
        let info = ab13dx(
            Dico::Continuous,
            Jobe::Identity,
            Jobd::Zero,
            0.0,
            &a,
            &b,
            &c,
            &d,
            &mut fpeak,
            &mut sigma_max,
        );
        assert_eq!(info, 0);
        assert!((sigma_max - 1.0).abs() < 1e-10);
        assert_eq!(fpeak, 0.0);
    }
}
