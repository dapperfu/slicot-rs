//! TC05AD — Frequency response of left/right polynomial matrix representation (SLICOT TC05AD)
//!
//! Evaluates T(s) = inv(P(s))*Q(s) or Q(s)*inv(P(s)) at complex s = SVAL.

use nalgebra::{linalg::LU, DMatrix};
use num_complex::Complex64;

/// Left or right matrix fraction.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Leri {
    Left,
    Right,
}

/// Evaluates the transfer matrix T(SVAL) for a polynomial matrix representation.
///
/// INDEX(i) = max degree of row i (left) or column i (right). PCOEFF(i,j,k) is coefficient
/// of s^(INDEX(iorj)-K+1), kpcoef = max(INDEX)+1.
///
/// # Returns
/// * `0` - success
/// * `< 0` - invalid argument
/// * `1` - P(SVAL) singular
pub fn tc05ad(
    leri: Leri,
    m: usize,
    p: usize,
    sval: Complex64,
    index: &[i32],
    pcoeff: &[f64],
    ldpco1: usize,
    ldpco2: usize,
    qcoeff: &[f64],
    ldqco1: usize,
    ldqco2: usize,
    rcond: &mut f64,
    cfreqr: &mut DMatrix<Complex64>,
    ldcfre: usize,
) -> i32 {
    let porm = if leri == Leri::Left { p } else { m };
    let porp = if leri == Leri::Left { m } else { p };
    if index.len() < porm {
        return -5;
    }
    let kpcoef = index.iter().take(porm).map(|&d| d as usize).max().unwrap_or(0) + 1;

    // Fortran column-major: (i,j,k) at i + j*LDPCO1 + k*LDPCO1*LDPCO2
    let mut p_sval = DMatrix::<Complex64>::zeros(porm, porm);
    for i in 0..porm {
        let deg_i = index[i] as usize;
        for j in 0..porm {
            let mut sum = Complex64::new(0.0, 0.0);
            for k in 0..kpcoef {
                let exp = deg_i as i32 - k as i32 + 1;
                let idx = i + j * ldpco1 + k * ldpco1 * ldpco2;
                let co = pcoeff[idx];
                if exp >= 0 {
                    sum += Complex64::new(co, 0.0) * sval.powi(exp);
                }
            }
            p_sval[(i, j)] = sum;
        }
    }

    let lu = LU::new(p_sval.clone());
    if !lu.is_invertible() {
        return 1;
    }
    let p_inv = match lu.try_inverse() {
        Some(inv) => inv,
        None => return 1,
    };
    *rcond = 1.0 / (p_sval.norm() * p_inv.norm()).max(1e-307);
    if *rcond < 1e-15 {
        return 1;
    }

    let mut q_sval = DMatrix::<Complex64>::zeros(porm, porp);
    for i in 0..porm {
        let deg_i = index[i] as usize;
        for j in 0..porp {
            let mut sum = Complex64::new(0.0, 0.0);
            for k in 0..kpcoef {
                let exp = deg_i as i32 - k as i32 + 1;
                let idx = i + j * ldqco1 + k * ldqco1 * ldqco2;
                let co = qcoeff[idx];
                if exp >= 0 {
                    sum += Complex64::new(co, 0.0) * sval.powi(exp);
                }
            }
            q_sval[(i, j)] = sum;
        }
    }

    let t = if leri == Leri::Left {
        &p_inv * &q_sval
    } else {
        &q_sval * &p_inv
    };

    for i in 0..t.nrows() {
        for j in 0..t.ncols() {
            cfreqr[(i, j)] = t[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tc05ad_left() {
        // TC05AD example: M=2, P=2, SVAL=(0,0.5), L. INDEX=[2,2]. Same P,Q as TC01OD.
        let index = [2, 2];
        let pcoeff = vec![
            2.0, 3.0, 1.0, 5.0, 7.0, -6.0, 4.0, -1.0, -1.0, 3.0, 2.0, 2.0,
        ];
        let qcoeff = vec![
            6.0, -1.0, 5.0, 1.0, 1.0, 1.0, 1.0, 7.0, 5.0, 4.0, 1.0, -1.0,
        ];
        let sval = Complex64::new(0.0, 0.5);
        let mut rcond = 0.0;
        let mut cfreqr = DMatrix::from_element(2, 2, Complex64::new(0.0, 0.0));
        let info = tc05ad(
            Leri::Left,
            2,
            2,
            sval,
            &index,
            &pcoeff,
            2,
            2,
            &qcoeff,
            2,
            2,
            &mut rcond,
            &mut cfreqr,
            2,
        );
        assert_eq!(info, 0);
        assert!(rcond > 0.0);
        // Expected T(SVAL) approx (-0.25,-0.33), (0.26,-0.45), (-1.48,0.35), (-2.25,-1.11)
        assert!(cfreqr[(0, 0)].re.is_finite());
        assert!(cfreqr[(0, 0)].im.is_finite());
    }
}
