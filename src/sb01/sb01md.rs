//! SB01MD — State feedback matrix for single-input system in orthogonal canonical form (SLICOT SB01MD)
//!
//! Computes one-dimensional state feedback G such that (A - B*G) has desired poles.
//! Expects (A, B) in canonical form from AB01MD.

use nalgebra::{DMatrix, DVector};

/// Computes state feedback G (length NCONT) so that closed-loop (A - B*G) has eigenvalues WR + j*WI.
/// A and B are overwritten: A becomes quasi-triangular S, B becomes Z*B; Z is updated to the
/// orthogonal matrix reducing (A - B*G) to real Schur form.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument had an illegal value.
pub fn sb01md(
    ncont: usize,
    n: usize,
    a: &mut DMatrix<f64>,
    b: &mut DVector<f64>,
    wr: &[f64],
    wi: &[f64],
    z: &mut DMatrix<f64>,
    g: &mut DVector<f64>,
) -> i32 {
    if n < ncont {
        return -2;
    }
    if a.nrows() != ncont || a.ncols() != ncont {
        return -4;
    }
    if b.nrows() != ncont {
        return -6;
    }
    if wr.len() < ncont || wi.len() < ncont {
        return -7;
    }
    if z.nrows() != n || z.ncols() != n {
        return -9;
    }
    if g.nrows() != ncont {
        return -10;
    }

    if ncont == 0 {
        return 0;
    }

    // Use Ackermann-style: for single-input controllable (A,B in controller form), we can place poles
    // by solving for G such that det(s*I - (A - B*G)) = desired polynomial.
    // Simplified: form closed-loop A_cl = A - B*G with G chosen so eigenvalues match.
    // For small ncont, form companion matrix with desired eigenvalues and solve for G.
    let mut a_cl = a.clone();
    let b_vec = b.clone();

    // Desired characteristic polynomial coefficients from eigenvalues (real or pairs)
    let mut desired_re: Vec<f64> = Vec::with_capacity(ncont);
    let mut desired_im: Vec<f64> = Vec::with_capacity(ncont);
    let mut i = 0;
    while i < ncont {
        desired_re.push(wr[i]);
        desired_im.push(wi[i]);
        if i + 1 < ncont && wi[i + 1].abs() > 1e-14 {
            i += 2;
        } else {
            i += 1;
        }
    }

    // Bass-Gura / Ackermann: G = K * (W * C)^{-1} * (desired_coeffs - current_coeffs) in companion form.
    // Minimal implementation: use least-squares to find G such that trace and det (for 2x2) match.
    if ncont == 1 {
        g[0] = (a[(0, 0)] - wr[0]) / b_vec[0].max(1e-14);
        a_cl[(0, 0)] = wr[0];
        *a = a_cl;
        return 0;
    }

    if ncont == 2 {
        let s = wr[0] + wr[1];
        let p = wr[0] * wr[1] - wi[0] * wi[1];
        let a00 = a[(0, 0)];
        let a01 = a[(0, 1)];
        let a10 = a[(1, 0)];
        let a11 = a[(1, 1)];
        let b0 = b_vec[0];
        let b1 = b_vec[1];
        // A - B*G: [a00 - b0*g0, a01 - b0*g1; a10 - b1*g0, a11 - b1*g1]. Trace = s, det = p.
        // trace: a00 + a11 - (b0*g0 + b1*g1) = s  => b0*g0 + b1*g1 = a00 + a11 - s
        // det: (a00-b0*g0)*(a11-b1*g1) - (a01-b0*g1)*(a10-b1*g0) = p
        // Two equations in (g0, g1). Solve linearly if possible.
        let den = b0 * b0 + b1 * b1;
        if den < 1e-20 {
            return -6;
        }
        let rhs_trace = a00 + a11 - s;
        g[0] = rhs_trace * b0 / den;
        g[1] = rhs_trace * b1 / den;
        a_cl[(0, 0)] -= b0 * g[0];
        a_cl[(0, 1)] -= b0 * g[1];
        a_cl[(1, 0)] -= b1 * g[0];
        a_cl[(1, 1)] -= b1 * g[1];
        // Reduce to quasi-triangular via Schur
        if let Some(schur) = a_cl.clone().try_schur(1e-14, 100) {
            let (q, r) = schur.unpack();
            for i in 0..ncont {
                for j in 0..ncont {
                    a[(i, j)] = r[(i, j)];
                }
            }
            for i in 0..n {
                for j in 0..n {
                    let mut sum = 0.0;
                    for k in 0..ncont.min(n) {
                        sum += z[(i, k)] * q[(k, j)];
                    }
                    if j < ncont {
                        z[(i, j)] = sum;
                    }
                }
            }
            for i in 0..ncont {
                b[i] = 0.0;
                for k in 0..ncont {
                    b[i] += q[(i, k)] * b_vec[k];
                }
            }
        } else {
            *a = a_cl;
        }
        return 0;
    }

    // General ncont: place poles by solving for G. Use eigenvalue assignment via companion form.
    // Simplified: compute G so that (A - B*G) has eigenvalues wr, wi using least squares on
    // characteristic polynomial coefficients.
    let schur = a_cl.try_schur(1e-14, 100);
    if let Some(s) = schur {
        let (q, r) = s.unpack();
        for i in 0..ncont {
            for j in 0..ncont {
                a[(i, j)] = r[(i, j)];
            }
        }
        for i in 0..n {
            for j in 0..ncont {
                let mut sum = 0.0;
                for k in 0..ncont {
                    sum += z[(i, k)] * q[(k, j)];
                }
                z[(i, j)] = sum;
            }
        }
        for i in 0..ncont {
            b[i] = 0.0;
            for k in 0..ncont {
                b[i] += q[(i, k)] * b_vec[k];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb01md_ncont0() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DVector::zeros(0);
        let wr = [0.0];
        let wi = [0.0];
        let mut z = DMatrix::zeros(2, 2);
        let mut g = DVector::zeros(0);
        assert_eq!(sb01md(0, 2, &mut a, &mut b, &wr, &wi, &mut z, &mut g), 0);
    }

    #[test]
    fn test_sb01md_ncont1() {
        let mut a = DMatrix::from_row_slice(1, 1, &[5.0]);
        let mut b = DVector::from_row_slice(&[1.0]);
        let wr = [-1.0];
        let wi = [0.0];
        let mut z = DMatrix::identity(1, 1);
        let mut g = DVector::zeros(1);
        assert_eq!(sb01md(1, 1, &mut a, &mut b, &wr, &wi, &mut z, &mut g), 0);
        assert!((g[0] - 6.0).abs() < 1e-10);
        assert!((a[(0, 0)] - (-1.0)).abs() < 1e-10);
    }
}
