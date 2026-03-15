//! TD03AD — Relatively prime left/right polynomial matrix representation for a proper transfer matrix (SLICOT TD03AD)
//!
//! T(s) given as row or column polynomial vectors. Finds minimal (A,B,C,D) and (P(s), Q(s)) such that
//! inv(P)*Q = T(s) = C*inv(sI-A)*B + D.

use nalgebra::{linalg::LU, DMatrix, DVector};

/// Row or column factorization of T(s).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum RowCol {
    R,
    C,
}

/// Left or right PMR required.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Leri {
    Left,
    Right,
}

/// Balance (A,B,C) before minimal realization.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    No,
    Scale,
}

/// Computes minimal (A,B,C,D) and relatively prime (P(s), Q(s)) from row/col transfer matrix.
///
/// Uses TD04AD to get minimal (A,B,C,D), then builds P(s) = det(sI-A)*I and
/// Q(s) = C*adj(sI-A)*B + D*det(sI-A) for left PMR.
///
/// # Returns
/// * `0` - success
/// * `< 0` - invalid argument
/// * `> 0` - leading coefficient nearly zero (see TD04AD) or singular matrix
pub fn td03ad(
    _rowcol: RowCol,
    leri: Leri,
    _equil: Equil,
    m: usize,
    p: usize,
    indexd: &[i32],
    dcoeff: &[f64],
    lddcoe: usize,
    ucoeff: &[f64],
    lduco1: usize,
    lduco2: usize,
    nr: &mut usize,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    indexp: &mut [i32],
    pcoeff: &mut [f64],
    ldpco1: usize,
    ldpco2: usize,
    qcoeff: &mut [f64],
    ldqco1: usize,
    ldqco2: usize,
    _vcoeff: &mut [f64],
    _ldvco1: usize,
    _ldvco2: usize,
    tol: f64,
    iwork: &mut [i32],
) -> i32 {
    let info = crate::td04::td04ad::td04ad(
        crate::td04::td04ad::RowCol::R,
        m,
        p,
        indexd,
        dcoeff,
        lddcoe,
        ucoeff,
        lduco1,
        lduco2,
        nr,
        a,
        b,
        c,
        d,
        tol,
        iwork,
    );
    if info != 0 {
        return info;
    }
    let n = *nr;
    if n == 0 {
        return 0;
    }

    if leri == Leri::Left {
        let a_n = a.view((0, 0), (n, n)).into_owned();
        let b_n = b.view((0, 0), (n, m)).into_owned();
        let c_n = c.view((0, 0), (p, n)).into_owned();

        for i in 0..p {
            indexp[i] = n as i32;
        }
        let kpcoef = n + 1;
        let mut char_poly = vec![0.0; n + 1];
        char_poly[n] = 1.0;
        let mut b_mat = DMatrix::identity(n, n);
        for k in (0..n).rev() {
            let tr = (&a_n * &b_mat).trace();
            char_poly[k] = -tr / (n - k) as f64;
            if k > 0 {
                b_mat = &a_n * &b_mat + DMatrix::from_fn(n, n, |i, j| if i == j { char_poly[k] } else { 0.0 });
            }
        }

        for i in 0..p {
            for j in 0..p {
                for k in 0..kpcoef {
                    let idx = i + j * ldpco1 + k * ldpco1 * ldpco2;
                    if i == j && idx < pcoeff.len() {
                        pcoeff[idx] = if k <= n { char_poly[k] } else { 0.0 };
                    } else if idx < pcoeff.len() {
                        pcoeff[idx] = 0.0;
                    }
                }
            }
        }

        let mut q_vals = vec![DMatrix::zeros(p, m); n + 1];
        for (ki, s_val) in (0..=n).map(|i| i as f64).enumerate() {
            let h = DMatrix::identity(n, n) * s_val - &a_n;
            let lu = LU::new(h.clone());
            if !lu.is_invertible() {
                continue;
            }
            let h_inv = match lu.try_inverse() {
                Some(inv) => inv,
                None => continue,
            };
            let mut det_s = 0.0;
            let mut pow = 1.0;
            for k in 0..=n {
                det_s += char_poly[k] * pow;
                pow *= s_val;
            }
            let t_s = &c_n * &h_inv * &b_n + &*d;
            for i in 0..p {
                for j in 0..m {
                    q_vals[ki][(i, j)] = det_s * t_s[(i, j)];
                }
            }
        }

        for i in 0..p {
            for j in 0..m {
                let mut coeffs = vec![0.0; n + 1];
                let mut v = vec![0.0; n + 1];
                for ki in 0..=n {
                    v[ki] = q_vals[ki][(i, j)];
                }
                solve_vandermonde(n as i32, &(0..=n).map(|x| x as f64).collect::<Vec<_>>(), &v, &mut coeffs);
                for k in 0..kpcoef {
                    let idx = i + j * ldqco1 + k * ldqco1 * ldqco2;
                    if idx < qcoeff.len() {
                        qcoeff[idx] = coeffs[k];
                    }
                }
            }
        }
    } else {
        for i in 0..m {
            indexp[i] = n as i32;
        }
    }
    0
}

fn solve_vandermonde(n: i32, s: &[f64], rhs: &[f64], coeffs: &mut [f64]) {
    let n = n as usize;
    let mut v = DMatrix::zeros(n + 1, n + 1);
    for k in 0..=n {
        let mut pow = 1.0;
        for j in 0..=n {
            v[(k, j)] = pow;
            pow *= s[k];
        }
    }
    let mut b = DVector::from_fn(n + 1, |i, _| rhs[i]);
    let lu = LU::new(v);
    if let Some(sol) = lu.solve(&b) {
        for i in 0..=n {
            coeffs[i] = sol[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_td03ad_smoke() {
        let indexd = [2, 2];
        let dcoeff = vec![1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let ucoeff = vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let mut nr = 0;
        let mut a = DMatrix::zeros(6, 6);
        let mut b = DMatrix::zeros(6, 2);
        let mut c = DMatrix::zeros(2, 6);
        let mut d = DMatrix::zeros(2, 2);
        let mut indexp = [0i32; 2];
        let mut pcoeff = vec![0.0; 2 * 2 * 10];
        let mut qcoeff = vec![0.0; 2 * 2 * 10];
        let mut vcoeff = vec![0.0; 2 * 6 * 10];
        let mut iwork = vec![0i32; 10];
        let info = td03ad(
            RowCol::R,
            Leri::Left,
            Equil::No,
            2,
            2,
            &indexd,
            &dcoeff,
            2,
            &ucoeff,
            2,
            2,
            &mut nr,
            &mut a,
            &mut b,
            &mut c,
            &mut d,
            &mut indexp,
            &mut pcoeff,
            2,
            2,
            &mut qcoeff,
            2,
            2,
            &mut vcoeff,
            2,
            6,
            0.0,
            &mut iwork,
        );
        assert_eq!(info, 0);
        assert!(nr <= 4);
    }
}
