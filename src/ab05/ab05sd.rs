//! AB05SD — Closed-loop system for output feedback u = alpha*F*y + v (SLICOT AB05SD).
//!
//! Computes (Ac, Bc, Cc, Dc) with E = (I - alpha*D*F)^{-1},
//! Ac = A + alpha*B*F*E*C, Bc = B + alpha*B*F*E*D, Cc = E*C, Dc = E*D.

use nalgebra::DMatrix;

/// Output feedback closed-loop. FBTYPE: b'I' = F identity, b'O' = general F.
/// JOBD: b'D' = D present, b'Z' = D zero.
/// Returns 0 on success; 1 if I - alpha*D*F singular; < 0 invalid argument.
pub fn ab05sd(
    fbtype: u8,
    jobd: u8,
    n: usize,
    m: usize,
    p: usize,
    alpha: f64,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    f: &DMatrix<f64>,
    rcond: &mut f64,
) -> i32 {
    *rcond = 1.0;
    if n == 0 && (m == 0 || p == 0) {
        return 0;
    }
    let unitf = fbtype == b'I' || fbtype == b'i';
    let outpf = fbtype == b'O' || fbtype == b'o';
    let ljobd = jobd == b'D' || jobd == b'd';
    if !unitf && !outpf {
        return -1;
    }
    if !ljobd && jobd != b'Z' && jobd != b'z' {
        return -2;
    }
    if unitf && p != m {
        return -5;
    }
    if a.nrows() != n || a.ncols() != n {
        return -7;
    }
    if b.nrows() != n || b.ncols() != m {
        return -9;
    }
    if c.nrows() != p || c.ncols() != n {
        return -11;
    }
    if ljobd && (d.nrows() != p || d.ncols() != m) {
        return -13;
    }
    if outpf && alpha != 0.0 && (f.nrows() != m || f.ncols() != p) {
        return -16;
    }
    if n == 0 || m == 0 || p == 0 || alpha == 0.0 {
        return 0;
    }
    if ljobd {
        let d_copy = d.clone();
        let mut e = if unitf {
            -alpha * d_copy
        } else {
            -alpha * &d_copy * f
        };
        for i in 0..p {
            e[(i, i)] += 1.0;
        }
        let norm_e = e.norm();
        let e_inv = match e.try_inverse() {
            Some(x) => x,
            None => {
                *rcond = 0.0;
                return 1;
            },
        };
        let norm_einv = e_inv.norm();
        *rcond = if norm_e > 0.0 && norm_einv > 0.0 {
            1.0 / (norm_e * norm_einv)
        } else {
            0.0
        };
        if *rcond < f64::EPSILON {
            return 1;
        }
        let cc = &e_inv * c.clone();
        let dc = &e_inv * d.clone();
        *c = cc;
        *d = dc;
    }
    if n == 0 {
        return 0;
    }
    if unitf {
        let delta_a = alpha * (&*b) * c.clone();
        *a += delta_a;
        if ljobd {
            let delta_b = alpha * (&*b) * d.clone();
            *b += delta_b;
        }
    } else {
        let bf = (&*b) * f;
        let delta_a = alpha * &bf * c.clone();
        *a += delta_a;
        if ljobd {
            let delta_b = alpha * &bf * d.clone();
            *b += delta_b;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab05sd_trivial() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut c = DMatrix::zeros(0, 0);
        let mut d = DMatrix::zeros(0, 0);
        let f = DMatrix::zeros(0, 0);
        let mut rcond = 0.0;
        assert_eq!(
            ab05sd(b'Z', b'Z', 0, 0, 0, 0.0, &mut a, &mut b, &mut c, &mut d, &f, &mut rcond),
            0
        );
    }

    #[test]
    fn test_ab05sd_alpha_zero() {
        let mut a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let f = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut rcond = 0.0;
        let info = ab05sd(b'I', b'D', 1, 1, 1, 0.0, &mut a, &mut b, &mut c, &mut d, &f, &mut rcond);
        assert_eq!(info, 0);
        assert!((a[(0, 0)] - 1.0).abs() < 1e-10);
    }
}
