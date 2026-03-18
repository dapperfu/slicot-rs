//! AB05RD — Closed-loop system for mixed output and state feedback (SLICOT AB05RD).
//!
//! u = alpha*F*y + beta*K*x + G*v, z = H*y.
//! First applies output feedback via AB05SD, then Ac = A1 + beta*B1*K,
//! Bc = B1*G, Cc = H*(C1 + beta*D1*K), Dc = H*D1*G.

use nalgebra::DMatrix;

use super::ab05sd::ab05sd;

/// Mixed output and state feedback closed-loop.
/// FBTYPE: b'I' = F identity, b'O' = general F. JOBD: b'D' = D present, b'Z' = D zero.
/// Returns 0 on success; 1 if I - alpha*D*F singular; < 0 invalid argument.
pub fn ab05rd(
    fbtype: u8,
    jobd: u8,
    n: usize,
    m: usize,
    p: usize,
    mv: usize,
    pz: usize,
    alpha: f64,
    beta: f64,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    f: &DMatrix<f64>,
    k: &DMatrix<f64>,
    g: &DMatrix<f64>,
    h: &DMatrix<f64>,
    rcond: &mut f64,
    bc: &mut DMatrix<f64>,
    cc: &mut DMatrix<f64>,
    dc: &mut DMatrix<f64>,
) -> i32 {
    *rcond = 1.0;
    if n == 0 && m == 0 && p == 0 && mv == 0 && pz == 0 {
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
        return -11;
    }
    if b.nrows() != n || b.ncols() != m {
        return -13;
    }
    if c.nrows() != p || c.ncols() != n {
        return -15;
    }
    if ljobd && (d.nrows() != p || d.ncols() != m) {
        return -17;
    }
    if outpf && alpha != 0.0 && (f.nrows() != m || f.ncols() != p) {
        return -19;
    }
    if beta != 0.0 && (k.nrows() != m || k.ncols() != n) {
        return -21;
    }
    if g.nrows() != m || g.ncols() != mv {
        return -23;
    }
    if h.nrows() != pz || h.ncols() != p {
        return -25;
    }
    if bc.nrows() != n || bc.ncols() != mv {
        return -28;
    }
    if cc.nrows() != pz || cc.ncols() != n {
        return -30;
    }
    if ljobd && (dc.nrows() != pz || dc.ncols() != mv) {
        return -32;
    }
    if n == 0 || (m == 0 && p == 0) || (mv == 0 && pz == 0) {
        return 0;
    }
    let info = ab05sd(fbtype, jobd, n, m, p, alpha, a, b, c, d, f, rcond);
    if info != 0 {
        return info;
    }
    if beta != 0.0 && n > 0 {
        *a += beta * (&*b) * k;
        if ljobd {
            *c += beta * (&*d) * k;
        }
    }
    *bc = (&*b) * g;
    if n > 0 {
        *cc = h * (&*c);
    }
    if ljobd {
        *dc = h * (&*d) * g;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab05rd_trivial() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut c = DMatrix::zeros(0, 0);
        let mut d = DMatrix::zeros(0, 0);
        let f = DMatrix::zeros(0, 0);
        let k = DMatrix::zeros(0, 0);
        let g = DMatrix::zeros(0, 0);
        let h = DMatrix::zeros(0, 0);
        let mut rcond = 0.0;
        let mut bc = DMatrix::zeros(0, 0);
        let mut cc = DMatrix::zeros(0, 0);
        let mut dc = DMatrix::zeros(0, 0);
        assert_eq!(
            ab05rd(
                b'Z', b'Z', 0, 0, 0, 0, 0, 0.0, 0.0,
                &mut a, &mut b, &mut c, &mut d,
                &f, &k, &g, &h,
                &mut rcond, &mut bc, &mut cc, &mut dc
            ),
            0
        );
    }

    #[test]
    fn test_ab05rd_alpha_beta_zero() {
        let mut a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let f = DMatrix::from_row_slice(1, 1, &[1.0]);
        let k = DMatrix::from_row_slice(1, 1, &[0.0]);
        let g = DMatrix::from_row_slice(1, 1, &[1.0]);
        let h = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut rcond = 0.0;
        let mut bc = DMatrix::zeros(1, 1);
        let mut cc = DMatrix::zeros(1, 1);
        let mut dc = DMatrix::zeros(1, 1);
        let info = ab05rd(
            b'I', b'D', 1, 1, 1, 1, 1, 0.0, 0.0,
            &mut a, &mut b, &mut c, &mut d,
            &f, &k, &g, &h,
            &mut rcond, &mut bc, &mut cc, &mut dc,
        );
        assert_eq!(info, 0);
        assert!((bc[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((cc[(0, 0)] - 1.0).abs() < 1e-10);
    }
}
