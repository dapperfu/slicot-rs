//! TB04AD — Transfer matrix of a given state-space representation (A,B,C,D) (SLICOT TB04AD)
//!
//! Finds T(s) expressed as row or column polynomial vectors over monic denominator polynomials.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum RowCol {
    /// T(s) as rows over common denominators.
    R,
    /// T(s) as columns over common denominators.
    C,
}

/// Faddeev-LeVerrier: characteristic polynomial phi(s) = s^n + p[1]*s^{n-1} + ... + p[n].
/// Returns p[0]=1, p[1..=n] in `coeff`, coeff[k] = coefficient of s^{n-k}.
fn char_poly(n: usize, a: &DMatrix<f64>, coeff: &mut [f64]) {
    if n == 0 {
        return;
    }
    coeff[0] = 1.0;
    let mut m = DMatrix::identity(n, n);
    for k in 1..=n {
        m = a * m;
        let tr = m.trace();
        coeff[k] = -tr / (k as f64);
        if k < n {
            for i in 0..n {
                m[(i, i)] += coeff[k];
            }
        }
    }
}

/// Computes transfer matrix T(s) of (A,B,C,D) as row or column polynomial vectors.
///
/// On exit A, B, C contain the transformed (block Hessenberg) representation; NR is its order.
/// INDEX(i) = degree of denominator for row/column i; DCOEFF and UCOEFF hold coefficients.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument had an illegal value.
pub fn tb04ad(
    rowcol: RowCol,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &DMatrix<f64>,
    nr: &mut usize,
    index: &mut [i32],
    dcoeff: &mut [f64],
    lddcoe: usize,
    ucoeff: &mut [f64],
    lduco1: usize,
    lduco2: usize,
    tol1: f64,
    tol2: f64,
    iwork: &mut [i32],
    _dwork: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if a.ncols() != n {
        return -6;
    }
    if b.nrows() != n || c.ncols() != n || d.nrows() != p || d.ncols() != m {
        return -8;
    }
    let porm = if rowcol == RowCol::R { p } else { m };
    let porp = if rowcol == RowCol::R { m } else { p };
    if index.len() < porm {
        return -14;
    }
    if lddcoe < porm {
        return -16;
    }
    if lduco1 < porm || lduco2 < porp {
        return -18;
    }
    if iwork.len() < n + m.max(p) {
        return -20;
    }

    if n == 0 || (rowcol == RowCol::R && m == 0) || (rowcol == RowCol::C && p == 0) {
        *nr = 0;
        for i in 0..porm {
            index[i] = 0;
        }
        return 0;
    }

    let mut ncont = 0_usize;
    let mut indcon = 0_usize;
    let mut nblk = vec![0i32; n];
    let mut tau = vec![0.0; n];

    let info = crate::tb01::tb01ud::tb01ud(
        crate::tb01::tb01ud::JobZ::No,
        a,
        b,
        c,
        &mut ncont,
        &mut indcon,
        &mut nblk,
        None,
        &mut tau,
        tol2,
    );
    if info != 0 {
        return info;
    }

    *nr = ncont.max(1).min(n);
    for (i, &blk) in nblk.iter().take(n).enumerate() {
        iwork[i] = blk;
    }
    for i in porm..iwork.len() {
        iwork[i] = 0;
    }

    // Common denominator representation: denominator = char poly of A(1:NR,1:NR)
    let nr_ = *nr;
    let mut phi = vec![0.0; nr_ + 1];
    char_poly(nr_, &a.view((0, 0), (nr_, nr_)).into_owned(), &mut phi);

    for i in 0..porm {
        index[i] = nr_ as i32;
    }
    let kdcoef = nr_ + 1;

    // DCOEFF: for each row (or column) the same denominator in this implementation
    // DCOEFF(I,K) = coefficient of s^(INDEX(I)-K+1), K=1..kdcoef => s^NR, s^{NR-1}, ..., s^0
    for i in 0..porm {
        for k in 0..kdcoef {
            let idx = i + k * lddcoe;
            if idx < dcoeff.len() {
                dcoeff[idx] = if k <= nr_ { phi[k] } else { 0.0 };
            }
        }
    }

    // H_k for adj(sI-A): H_0 = I, H_k = A*H_{k-1} + phi[k]*I
    let a_nr = a.view((0, 0), (nr_, nr_)).into_owned();
    let b_nr = b.view((0, 0), (nr_, m)).into_owned();
    let c_nr = c.view((0, 0), (p, nr_)).into_owned();

    for i in 0..p {
        for j in 0..m {
            // N_ij(s) = C_i*adj(sI-A)*B_j + D_ij*phi(s). adj = sum_{k=0}^{NR-1} H_k s^{NR-1-k}, H_0=I, H_k = A*H_{k-1}+phi[k]*I.
            let mut num_coeff = vec![0.0; kdcoef];
            num_coeff[0] = d[(i, j)] * phi[0]; // s^NR
            let mut h = DMatrix::identity(nr_, nr_);
            for k in 1..=nr_ {
                let cib = c_nr.row(i) * &h * b_nr.column(j);
                num_coeff[k] = cib[(0, 0)] + d[(i, j)] * phi[k];
                if k < nr_ {
                    h = &a_nr * &h + DMatrix::from_fn(nr_, nr_, |ii, jj| if ii == jj { phi[k] } else { 0.0 });
                }
            }
            for k in 0..kdcoef {
                let idx = i + j * lduco1 + k * lduco1 * lduco2;
                if idx < ucoeff.len() {
                    ucoeff[idx] = num_coeff[k];
                }
            }
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb04ad_smoke() {
        // N=3, M=2, P=2 from SLICOT example
        let mut a = DMatrix::from_row_slice(3, 3, &[-1.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 0.0, -3.0]);
        let mut b = DMatrix::from_row_slice(3, 2, &[0.0, 1.0, -1.0, 1.0, 1.0, 0.0]);
        let mut c = DMatrix::from_row_slice(2, 3, &[0.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        let d = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut nr = 0;
        let mut index = [0i32; 2];
        let mut dcoeff = vec![0.0; 2 * 4];
        let mut ucoeff = vec![0.0; 2 * 2 * 4];
        let mut iwork = vec![0i32; 3 + 2];
        let mut dwork = vec![0.0; 100];
        let info = tb04ad(
            RowCol::R,
            &mut a,
            &mut b,
            &mut c,
            &d,
            &mut nr,
            &mut index,
            &mut dcoeff,
            2,
            &mut ucoeff,
            2,
            2,
            0.0,
            0.0,
            &mut iwork,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert!(nr >= 1 && nr <= 3);
        assert!(index[0] >= 0 && index[1] >= 0);
    }

    #[test]
    fn test_tb04ad_zero_state() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 1);
        let mut c = DMatrix::zeros(1, 0);
        let d = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut nr = 99;
        let mut index = [99i32];
        let mut dcoeff = vec![0.0; 2];
        let mut ucoeff = vec![0.0; 2];
        let mut iwork = vec![0i32; 1];
        let mut dwork = vec![0.0; 1];
        let info = tb04ad(
            RowCol::R,
            &mut a,
            &mut b,
            &mut c,
            &d,
            &mut nr,
            &mut index,
            &mut dcoeff,
            1,
            &mut ucoeff,
            1,
            1,
            0.0,
            0.0,
            &mut iwork,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert_eq!(nr, 0);
        assert_eq!(index[0], 0);
    }
}
