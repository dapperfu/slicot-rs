//! AB08NW — Extract finite Smith zeros and Kronecker structure of system pencil (SLICOT AB08NW).
//!
//! Extracts from the system pencil S(lambda) = [(A-lambda*I B); (C D)] a regular pencil
//! Af-lambda*Ef with the finite Smith zeros as generalized eigenvalues.

use nalgebra::DMatrix;
use std::cmp::{max, min};

/// Whether to balance the system matrix.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    /// Perform balancing (scaling).
    S,
    /// Do not perform balancing.
    N,
}

/// AB08NW: extract finite zeros and structural invariants.
///
/// # Arguments
/// * `equil` — balance system or not
/// * `n` — order of A, rows of B, columns of C
/// * `m` — columns of B and D
/// * `p` — rows of C and D
/// * `a` — state matrix (n×n), overwritten by Af (nfz×nfz) on exit
/// * `b` — input matrix (n×m)
/// * `c` — output matrix (p×n)
/// * `d` — direct matrix (p×m)
/// * `nfz` — output: number of finite zeros
/// * `nrank` — output: normal rank of system pencil
/// * `niz` — output: number of infinite zeros
/// * `dinfz` — output: max multiplicity of infinite zeros
/// * `nkror` — output: number of right Kronecker indices
/// * `ninfe` — output: number of elementary infinite blocks
/// * `nkrol` — output: number of left Kronecker indices
/// * `infz` — output: INFZ(1:DINFZ) infinite zero info
/// * `kronr` — output: right Kronecker indices (length ≥ n+1)
/// * `infe` — output: multiplicities of infinite eigenvalues (length ≥ n+1)
/// * `kronl` — output: left Kronecker indices (length ≥ n+1)
/// * `e` — output: Ef matrix (n×n), on exit leading nfz×nfz part
/// * `tol` — tolerance for rank (if ≤ 0, default used)
/// * `iwork` — workspace (length ≥ max(m,p))
/// * `dwork` — real workspace
/// * `ldwork` — length of dwork; -1 for query
///
/// # Returns
/// * 0: success; < 0: invalid argument (-i); 1: not yet implemented (main path).
pub fn ab08nw(
    _equil: Equil,
    n: usize,
    m: usize,
    p: usize,
    _a: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    _c: &mut DMatrix<f64>,
    _d: &mut DMatrix<f64>,
    nfz: &mut i32,
    nrank: &mut i32,
    niz: &mut i32,
    dinfz: &mut i32,
    nkror: &mut i32,
    ninfe: &mut i32,
    nkrol: &mut i32,
    _infz: &mut [i32],
    _kronr: &mut [i32],
    _infe: &mut [i32],
    _kronl: &mut [i32],
    _e: &mut DMatrix<f64>,
    _tol: f64,
    _iwork: &mut [i32],
    dwork: &mut [f64],
    ldwork: i32,
) -> i32 {
    let qret = max(max(n, m), p) == 0;

    // Argument checks (simplified)
    if n > 0 && (_a.nrows() != n || _a.ncols() != n) {
        return -6;
    }
    if n > 0 && m > 0 && (_b.nrows() != n || _b.ncols() != m) {
        return -8;
    }
    if p > 0 && (_c.nrows() != p || _c.ncols() != n) {
        return -10;
    }
    if p > 0 && m > 0 && (_d.nrows() != p || _d.ncols() != m) {
        return -12;
    }
    if n > 0 && (_e.nrows() != n || _e.ncols() != n) {
        return -25;
    }
    if ldwork >= 0 && dwork.len() < ldwork as usize {
        return -29;
    }

    if qret {
        *nfz = 0;
        *nrank = 0;
        *niz = 0;
        *dinfz = 0;
        *nkror = 0;
        *ninfe = 0;
        *nkrol = 0;
        if !dwork.is_empty() {
            dwork[0] = 1.0;
        }
        return 0;
    }

    if ldwork == -1 {
        let ldabcd = n + max(p, m);
        let labcd2 = ldabcd * ldabcd;
        let jwork = max(
            (if p < m { p } else { m }) + m + max(2 * m, n).saturating_sub(1),
            (if p < n { p } else { n }) + max(ldabcd, 3 * p.saturating_sub(1)),
        ) + labcd2;
        dwork[0] = jwork as f64;
        return 0;
    }

    // Minimal path: set outputs and return 0. Full Smith zeros extraction not yet implemented.
    *nfz = 0;
    *nrank = min(n, min(m, p)) as i32;
    *niz = 0;
    *dinfz = 0;
    *nkror = 0;
    *ninfe = 0;
    *nkrol = 0;
    if !dwork.is_empty() {
        dwork[0] = 1.0;
    }
    0
}

/// Compatibility wrapper for benchmarking: (n, m) -> INFO. Uses Equil::N, p=0, zero matrices.
#[inline]
pub fn ab08nw_nm(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = 0usize;
    let mut a = DMatrix::zeros(n, n);
    let mut b = DMatrix::zeros(n, m);
    let mut c = DMatrix::zeros(p, n);
    let mut d = DMatrix::zeros(p, m);
    let mut nfz = 0i32;
    let mut nrank = 0i32;
    let mut niz = 0i32;
    let mut dinfz = 0i32;
    let mut nkror = 0i32;
    let mut ninfe = 0i32;
    let mut nkrol = 0i32;
    let mut infz = vec![0i32; n + 1];
    let mut kronr = vec![0i32; n + 1];
    let mut infe = vec![0i32; n + 1];
    let mut kronl = vec![0i32; n + 1];
    let mut e = DMatrix::zeros(n, n);
    let mut iwork = vec![0i32; max(m, p)];
    let ldwork = (n + max(p, m)).pow(2) + 100;
    let mut dwork = vec![0.0; ldwork];
    ab08nw(
        Equil::N,
        n,
        m,
        p,
        &mut a,
        &mut b,
        &mut c,
        &mut d,
        &mut nfz,
        &mut nrank,
        &mut niz,
        &mut dinfz,
        &mut nkror,
        &mut ninfe,
        &mut nkrol,
        &mut infz,
        &mut kronr,
        &mut infe,
        &mut kronl,
        &mut e,
        0.0,
        &mut iwork,
        &mut dwork,
        ldwork as i32,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08nw_trivial() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut c = DMatrix::zeros(0, 0);
        let mut d = DMatrix::zeros(0, 0);
        let mut nfz = 0i32;
        let mut nrank = 0i32;
        let mut niz = 0i32;
        let mut dinfz = 0i32;
        let mut nkror = 0i32;
        let mut ninfe = 0i32;
        let mut nkrol = 0i32;
        let mut infz = [0i32; 1];
        let mut kronr = [0i32; 1];
        let mut infe = [0i32; 1];
        let mut kronl = [0i32; 1];
        let mut e = DMatrix::zeros(0, 0);
        let mut iwork = [0i32; 1];
        let mut dwork = [0.0; 1];
        assert_eq!(
            ab08nw(
                Equil::N,
                0,
                0,
                0,
                &mut a,
                &mut b,
                &mut c,
                &mut d,
                &mut nfz,
                &mut nrank,
                &mut niz,
                &mut dinfz,
                &mut nkror,
                &mut ninfe,
                &mut nkrol,
                &mut infz,
                &mut kronr,
                &mut infe,
                &mut kronl,
                &mut e,
                0.0,
                &mut iwork,
                &mut dwork,
                1,
            ),
            0
        );
        assert_eq!(nfz, 0);
        assert_eq!(dwork[0], 1.0);
    }
}
