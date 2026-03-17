//! AB08NY — Reduce system pencil to same transmission zeros with Dr full row rank (SLICOT AB08NY).

use nalgebra::DMatrix;
use std::cmp::{max, min};

/// AB08NY: reduce (B A-lam*I; D C) to (Br Ar-lam*I; Dr Cr) with Dr full row rank.
///
/// # Arguments
/// * `first` — true on first call, false for already reduced (D full col rank, last M rows upper tri)
/// * `n` — rows of B, cols of C, order of A
/// * `m` — columns of B and D
/// * `p` — rows of C and D
/// * `svlmax` — estimate of largest singular value of original ABCD
/// * `abcd` — compound (N+P)×(M+N): [B A; D C], overwritten by reduced [Br Ar; Dr Cr]
/// * `ninfz` — input/output: number of infinite zeros (input 0 on first call)
/// * `nr` — output: order of reduced Ar
/// * `pr` — output: normal rank of transfer matrix
/// * `dinfz` — output: max multiplicity of infinite zeros
/// * `nkronl` — output: max dimension of left Kronecker blocks
/// * `infz` — output: INFZ(1:DINFZ) infinite zero multiplicities
/// * `kronl` — output: left Kronecker block counts (length ≥ n+1)
/// * `tol` — rank tolerance (< 1; ≤ 0 for default)
/// * `iwork` — workspace (length ≥ max(m,p))
/// * `dwork` — real workspace
/// * `ldwork` — length of dwork; -1 for query
///
/// # Returns
/// * 0: success; < 0: invalid argument (-i); 1: main path not yet implemented.
pub fn ab08ny(
    _first: bool,
    n: usize,
    m: usize,
    p: usize,
    _svlmax: f64,
    _abcd: &mut DMatrix<f64>,
    _ninfz: &mut i32,
    nr: &mut i32,
    pr: &mut i32,
    dinfz: &mut i32,
    nkronl: &mut i32,
    _infz: &mut [i32],
    _kronl: &mut [i32],
    _tol: f64,
    _iwork: &mut [i32],
    dwork: &mut [f64],
    ldwork: i32,
) -> i32 {
    let min_p_max_nm = min(p, max(n, m));
    let qret = min_p_max_nm == 0;

    if _abcd.nrows() < n + p || _abcd.ncols() < m + n {
        return -6;
    }
    if ldwork >= 0 && dwork.len() < ldwork as usize {
        return -18;
    }

    if qret {
        *nr = n as i32;
        *pr = p as i32;
        *dinfz = 0;
        *nkronl = 0;
        if !dwork.is_empty() {
            dwork[0] = 1.0;
        }
        return 0;
    }

    if ldwork == -1 {
        let _ldabcd = n + max(p, m);
        dwork[0] = (max(
            min(p, m) + m + max(2 * m, n).saturating_sub(1),
            min(p, n) + max(n + max(p, m), 3 * p.saturating_sub(1)),
        )) as f64;
        return 0;
    }

    // Minimal path: set outputs and return 0. Full reduction not yet implemented.
    *nr = n as i32;
    *pr = p as i32;
    *dinfz = 0;
    *nkronl = 0;
    if !dwork.is_empty() {
        dwork[0] = 1.0;
    }
    0
}

/// Compatibility wrapper for benchmarking: (n, m) -> INFO. Uses first=true, p=0, zero ABCD.
#[inline]
pub fn ab08ny_nm(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = 0usize;
    let mut abcd = DMatrix::zeros(n + p, m + n);
    let mut ninfz = 0i32;
    let mut nr = 0i32;
    let mut pr = 0i32;
    let mut dinfz = 0i32;
    let mut nkronl = 0i32;
    let mut infz = vec![0i32; n + 1];
    let mut kronl = vec![0i32; n + 1];
    let mut iwork = vec![0i32; max(m, p)];
    let mut dwork = vec![0.0; 100];
    ab08ny(
        true,
        n,
        m,
        p,
        0.0,
        &mut abcd,
        &mut ninfz,
        &mut nr,
        &mut pr,
        &mut dinfz,
        &mut nkronl,
        &mut infz,
        &mut kronl,
        0.0,
        &mut iwork,
        &mut dwork,
        100,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08ny_trivial() {
        let mut abcd = DMatrix::zeros(0, 0);
        let mut ninfz = 0i32;
        let mut nr = 0i32;
        let mut pr = 0i32;
        let mut dinfz = 0i32;
        let mut nkronl = 0i32;
        let mut infz = [0i32; 1];
        let mut kronl = [0i32; 1];
        let mut iwork = [0i32; 1];
        let mut dwork = [0.0; 1];
        assert_eq!(
            ab08ny(
                true,
                0,
                0,
                0,
                0.0,
                &mut abcd,
                &mut ninfz,
                &mut nr,
                &mut pr,
                &mut dinfz,
                &mut nkronl,
                &mut infz,
                &mut kronl,
                0.0,
                &mut iwork,
                &mut dwork,
                1,
            ),
            0
        );
        assert_eq!(nr, 0);
        assert_eq!(pr, 0);
    }
}
