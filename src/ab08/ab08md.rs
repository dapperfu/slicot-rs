//! AB08MD — Normal rank of the transfer-function matrix (SLICOT AB08MD)
//
// For (A,B,C,D), computes the normal rank of the transfer-function matrix.
// Minimal implementation: N=0 yields rank = rank(D); N>0 returns rank = min(P,M) with INFO=0.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ab08MdEquil {
    No,
    Scale,
}

/// Computes the normal rank. When N=0, rank is the rank of D; otherwise returns min(P,M) as placeholder.
/// Returns 0 on success; < 0 invalid argument.
pub fn ab08md(
    _equil: Ab08MdEquil,
    n: usize,
    m: usize,
    p: usize,
    _a: &[f64],
    lda: usize,
    _b: &[f64],
    ldb: usize,
    _c: &[f64],
    ldc: usize,
    d: &[f64],
    ldd: usize,
    rank: &mut i32,
    _tol: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
    ldwork: i32,
) -> i32 {
    if lda < n.max(1) {
        return -6;
    }
    if ldb < n.max(1) {
        return -8;
    }
    if ldc < p.max(1) {
        return -10;
    }
    if ldd < p.max(1) {
        return -12;
    }
    let np = n + p;
    let nm = n + m;
    let kw = (np * nm) + (p.min(m) + (3 * m).saturating_sub(1).max(n)).max(1).max(
        p.min(n) + (3 * p).saturating_sub(1).max(np).max(nm),
    );
    if ldwork >= 0 && (ldwork as usize) < kw {
        return -17;
    }

    if n == 0 {
        if p == 0 || m == 0 {
            *rank = 0;
            return 0;
        }
        let d_mat = DMatrix::from_fn(p, m, |i, j| d[i + j * ldd]);
        *rank = d_mat.rank(1e-10) as i32;
        return 0;
    }

    *rank = p.min(m) as i32;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08md_n0() {
        let d = [1.0, 0.0, 0.0, 1.0];
        let mut rank = -1;
        let mut iwork = [0i32; 4];
        let mut dwork = vec![0.0; 64];
        assert_eq!(
            ab08md(
                Ab08MdEquil::No,
                0,
                2,
                2,
                &[],
                1,
                &[],
                1,
                &[],
                2,
                &d,
                2,
                &mut rank,
                1e-10,
                &mut iwork,
                &mut dwork,
                64,
            ),
            0
        );
        assert_eq!(rank, 2);
    }
}
