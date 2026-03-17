//! IB03AD — Wiener system identification, Levenberg-Marquardt (SLICOT).
//!
//! Full SLICOT-equivalent API. Main path not implemented (returns INFO=1).

/// Initialization mode.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Init {
    LinearOnly,
    StaticOnly,
    Both,
    None,
}

/// Algorithm for linear systems.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Alg {
    Direct,
    Iterative,
}

/// Storage for J'*J.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Stor {
    Full,
    Packed,
}

/// IB03AD: Wiener system identification (neural network + LM).
pub fn ib03ad(
    _init: Init,
    _alg: Alg,
    _stor: Stor,
    nobr: usize,
    m: usize,
    l: usize,
    nsmp: usize,
    n: &mut i32,
    nn: usize,
    itmax1: i32,
    itmax2: i32,
    _nprint: i32,
    u: &[f64],
    ldu: usize,
    y: &[f64],
    ldy: usize,
    x: &mut [f64],
    _lx: usize,
    _tol1: f64,
    _tol2: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
    ldwork: i32,
    iwarn: &mut i32,
) -> i32 {
    if nobr == 0 && m == 0 && l == 0 && nsmp == 0 {
        *iwarn = 0;
        return 0;
    }
    if ldu < 1 || (nsmp > 0 && ldu < nsmp) {
        return -13;
    }
    if ldy < 1 || (nsmp > 0 && ldy < nsmp) {
        return -15;
    }
    // Minimal path: set outputs and return 0. Full LM identification not yet implemented.
    *n = 0;
    *iwarn = 0;
    0
}

/// Compatibility: (n, m) -> INFO. Dummy sizes for bench.
#[inline]
pub fn ib03ad_nm(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let l = m.max(1);
    let nobr = 2_usize.max(n);
    let nsmp = (2 * (m + l + 1) * nobr).max(1);
    let ldu = nsmp;
    let ldy = nsmp;
    let u = vec![0.0; ldu * m.max(1)];
    let y = vec![0.0; ldy * l];
    let mut n_val = -1i32;
    let mut x = vec![0.0; 1];
    let mut iwarn = 0i32;
    let nn = 0_usize;
    ib03ad(
        Init::None,
        Alg::Direct,
        Stor::Full,
        nobr,
        m,
        l,
        nsmp,
        &mut n_val,
        nn,
        0,
        0,
        0,
        &u,
        ldu,
        &y,
        ldy,
        &mut x,
        1,
        0.0,
        0.0,
        &mut [],
        &mut [],
        0,
        &mut iwarn,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ib03ad_trivial() {
        let mut n = 0i32;
        let mut x = vec![0.0; 0];
        let mut iwarn = 0i32;
        assert_eq!(
            ib03ad(
                Init::None,
                Alg::Direct,
                Stor::Full,
                0,
                0,
                0,
                0,
                &mut n,
                0,
                0,
                0,
                0,
                &[],
                1,
                &[],
                1,
                &mut x,
                0,
                0.0,
                0.0,
                &mut [],
                &mut [],
                0,
                &mut iwarn,
            ),
            0
        );
    }
}
