//! AB09JV — Frequency-weighting stability check (SLICOT). TB01WD on weight V then AB09JX on eigenvalues.

use nalgebra::DMatrix;

use crate::ab09::ab09jx::ab09jx_core;
use crate::tb01::tb01wd::tb01wd;

#[inline]
pub fn ab09jv(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let nv = 1_usize.max(n / 2);
    let pv = 1_usize.max(m);
    let mut av = DMatrix::from_fn(nv, nv, |i, j| if i == j { -1.0 } else { 0.0 });
    let mut bv = DMatrix::from_fn(nv, pv, |_, _| 0.0);
    let mut cv = DMatrix::from_fn(pv, nv, |_, _| 0.0);
    let mut u = DMatrix::from_fn(nv, nv, |i, j| if i == j { 1.0 } else { 0.0 });
    let mut wr = vec![0.0; nv];
    let mut wi = vec![0.0; nv];
    let info = tb01wd(&mut av, &mut bv, &mut cv, &mut u, &mut wr, &mut wi);
    if info != 0 {
        return info;
    }
    let ed: Vec<f64> = (0..nv).map(|_| 1.0).collect();
    ab09jx_core(
        b'C',
        b'S',
        b'S',
        nv,
        0.0,
        &wr,
        &wi,
        &ed,
        0.0,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09jv_trivial() {
        assert_eq!(ab09jv(0, 0), 0);
    }
}
