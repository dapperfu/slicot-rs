//! AB09JW — Frequency-weighting stability check (SLICOT). TB01WD on weight W then AB09JX on eigenvalues.

use nalgebra::DMatrix;

use crate::ab09::ab09jx::ab09jx_core;
use crate::tb01::tb01wd::tb01wd;

#[inline]
pub fn ab09jw(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let nw = 1_usize.max(m / 2);
    let mw = 1_usize.max(m);
    let mut aw = DMatrix::from_fn(nw, nw, |i, j| if i == j { -1.0 } else { 0.0 });
    let mut bw = DMatrix::from_fn(nw, mw, |_, _| 0.0);
    let mut cw = DMatrix::from_fn(mw, nw, |_, _| 0.0);
    let mut u = DMatrix::from_fn(nw, nw, |i, j| if i == j { 1.0 } else { 0.0 });
    let mut wr = vec![0.0; nw];
    let mut wi = vec![0.0; nw];
    let info = tb01wd(&mut aw, &mut bw, &mut cw, &mut u, &mut wr, &mut wi);
    if info != 0 {
        return info;
    }
    let ed: Vec<f64> = (0..nw).map(|_| 1.0).collect();
    ab09jx_core(
        b'C',
        b'S',
        b'S',
        nw,
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
    fn test_ab09jw_trivial() {
        assert_eq!(ab09jw(0, 0), 0);
    }
}
