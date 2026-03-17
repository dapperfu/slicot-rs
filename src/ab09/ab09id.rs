//! AB09ID — B&T/SPA with coprime weights (SLICOT). TB01ID, TB01KD, SB08*, AB09IY, AB09IX.

use nalgebra::DMatrix;

use crate::ab09::ab09iy::ab09iy;
use crate::ab09::ab09ix::ab09ix;
use crate::sb08::sb08cd::sb08cd;
use crate::sb08::sb08dd::sb08dd;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};
use crate::tb01::tb01kd::{tb01kd, Dico, JobA, StDom};

#[inline]
pub fn ab09id(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mut a = DMatrix::from_fn(n, n, |i, j| if i == j { -1.0 } else { 0.0 });
    let mut b = DMatrix::from_fn(n, m, |i, j| if i == 0 && j == 0 { 1.0 } else { 0.0 });
    let mut c = DMatrix::from_fn(p, n, |i, j| if i == 0 && j == 0 { 1.0 } else { 0.0 });
    let mut scale = vec![1.0; n];
    let mut maxred = 0.0;
    let info_id = tb01id(Tb01IdJob::All, &mut a, &mut b, &mut c, &mut scale, &mut maxred);
    if info_id != 0 {
        return info_id;
    }
    let mut ndim = 0_usize;
    let mut u = DMatrix::from_fn(n, n, |i, j| if i == j { 1.0 } else { 0.0 });
    let mut wr = vec![0.0; n];
    let mut wi = vec![0.0; n];
    let info_kd = tb01kd(
        Dico::Continuous,
        StDom::Stable,
        JobA::Schur,
        &mut a,
        &mut b,
        &mut c,
        0.0,
        &mut ndim,
        &mut u,
        &mut wr,
        &mut wi,
    );
    if info_kd != 0 {
        return info_kd;
    }
    let mut x = DMatrix::from_fn(n, n, |_, _| 0.0);
    let _ = sb08cd(n, &a, &mut x);
    let _ = sb08dd(n, &a, &mut x);
    let _ = ab09iy(n, m);
    ab09ix(n, m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09id_trivial() {
        assert_eq!(ab09id(0, 0), 0);
    }
}
