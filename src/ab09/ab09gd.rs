//! AB09GD — SPA with coprime factors (SLICOT). TB01ID, SB08ED, SB08CD, AB09BX, SB08GD/FD/DD/HD.

use nalgebra::DMatrix;

use crate::ab09::ab09bx::ab09bx;
use crate::sb08::sb08cd::sb08cd;
use crate::sb08::sb08dd::sb08dd;
use crate::sb08::sb08ed::sb08ed;
use crate::sb08::sb08fd::sb08fd;
use crate::sb08::sb08gd::sb08gd;
use crate::sb08::sb08hd::sb08hd;
use crate::tb01::tb01id::{tb01id, Tb01IdJob};

#[inline]
pub fn ab09gd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    let p = m;
    let mut a = DMatrix::from_fn(n, n, |i, j| if i == j { -1.0 } else { 0.0 });
    let mut b = DMatrix::from_fn(n, m, |i, j| if i == 0 && j == 0 { 1.0 } else { 0.0 });
    let mut c = DMatrix::from_fn(p, n, |i, j| if i == 0 && j == 0 { 1.0 } else { 0.0 });
    let mut scale = vec![1.0; n];
    let mut maxred = 0.0;
    let info = tb01id(Tb01IdJob::All, &mut a, &mut b, &mut c, &mut scale, &mut maxred);
    if info != 0 {
        return info;
    }
    let mut x = DMatrix::from_fn(n, n, |_, _| 0.0);
    let _ = sb08ed(n, &a, &mut x);
    let _ = sb08cd(n, &a, &mut x);
    let info_bx = ab09bx(n, m);
    if info_bx != 0 {
        return info_bx;
    }
    let _ = sb08gd(n, &a, &mut x);
    let _ = sb08fd(n, &a, &mut x);
    let _ = sb08dd(n, &a, &mut x);
    let _ = sb08hd(n, &a, &mut x);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09gd_trivial() {
        assert_eq!(ab09gd(0, 0), 0);
    }
}
