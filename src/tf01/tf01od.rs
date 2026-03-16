//! TF01OD — Block Hankel expansion of a multivariable parameter sequence (SLICOT TF01OD)
//!
//! T(r,c) block = M(r+c-1). Row block r, column block c.

use nalgebra::DMatrix;

/// H is NH1×((NR+NC-1)*NH2); each M(k) is NH1×NH2 at H(0..NH1, (k-1)*NH2..k*NH2).
/// T is (NH1*NR)×(NH2*NC); T(i,j) = M(row_block + col_block - 1)(i mod NH1, j mod NH2).
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tf01od(
    nh1: usize,
    nh2: usize,
    nr: usize,
    nc: usize,
    h: &DMatrix<f64>,
    t: &mut DMatrix<f64>,
) -> i32 {
    if h.nrows() < nh1 || h.ncols() < (nr + nc - 1) * nh2 {
        return -6;
    }
    if t.nrows() < nh1 * nr || t.ncols() < nh2 * nc {
        return -8;
    }
    for r in 0..nr {
        for c in 0..nc {
            let k = r + c;
            for i in 0..nh1 {
                for j in 0..nh2 {
                    let hi = i;
                    let hj = k * nh2 + j;
                    t[(r * nh1 + i, c * nh2 + j)] = if hj < h.ncols() { h[(hi, hj)] } else { 0.0 };
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
    fn test_tf01od_smoke() {
        // NH1=2, NH2=2, NR=2, NC=2 => H is 2 x (2+2-1)*2 = 2x6
        let h = DMatrix::from_row_slice(2, 6, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);
        let mut t = DMatrix::zeros(4, 4);
        assert_eq!(tf01od(2, 2, 2, 2, &h, &mut t), 0);
        assert!((t[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((t[(0, 2)] - 3.0).abs() < 1e-10);
    }
}
