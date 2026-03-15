//! SB01BX — Choose closest real or complex conjugate eigenvalue(s) to a given value (SLICOT SB01BX)
//!
//! Reorders WR, WI so the selected eigenvalue(s) are in the last position(s).
//! Outputs S and P: for real, S = P = selected eigenvalue; for complex pair, S = sum, P = product.

/// Selects a real eigenvalue or a pair of complex conjugate eigenvalues at minimal
/// "distance" to (xr, xi). Uses |x| + |y| for complex distance.
/// On exit, WR and WI are reordered so the selected eigenvalue(s) are last;
/// (s, p) are the selected real eigenvalue or (sum, product) of the complex pair.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid (e.g. N < 1).
pub fn sb01bx(
    reig: bool,
    n: usize,
    xr: f64,
    xi: f64,
    wr: &mut [f64],
    wi: &mut [f64],
    s: &mut f64,
    p: &mut f64,
) -> i32 {
    if n < 1 {
        return -2;
    }
    if wr.len() < n || wi.len() < n {
        return -5;
    }
    if reig {
        let mut best = 0;
        let mut best_dist = (wr[0] - xr).abs();
        for i in 1..n {
            let d = (wr[i] - xr).abs();
            if d < best_dist {
                best_dist = d;
                best = i;
            }
        }
        *s = wr[best];
        *p = wr[best];
        // move selected to end
        if best != n - 1 {
            wr.swap(best, n - 1);
            wi.swap(best, n - 1);
        }
    } else {
        // complex pair: find pair with minimal "distance" (|wr-xr| + |wi-xi|) for first of pair
        let mut best = 0;
        let mut best_dist = (wr[0] - xr).abs() + (wi[0] - xi).abs();
        let mut i = 0;
        while i < n {
            let d = (wr[i] - xr).abs() + (wi[i] - xi).abs();
            if d < best_dist {
                best_dist = d;
                best = i;
            }
            if i + 1 < n && wi[i + 1] != 0.0 {
                i += 2;
            } else {
                i += 1;
            }
        }
        let idx2 = if best + 1 < n && wi[best + 1] != 0.0 {
            best + 1
        } else {
            best
        };
        let w1r = wr[best];
        let w1i = wi[best];
        let (w2r, w2i) = if idx2 != best {
            (wr[idx2], wi[idx2])
        } else {
            (w1r, -w1i)
        };
        *s = w1r + w2r;
        *p = w1r * w2r - w1i * w2i;
        // move selected pair to last two positions
        if n >= 2 && (best != n - 2 || idx2 != n - 1) {
            let (last_r0, last_i0) = (wr[n - 2], wi[n - 2]);
            let (last_r1, last_i1) = (wr[n - 1], wi[n - 1]);
            wr[n - 2] = w1r;
            wi[n - 2] = w1i;
            wr[n - 1] = w2r;
            wi[n - 1] = w2i;
            wr[best] = last_r0;
            wi[best] = last_i0;
            if idx2 != best {
                wr[idx2] = last_r1;
                wi[idx2] = last_i1;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb01bx_real_n1() {
        let mut wr = [3.0];
        let mut wi = [0.0];
        let mut s = 0.0;
        let mut p = 0.0;
        assert_eq!(sb01bx(true, 1, 2.0, 0.0, &mut wr, &mut wi, &mut s, &mut p), 0);
        assert_eq!(s, 3.0);
        assert_eq!(p, 3.0);
    }

    #[test]
    fn test_sb01bx_real_closest() {
        let mut wr = [1.0, 5.0, 2.0];
        let mut wi = [0.0, 0.0, 0.0];
        let mut s = 0.0;
        let mut p = 0.0;
        assert_eq!(sb01bx(true, 3, 2.0, 0.0, &mut wr, &mut wi, &mut s, &mut p), 0);
        assert_eq!(s, 2.0);
        assert_eq!(p, 2.0);
        assert_eq!(wr[2], 2.0);
    }

    #[test]
    fn test_sb01bx_invalid_n() {
        let mut wr = [1.0];
        let mut wi = [0.0];
        let mut s = 0.0;
        let mut p = 0.0;
        assert!(sb01bx(true, 0, 0.0, 0.0, &mut wr, &mut wi, &mut s, &mut p) < 0);
    }
}
