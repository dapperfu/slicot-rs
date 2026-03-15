//! MC03MD — P(x) = P1(x)*P2(x) + alpha*P3(x) (SLICOT MC03MD)
//
// Real polynomial matrix operation. P1(rp1×cp1), P2(cp1×cp2), P3(rp1×cp2).
// Coefficients in increasing powers; stored as (row, col, power).

/// P1, P2, P3 are (ldp1, ldp2, deg+1). P3 overwritten with result; DP3 set to output degree.
pub fn mc03md(
    rp1: i32,
    cp1: i32,
    cp2: i32,
    dp1: i32,
    dp2: i32,
    dp3: &mut i32,
    alpha: f64,
    p1: &[f64],
    ldp11: usize,
    ldp12: usize,
    p2: &[f64],
    ldp21: usize,
    ldp22: usize,
    p3: &mut [f64],
    ldp31: usize,
    ldp32: usize,
    dwork: &mut [f64],
) -> i32 {
    let (rp1, cp1, cp2) = (rp1 as usize, cp1 as usize, cp2 as usize);
    let (dp1, dp2) = (dp1 as i32, dp2 as i32);
    if cp1 > dwork.len() {
        return -18;
    }
    let out_deg = if dp1 >= 0 && dp2 >= 0 {
        dp1 + dp2
    } else if dp1 >= 0 {
        dp1
    } else if dp2 >= 0 {
        dp2
    } else {
        -1
    };
    let n_out = (out_deg.max(*dp3).max(0) + 1) as usize;
    if ldp31 < rp1.max(1) || ldp32 < cp2.max(1) {
        return -16;
    }
    let base3 = ldp31 * ldp32;
    for i in 0..rp1 {
        for j in 0..cp2 {
            for k in 0..n_out {
                let idx3 = i + j * ldp31 + k * base3;
                let mut s = if (*dp3 >= 0) && (k <= *dp3 as usize) {
                    alpha * p3[idx3]
                } else {
                    0.0
                };
                if dp1 >= 0 && dp2 >= 0 {
                    let r = 0_i32.max(k as i32 - dp2);
                    let s_ = (k as i32).min(dp1);
                    for kk in r..=s_ {
                        let k1 = kk as usize;
                        let k2 = k - k1;
                        if k2 <= dp2 as usize {
                            for c in 0..cp1 {
                                let a = p1[i + c * ldp11 + k1 * ldp11 * ldp12];
                                let b = p2[c + j * ldp21 + k2 * ldp21 * ldp22];
                                s += a * b;
                            }
                        }
                    }
                }
                p3[idx3] = s;
            }
        }
    }
    *dp3 = out_deg;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc03md_simple() {
        let p1 = [1.0, 1.0];
        let p2 = [1.0, 1.0];
        let mut p3 = [0.0, 0.0, 1.0];
        let mut dp3 = 1;
        let mut dwork = [0.0; 2];
        assert_eq!(
            mc03md(
                1, 1, 1, 1, 1, &mut dp3, 0.0,
                &p1, 1, 1, &p2, 1, 1, &mut p3, 1, 1, &mut dwork
            ),
            0
        );
        assert_eq!(dp3, 2);
        assert!((p3[0] - 1.0).abs() < 1e-10);
        assert!((p3[1] - 2.0).abs() < 1e-10);
        assert!((p3[2] - 1.0).abs() < 1e-10);
    }
}
