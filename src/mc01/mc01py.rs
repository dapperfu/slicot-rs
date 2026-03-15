//! MC01PY — Coefficients of real polynomial from zeros (SLICOT MC01PY)
//
// P(x) = (x-r1)(x-r2)...(x-rK). Output P in decreasing powers of x. Complex conjugate pairs must be consecutive.

/// REZ, IMZ are real and imaginary parts of zeros. P output: P[0]=x^K coef, P[K]=constant.
pub fn mc01py(k: i32, rez: &[f64], imz: &[f64], p: &mut [f64], dwork: &mut [f64]) -> i32 {
    let k = k as usize;
    if k > rez.len() || k > imz.len() {
        return -2;
    }
    if k + 1 > p.len() {
        return -4;
    }
    if k > 0 && k > dwork.len() {
        return -5;
    }
    for i in 1..k {
        if imz[i].abs() > 1e-15 {
            let conj = (rez[i] - rez[i - 1]).abs() < 1e-15 && (imz[i] + imz[i - 1]).abs() < 1e-15;
            if !conj {
                return i as i32 + 1;
            }
        }
    }
    if k == 0 {
        p[0] = 1.0;
        return 0;
    }
    p[0] = 1.0;
    let mut deg = 0_usize;
    let mut i = 0_usize;
    while i < k {
        let (cr, ci) = (rez[i], imz[i]);
        if ci.abs() < 1e-15 {
            if deg == 0 {
                p[0] = 1.0;
                p[1] = -cr;
                deg = 1;
            } else {
                dwork[0] = -cr;
                for j in (0..=deg).rev() {
                    dwork[j + 1] = p[j];
                }
                deg += 1;
                for j in 0..=deg {
                    p[j] = dwork[j];
                }
            }
            i += 1;
        } else {
            let a = -2.0 * cr;
            let b = cr * cr + ci * ci;
            for j in 0..=deg + 2 {
                dwork[j] = 0.0;
            }
            for j in 0..=deg {
                dwork[j + 2] += p[j];
                dwork[j + 1] += a * p[j];
                dwork[j] += b * p[j];
            }
            deg += 2;
            for j in 0..=deg {
                p[j] = dwork[j];
            }
            i += 2;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01py_real_zero() {
        let mut p = [0.0; 2];
        let mut dwork = [0.0; 1];
        assert_eq!(mc01py(1, &[2.0], &[0.0], &mut p, &mut dwork), 0);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - (-2.0)).abs() < 1e-10);
    }
}
