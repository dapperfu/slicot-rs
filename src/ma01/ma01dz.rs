//! MA01DZ — Approximate symmetric chordal metric for two rational complex numbers (SLICOT MA01DZ)
//
// A1 = (AR1+i*AI1)/B1, A2 = (AR2+i*AI2)/B2. D = min(|A1-A2|, |1/A1-1/A2|). Output (D1, D2) = numerator/denominator of D.

fn dlapy2(x: f64, y: f64) -> f64 {
    x.hypot(y)
}

/// Computes chordal metric for (AR1, AI1, B1) and (AR2, AI2, B2). D1, D2 define D (D2=0 or 1).
/// IWARN: 0 = ok, 1 = A1 or A2 is not a number (NaN); D1=D2=0.
pub fn ma01dz(
    ar1: f64,
    ai1: f64,
    b1: f64,
    ar2: f64,
    ai2: f64,
    b2: f64,
    eps: f64,
    safemn: f64,
    d1: &mut f64,
    d2: &mut f64,
    iwarn: &mut i32,
) -> i32 {
    const ZERO: f64 = 0.0;
    const ONE: f64 = 1.0;
    const TWO: f64 = 2.0;
    const FOUR: f64 = 4.0;

    *iwarn = 0;

    let par = FOUR - TWO * eps;
    let mut big = par / safemn;
    if big * safemn > par {
        big = ONE / safemn;
    }

    let mx1 = ar1.abs().max(ai1.abs());
    let mx2 = ar2.abs().max(ai2.abs());

    if b1 == ZERO {
        if mx1 == ZERO {
            *d1 = ZERO;
            *d2 = ZERO;
            *iwarn = 1;
        } else {
            if b2 == ZERO {
                *d1 = ZERO;
                if mx2 == ZERO {
                    *d2 = ZERO;
                    *iwarn = 1;
                } else {
                    *d2 = ONE;
                }
            } else if b2 > ONE {
                if mx2 > b2 / big {
                    *d1 = b2 / dlapy2(ar2, ai2);
                    *d2 = ONE;
                } else {
                    *d1 = ONE;
                    *d2 = ZERO;
                }
            } else if mx2 > ZERO {
                *d1 = b2 / dlapy2(ar2, ai2);
                *d2 = ONE;
            } else {
                *d1 = ONE;
                *d2 = ZERO;
            }
        }
        return 0;
    }

    if b2 == ZERO {
        if mx2 == ZERO {
            *d1 = ZERO;
            *d2 = ZERO;
            *iwarn = 1;
        } else {
            if b1 > ONE {
                if mx1 > b1 / big {
                    *d1 = b1 / dlapy2(ar1, ai1);
                    *d2 = ONE;
                } else {
                    *d1 = ONE;
                    *d2 = ZERO;
                }
            } else if mx1 > ZERO {
                *d1 = b1 / dlapy2(ar1, ai1);
                *d2 = ONE;
            } else {
                *d1 = ONE;
                *d2 = ZERO;
            }
        }
        return 0;
    }

    let (zer1, inf1, ap1) = if b1 >= ONE {
        let ap = dlapy2(ar1 / b1, ai1 / b1);
        (ap < ONE / big, false, ap)
    } else {
        let inf = mx1 > b1 * big;
        let ap = if inf { ZERO } else { dlapy2(ar1 / b1, ai1 / b1) };
        (false, inf, ap)
    };

    let (zer2, inf2, ap2) = if b2 >= ONE {
        let ap = dlapy2(ar2 / b2, ai2 / b2);
        (ap < ONE / big, false, ap)
    } else {
        let inf = mx2 > b2 * big;
        let ap = if inf { ZERO } else { dlapy2(ar2 / b2, ai2 / b2) };
        (false, inf, ap)
    };

    *d2 = ONE;

    if zer1 && zer2 {
        *d1 = ZERO;
    } else if zer1 {
        if !inf2 {
            *d1 = ap2;
        } else {
            *d1 = ONE;
            *d2 = ZERO;
        }
    } else if zer2 {
        if !inf1 {
            *d1 = ap1;
        } else {
            *d1 = ONE;
            *d2 = ZERO;
        }
    } else if inf1 {
        if inf2 {
            *d1 = ZERO;
        } else {
            *d1 = b2 / dlapy2(ar2, ai2);
        }
    } else if inf2 {
        *d1 = b1 / dlapy2(ar1, ai1);
    } else {
        let pr1 = ar1 / b1;
        let pi1 = ai1 / b1;
        let pr2 = ar2 / b2;
        let pi2 = ai2 / b2;
        let d_direct = dlapy2(pr1 - pr2, pi1 - pi2);
        let re_inv = (pr1 / ap1) / ap1 - (pr2 / ap2) / ap2;
        let im_inv = (pi2 / ap2) / ap2 - (pi1 / ap1) / ap1;
        let d_inv = dlapy2(re_inv, im_inv);
        *d1 = d_direct.min(d_inv);
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma01dz_same() {
        let eps = 1e-15_f64;
        let safemn = 1e-308_f64;
        let mut d1 = -1.0_f64;
        let mut d2 = -1.0_f64;
        let mut iwarn = -1_i32;
        ma01dz(1.0, 0.0, 1.0, 1.0, 0.0, 1.0, eps, safemn, &mut d1, &mut d2, &mut iwarn);
        assert_eq!(iwarn, 0);
        assert!(d1 >= 0.0 && d1 < 1e-14);
        assert_eq!(d2, 1.0);
    }
}
