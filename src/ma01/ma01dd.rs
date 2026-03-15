//! MA01DD — Approximate symmetric chordal metric for two complex numbers (SLICOT MA01DD)
//
// D = min(|A1 - A2|, |1/A1 - 1/A2|). Aj = ARj + i*AIj.

fn dlapy2(x: f64, y: f64) -> f64 {
    x.hypot(y)
}

/// Computes the chordal metric D for (AR1, AI1) and (AR2, AI2). Returns 0.
pub fn ma01dd(
    ar1: f64,
    ai1: f64,
    ar2: f64,
    ai2: f64,
    eps: f64,
    safemn: f64,
    d: &mut f64,
) -> i32 {
    const ZERO: f64 = 0.0;
    const ONE: f64 = 1.0;
    const TWO: f64 = 2.0;
    const FOUR: f64 = 4.0;

    let par = FOUR - TWO * eps;
    let mut big = par / safemn;
    if big * safemn > par {
        big = ONE / safemn;
    }

    let mx1 = ar1.abs().max(ai1.abs());
    let mx2 = ar2.abs().max(ai2.abs());
    let mx = mx1.max(mx2);

    if mx == ZERO {
        *d = ZERO;
        return 0;
    }

    let d1 = if mx < big {
        if mx2 == ZERO {
            dlapy2(ar1, ai1)
        } else if mx1 == ZERO {
            dlapy2(ar2, ai2)
        } else {
            dlapy2(ar1 - ar2, ai1 - ai2)
        }
    } else {
        big
    };

    let d2 = if mx > ONE / big {
        let ap1 = dlapy2(ar1, ai1);
        let ap2 = dlapy2(ar2, ai2);
        if mx1 <= big && mx2 <= big {
            let re = (ar1 / ap1) / ap1 - (ar2 / ap2) / ap2;
            let im = (ai2 / ap2) / ap2 - (ai1 / ap1) / ap1;
            dlapy2(re, im)
        } else if mx1 <= big {
            ONE / ap1
        } else if mx2 <= big {
            ONE / ap2
        } else {
            ZERO
        }
    } else {
        big
    };

    *d = d1.min(d2);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma01dd_same() {
        let eps = 1e-15_f64;
        let safemn = 1e-308_f64;
        let mut d = -1.0;
        assert_eq!(ma01dd(1.0, 0.0, 1.0, 0.0, eps, safemn, &mut d), 0);
        assert!(d >= 0.0 && d < 1e-14);
    }

    #[test]
    fn test_ma01dd_diff() {
        let eps = 1e-15_f64;
        let safemn = 1e-308_f64;
        let mut d = -1.0;
        assert_eq!(ma01dd(1.0, 0.0, 2.0, 0.0, eps, safemn, &mut d), 0);
        assert!(d > 0.0 && d < 2.0);
    }
}
