//! MA01AD — Complex square root in real arithmetic (SLICOT MA01AD)
//
// Computes Y = sqrt(X) with YR >= 0 and sign(YI) = sign(XI).

/// Computes the complex square root (yr, yi) of (xr, xi). Returns 0. yr >= 0, sign(yi) = sign(xi).
pub fn ma01ad(xr: f64, xi: f64, yr: &mut f64, yi: &mut f64) -> i32 {
    let half = 0.5;
    let r = (xr * xr + xi * xi).sqrt();
    let mut s = (half * (r + xr.abs())).sqrt();
    if s == 0.0 {
        *yr = 0.0;
        *yi = 0.0;
        return 0;
    }
    if xr >= 0.0 {
        *yr = s;
    }
    if xi < 0.0 {
        s = -s;
    }
    if xr <= 0.0 {
        *yi = s;
        if xr < 0.0 {
            *yr = half * (xi / s);
        }
    } else {
        *yi = half * (xi / *yr);
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma01ad_real_positive() {
        let mut yr = 0.0;
        let mut yi = 0.0;
        assert_eq!(ma01ad(4.0, 0.0, &mut yr, &mut yi), 0);
        assert!((yr - 2.0).abs() < 1e-15);
        assert!((yi - 0.0).abs() < 1e-15);
    }
}
