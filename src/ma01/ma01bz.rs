//! MA01BZ — General product of K complex scalars avoiding overflow/underflow (SLICOT MA01BZ)
//
// Product of A(i)^S(i) with S(i)=±1 (multiply or divide). Returns ALPHA/BETA * BASE^SCAL in complex.

const ONE: f64 = 1.0;
const ZERO: f64 = 0.0;

fn cabs(ar: f64, ai: f64) -> f64 {
    ar.hypot(ai)
}

fn cmul(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    (ar * br - ai * bi, ar * bi + ai * br)
}

fn cdiv(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    let n = br * br + bi * bi;
    if n == 0.0 {
        return (ZERO, ZERO);
    }
    ((ar * br + ai * bi) / n, (ai * br - ar * bi) / n)
}

/// Computes the scaled product of complex scalars. S[i] must be 1 or -1 (multiply or divide by A(i)).
/// A is stored as (a_re, a_im); a_re and a_im have stride INCA (so a_re[i*inca], a_im[i*inca] for i in 0..k).
/// Returns 0 on success; < 0 invalid argument.
pub fn ma01bz(
    base: f64,
    k: usize,
    s: &[i32],
    a_re: &[f64],
    a_im: &[f64],
    inca: i32,
    alpha_re: &mut f64,
    alpha_im: &mut f64,
    beta_re: &mut f64,
    beta_im: &mut f64,
    scal: &mut i32,
) -> i32 {
    if k == 0 {
        return 0;
    }
    if inca == 0 {
        return -5;
    }
    let inca_u = inca.unsigned_abs() as usize;
    if s.len() < k {
        return -3;
    }
    let need = 1usize.saturating_add((k - 1).saturating_mul(inca_u));
    if a_re.len() < need || a_im.len() < need {
        return -4;
    }

    let cbase_re = base;
    let cbase_im = ZERO;

    *alpha_re = ONE;
    *alpha_im = ZERO;
    *beta_re = ONE;
    *beta_im = ZERO;
    *scal = 0;

    for i in 0..k {
        let idx = i * inca_u;
        let ar = a_re[idx];
        let ai = a_im[idx];

        if s[i] == 1 {
            let (r, im) = cmul(*alpha_re, *alpha_im, ar, ai);
            *alpha_re = r;
            *alpha_im = im;
        } else {
            if ar == ZERO && ai == ZERO {
                *beta_re = ZERO;
                *beta_im = ZERO;
            } else {
                let (r, im) = cdiv(*alpha_re, *alpha_im, ar, ai);
                *alpha_re = r;
                *alpha_im = im;
            }
        }

        let abs_a = cabs(*alpha_re, *alpha_im);
        if abs_a == ZERO {
            *alpha_re = ZERO;
            *alpha_im = ZERO;
            *scal = 0;
            if cabs(*beta_re, *beta_im) == ZERO {
                return 0;
            }
        } else {
            while cabs(*alpha_re, *alpha_im) < ONE {
                let (r, im) = cmul(*alpha_re, *alpha_im, cbase_re, cbase_im);
                *alpha_re = r;
                *alpha_im = im;
                *scal -= 1;
            }
            while cabs(*alpha_re, *alpha_im) >= base {
                let (r, im) = cdiv(*alpha_re, *alpha_im, cbase_re, cbase_im);
                *alpha_re = r;
                *alpha_im = im;
                *scal += 1;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma01bz_single_multiply() {
        let base = 2.0_f64;
        let s = [1];
        let a_re = [4.0];
        let a_im = [0.0];
        let mut ar = 0.0;
        let mut ai = 0.0;
        let mut br = 0.0;
        let mut bi = 0.0;
        let mut scal = 0;
        assert_eq!(
            ma01bz(base, 1, &s, &a_re, &a_im, 1, &mut ar, &mut ai, &mut br, &mut bi, &mut scal),
            0
        );
        assert!((ar - 1.0).abs() < 1e-14 && ai.abs() < 1e-14);
        assert_eq!(scal, 2); // 4 = base^2
    }
}
