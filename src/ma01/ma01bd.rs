//! MA01BD — General product of K real scalars without overflow/underflow (SLICOT MA01BD)
//
// Computes product of A(i)^S(i) with S(i)=±1, returning ALPHA/BETA * BASE^SCAL.

/// Computes the scaled product. S[i] must be 1 or -1 (multiply or divide by A).
/// Returns 0 on success; < 0 invalid argument (e.g. invalid UPLO).
pub fn ma01bd(
    base: f64,
    lgbas: f64,
    k: usize,
    s: &[i32],
    a: &[f64],
    inca: i32,
    alpha: &mut f64,
    beta: &mut f64,
    scal: &mut i32,
) -> i32 {
    if k == 0 {
        return 0;
    }
    if inca == 0 {
        return -6;
    }
    if inca < 0 {
        return -6;
    }
    let inca = inca as usize;
    if s.len() < k {
        return -4;
    }
    let need = 1usize.saturating_add((k - 1).saturating_mul(inca));
    if a.len() < need {
        return -5;
    }

    *alpha = 1.0;
    *beta = 1.0;
    *scal = 0;

    for i in 0..k {
        let idx = 1 + i * inca;
        let idx = idx.saturating_sub(1);
        let temp = if idx < a.len() { a[idx] } else { 0.0 };
        let (temp, sl) = if temp != 0.0 {
            let sl = (temp.abs().ln() / lgbas).floor() as i32;
            let t = temp / base / base.powi(sl - 1);
            (t, sl)
        } else {
            (0.0, 0i32)
        };
        if s[i] == 1 {
            *alpha *= temp;
            *scal += sl;
        } else {
            *beta *= temp;
            *scal -= sl;
        }
        if (i + 1) % 10 == 0 {
            if *alpha != 0.0 {
                let sl = ((*alpha).abs().ln() / lgbas).floor() as i32;
                *scal += sl;
                *alpha /= base * base.powi(sl - 1);
            }
            if *beta != 0.0 {
                let sl = ((*beta).abs().ln() / lgbas).floor() as i32;
                *scal -= sl;
                *beta /= base * base.powi(sl - 1);
            }
        }
    }

    if *beta != 0.0 {
        *alpha /= *beta;
        *beta = 1.0;
    }
    if *alpha == 0.0 {
        *scal = 0;
    } else {
        let sl = ((*alpha).abs().ln() / lgbas).floor() as i32;
        *alpha /= base * base.powi(sl - 1);
        *scal += sl;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma01bd_single() {
        let base = 2.0_f64;
        let lgbas = base.ln();
        let s = [1];
        let a = [8.0]; // 2^3
        let mut alpha = 0.0;
        let mut beta = 0.0;
        let mut scal = 0;
        assert_eq!(ma01bd(base, lgbas, 1, &s, &a, 1, &mut alpha, &mut beta, &mut scal), 0);
        assert!((alpha - 1.0).abs() < 1e-14);
        assert!((beta - 1.0).abs() < 1e-14);
        assert_eq!(scal, 3);
    }
}
