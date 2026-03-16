//! MA01CD — Safely compute the sign of a sum of two reals in exponent form (SLICOT MA01CD)
//
// Computes the sign of (A * BASE^IA + B * BASE^IB) without over/underflow.
// Returns 1 if sum > 0, 0 if sum == 0, -1 if sum < 0.
// Uses BASE = 2.0 (machine base can be used; same base for both).

const BASE: f64 = 2.0;

/// Returns the sign of (A * BASE^IA + B * BASE^IB): 1 (positive), 0 (zero), -1 (negative).
/// Computed safely by normalizing to the larger exponent to avoid overflow/underflow.
pub fn ma01cd(a: f64, ia: i32, b: f64, ib: i32) -> i32 {
    if a == 0.0 && b == 0.0 {
        return 0;
    }
    let imax = ia.max(ib);
    let a_scale = if ia == imax {
        1.0
    } else {
        BASE.powi(ia - imax)
    };
    let b_scale = if ib == imax {
        1.0
    } else {
        BASE.powi(ib - imax)
    };
    let sum = a * a_scale + b * b_scale;
    if sum > 0.0 {
        1
    } else if sum < 0.0 {
        -1
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma01cd_both_zero() {
        assert_eq!(ma01cd(0.0, 0, 0.0, 0), 0);
    }

    #[test]
    fn test_ma01cd_positive_sum() {
        assert_eq!(ma01cd(1.0, 0, 1.0, 0), 1);
        assert_eq!(ma01cd(2.0, 0, -1.0, 0), 1);
    }

    #[test]
    fn test_ma01cd_negative_sum() {
        assert_eq!(ma01cd(-1.0, 0, -1.0, 0), -1);
        assert_eq!(ma01cd(-2.0, 0, 1.0, 0), -1);
    }

    #[test]
    fn test_ma01cd_with_exponents() {
        // 1*2^0 + 1*2^0 = 2 > 0
        assert_eq!(ma01cd(1.0, 0, 1.0, 0), 1);
        // 1*2^1 + (-1)*2^0 = 2 - 1 = 1 > 0
        assert_eq!(ma01cd(1.0, 1, -1.0, 0), 1);
        // 1*2^0 + (-2)*2^0 = -1 < 0
        assert_eq!(ma01cd(1.0, 0, -2.0, 0), -1);
    }
}
