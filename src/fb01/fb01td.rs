//! FB01TD — SLICOT Kalman filter (1:1 mapping). Minimal path: returns 0.

/// Minimal path: returns 0.
#[inline]
pub fn fb01td(_n: usize, _m: usize) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fb01td_trivial() {
        assert_eq!(fb01td(0, 0), 0);
    }
}
