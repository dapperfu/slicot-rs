//! AB09JX — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab09jx(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09jx_trivial() {
        assert_eq!(ab09jx(0, 0), 0);
    }
}
