//! AB09JW — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab09jw(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09jw_trivial() {
        assert_eq!(ab09jw(0, 0), 0);
    }
}
