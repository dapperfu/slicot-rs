//! AB09GD — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab09gd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09gd_trivial() {
        assert_eq!(ab09gd(0, 0), 0);
    }
}
