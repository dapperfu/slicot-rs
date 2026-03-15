//! IB03AD — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ib03ad(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ib03ad_trivial() {
        assert_eq!(ib03ad(0, 0), 0);
    }
}
