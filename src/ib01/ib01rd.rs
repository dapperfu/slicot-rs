//! IB01RD — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ib01rd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ib01rd_trivial() {
        assert_eq!(ib01rd(0, 0), 0);
    }
}
