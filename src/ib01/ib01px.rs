//! IB01PX — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ib01px(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ib01px_trivial() {
        assert_eq!(ib01px(0, 0), 0);
    }
}
