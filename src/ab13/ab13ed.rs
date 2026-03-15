//! AB13ED — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab13ed(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13ed_trivial() {
        assert_eq!(ab13ed(0, 0), 0);
    }
}
