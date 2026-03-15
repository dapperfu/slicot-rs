//! AB09IY — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab09iy(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09iy_trivial() {
        assert_eq!(ab09iy(0, 0), 0);
    }
}
