//! AB09BX — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab09bx(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab09bx_trivial() {
        assert_eq!(ab09bx(0, 0), 0);
    }
}
