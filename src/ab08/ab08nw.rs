//! AB08NW — Minimal stub (SLICOT AB08NW). Returns INFO=0 for trivial.

pub fn ab08nw(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08nw_trivial() {
        assert_eq!(ab08nw(0, 0), 0);
    }
}
