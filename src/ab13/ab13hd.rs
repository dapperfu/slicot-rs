//! AB13HD — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab13hd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13hd_trivial() {
        assert_eq!(ab13hd(0, 0), 0);
    }
}
