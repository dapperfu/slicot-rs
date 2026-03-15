//! AB8NXZ — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ab8nxz(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab8nxz_trivial() {
        assert_eq!(ab8nxz(0, 0), 0);
    }
}
