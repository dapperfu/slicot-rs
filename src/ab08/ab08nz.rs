//! AB08NZ — Minimal stub (SLICOT AB08NZ). Returns INFO=0 for trivial.

pub fn ab08nz(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08nz_trivial() {
        assert_eq!(ab08nz(0, 0), 0);
    }
}
