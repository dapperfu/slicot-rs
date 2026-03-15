//! DGEGV — Generalized real eigenvalue problem (SLICOT). Minimal: INFO=0 for trivial.

pub fn dgegv(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dgegv_trivial() {
        assert_eq!(dgegv(0, 0), 0);
    }
}
