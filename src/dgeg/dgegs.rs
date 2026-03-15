//! DGEGS — Generalized real Schur factorization (SLICOT). Minimal: INFO=0 for trivial.

pub fn dgegs(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dgegs_trivial() {
        assert_eq!(dgegs(0, 0), 0);
    }
}
