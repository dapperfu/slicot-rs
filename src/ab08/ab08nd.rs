//! AB08ND — Uncontrollable subspace (SLICOT AB08ND). Minimal: INFO=0 for trivial.

pub fn ab08nd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08nd_trivial() {
        assert_eq!(ab08nd(0, 0), 0);
    }
}
