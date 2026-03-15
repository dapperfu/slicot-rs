//! AB08NY — Minimal stub (SLICOT AB08NY). Returns INFO=0 for trivial.

pub fn ab08ny(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08ny_trivial() {
        assert_eq!(ab08ny(0, 0), 0);
    }
}
