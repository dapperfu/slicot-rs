//! AG08BD — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ag08bd(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ag08bd_trivial() {
        assert_eq!(ag08bd(0, 0), 0);
    }
}
