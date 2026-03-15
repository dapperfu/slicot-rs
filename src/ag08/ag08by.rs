//! AG08BY — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ag08by(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ag08by_trivial() {
        assert_eq!(ag08by(0, 0), 0);
    }
}
