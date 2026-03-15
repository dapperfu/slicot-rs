//! IB01OD — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ib01od(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ib01od_trivial() {
        assert_eq!(ib01od(0, 0), 0);
    }
}
