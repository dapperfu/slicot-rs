//! IB01PY — Minimal stub (SLICOT). Returns INFO=0 for trivial.

pub fn ib01py(n: usize, m: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ib01py_trivial() {
        assert_eq!(ib01py(0, 0), 0);
    }
}
