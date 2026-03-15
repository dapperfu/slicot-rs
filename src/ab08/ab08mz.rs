//! AB08MZ — Normal rank of transfer-function matrix (complex) (SLICOT AB08MZ)
//!
//! Minimal: returns INFO=0 for trivial dimensions.

/// Returns 0 when N=M=P=0; otherwise placeholder INFO=0.
pub fn ab08mz(n: usize, m: usize, _p: usize) -> i32 {
    if n == 0 && m == 0 {
        return 0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab08mz_trivial() {
        assert_eq!(ab08mz(0, 0, 0), 0);
    }
}
