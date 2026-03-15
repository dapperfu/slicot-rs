//! MC03NY — Stub (SLICOT MC03NY) — polynomial matrix pencil operations

/// Stub: returns 0.
pub fn mc03ny(
    _n: i32,
    _a: &[f64],
    _lda: usize,
    _e: &[f64],
    _lde: usize,
) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc03ny_stub() {
        assert_eq!(mc03ny(1, &[], 1, &[], 1), 0);
    }
}
