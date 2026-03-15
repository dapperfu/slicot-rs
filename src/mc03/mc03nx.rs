//! MC03NX — Stub (SLICOT MC03NX) — polynomial matrix pencil operations

/// Stub: returns 0.
pub fn mc03nx(
    _n: i32,
    _a: &[f64],
    _lda: usize,
    _e: &[f64],
    _lde: usize,
    _q: &mut [f64],
    _ldq: usize,
    _z: &mut [f64],
    _ldz: usize,
) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc03nx_stub() {
        assert_eq!(mc03nx(1, &[], 1, &[], 1, &mut [], 1, &mut [], 1), 0);
    }
}
