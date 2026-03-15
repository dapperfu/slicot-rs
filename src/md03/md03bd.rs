//! MD03BD — Stub (SLICOT MD03BD)

pub fn md03bd(
    _m: i32,
    _n: i32,
    _x: &[f64],
    _f: &mut [f64],
    _fjac: &mut [f64],
    _ldfjac: usize,
    _dwork: &mut [f64],
) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_md03bd_stub() {
        let mut f = [0.0; 2];
        let mut fjac = [0.0; 4];
        let mut dwork = [0.0; 20];
        assert_eq!(md03bd(1, 1, &[], &mut f, &mut fjac, 1, &mut dwork), 0);
    }
}
