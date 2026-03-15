//! MD03BF — Stub (SLICOT MD03BF)

pub fn md03bf(
    _n: i32,
    _x: &[f64],
    _f: &mut f64,
    _grad: &mut [f64],
    _dwork: &mut [f64],
) -> i32 {
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_md03bf_stub() {
        let mut f = 0.0;
        let mut grad = [0.0; 2];
        let mut dwork = [0.0; 10];
        assert_eq!(md03bf(1, &[], &mut f, &mut grad, &mut dwork), 0);
    }
}
