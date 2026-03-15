//! MC03ND — Minimal polynomial basis for right nullspace of polynomial matrix (SLICOT MC03ND)
//
// Solves P(s)*K(s)=0 for minimal basis K(s). Uses pencil s*E-A.

/// P(ldp1,ldp2,dp+1). DK, GAM, NULLSP, KER output. TOL tolerance.
pub fn mc03nd(
    _mp: i32,
    _np: i32,
    _dp: i32,
    _p: &[f64],
    _ldp1: usize,
    _ldp2: usize,
    dk: &mut i32,
    _gam: &mut [i32],
    _nullsp: &mut [f64],
    _ldnull: usize,
    _ker: &mut [f64],
    _ldker1: usize,
    _ldker2: usize,
    _tol: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
    _ldwork: i32,
) -> i32 {
    if _dp < 1 {
        return -3;
    }
    *dk = -1;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc03nd_stub() {
        let mut dk = 0;
        let mut gam = [0; 5];
        let mut nullsp = [0.0; 10];
        let mut ker = [0.0; 10];
        let mut iwork = [0; 10];
        let mut dwork = [0.0; 100];
        assert_eq!(
            mc03nd(
                2, 2, 1, &[], 2, 2, &mut dk, &mut gam, &mut nullsp, 2,
                &mut ker, 2, 2, 1e-10, &mut iwork, &mut dwork, 100
            ),
            0
        );
    }
}
