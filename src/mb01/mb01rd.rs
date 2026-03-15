//! MB01RD — R = alpha*R + beta*op(A)*X*op(A)' (SLICOT MB01RD)
// Same formula as MB01RU; full version. We delegate to mb01ru.

use super::mb01ru::{Mb01RuTrans, Mb01RuUplo};

/// Same operation as MB01RU. Delegates to mb01ru.
pub fn mb01rd(
    uplo: Mb01RuUplo,
    trans: Mb01RuTrans,
    m: usize,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    a: &[f64],
    lda: usize,
    x: &[f64],
    ldx: usize,
    dwork: &mut [f64],
) -> i32 {
    super::mb01ru::mb01ru(uplo, trans, m, n, alpha, beta, r, ldr, a, lda, x, ldx, dwork)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::mb01ru::{Mb01RuTrans, Mb01RuUplo};

    #[test]
    fn test_mb01rd_upper_notrans() {
        let m = 2;
        let n = 2;
        let mut r = vec![1.0, 0.0, 0.0, 1.0];
        let a = [1.0, 0.0, 0.0, 1.0];
        let x = [1.0, 0.0, 0.0, 1.0];
        let mut dwork = vec![0.0; m * n];
        assert_eq!(
            mb01rd(
                Mb01RuUplo::Upper,
                Mb01RuTrans::NoTrans,
                m,
                n,
                1.0,
                1.0,
                &mut r,
                2,
                &a,
                2,
                &x,
                2,
                &mut dwork,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
