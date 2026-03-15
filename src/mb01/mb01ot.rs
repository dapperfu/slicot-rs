//! MB01OT — R := alpha*R + beta*E*T' + beta*T*E' (SLICOT MB01OT)
// E and T upper triangular.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OtUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OtTrans {
    NoTrans,
    Trans,
}

/// Overwrites the triangle of R.
pub fn mb01ot(
    uplo: Mb01OtUplo,
    trans: Mb01OtTrans,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    e: &[f64],
    lde: usize,
    t: &[f64],
    ldt: usize,
) -> i32 {
    if !matches!(uplo, Mb01OtUplo::Upper | Mb01OtUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01OtTrans::NoTrans | Mb01OtTrans::Trans) {
        return -2;
    }
    if ldr < n.max(1) || lde < n.max(1) || ldt < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01OtUplo::Upper && i <= j) || (uplo == Mb01OtUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let e_mat = DMatrix::from_fn(n, n, |i, j| if i <= j { e[i + j * lde] } else { 0.0 });
    let t_mat = DMatrix::from_fn(n, n, |i, j| if i <= j { t[i + j * ldt] } else { 0.0 });
    let update = match trans {
        Mb01OtTrans::NoTrans => &e_mat * t_mat.transpose() + &t_mat * e_mat.transpose(),
        Mb01OtTrans::Trans => e_mat.transpose() * &t_mat + t_mat.transpose() * &e_mat,
    };

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01OtUplo::Upper && i <= j) || (uplo == Mb01OtUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01OtUplo::Upper && i <= j) || (uplo == Mb01OtUplo::Lower && i >= j) {
                r[i + j * ldr] = r_full[(i, j)];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01ot_upper_notrans() {
        let n = 2;
        let mut r = vec![0.0, 0.0, 0.0, 0.0];
        let e = [1.0, 0.0, 0.0, 1.0];
        let t = [1.0, 0.0, 0.0, 1.0];
        assert_eq!(
            mb01ot(
                Mb01OtUplo::Upper,
                Mb01OtTrans::NoTrans,
                n,
                0.0,
                1.0,
                &mut r,
                2,
                &e,
                2,
                &t,
                2,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
