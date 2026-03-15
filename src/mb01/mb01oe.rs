//! MB01OE — R := alpha*R + beta*H*E' + beta*E*H' (or H'*E + E'*H) (SLICOT MB01OE)
// H upper Hessenberg, E upper triangular.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OeUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OeTrans {
    NoTrans, // R := alpha*R + beta*H*E' + beta*E*H'
    Trans,   // R := alpha*R + beta*H'*E + beta*E'*H
}

/// Overwrites the triangle of R.
pub fn mb01oe(
    uplo: Mb01OeUplo,
    trans: Mb01OeTrans,
    n: usize,
    alpha: f64,
    beta: f64,
    r: &mut [f64],
    ldr: usize,
    h: &[f64],
    ldh: usize,
    e: &[f64],
    lde: usize,
) -> i32 {
    if !matches!(uplo, Mb01OeUplo::Upper | Mb01OeUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01OeTrans::NoTrans | Mb01OeTrans::Trans) {
        return -2;
    }
    if ldr < n.max(1) || ldh < n.max(1) || lde < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }
    if beta == 0.0 {
        if alpha != 1.0 {
            for j in 0..n {
                for i in 0..n {
                    if (uplo == Mb01OeUplo::Upper && i <= j) || (uplo == Mb01OeUplo::Lower && i >= j) {
                        r[i + j * ldr] *= alpha;
                    }
                }
            }
        }
        return 0;
    }

    let h_mat = DMatrix::from_fn(n, n, |i, j| {
        if i <= j + 1 {
            h[i + j * ldh]
        } else {
            0.0
        }
    });
    let e_mat = DMatrix::from_fn(n, n, |i, j| if i <= j { e[i + j * lde] } else { 0.0 });
    let update = match trans {
        Mb01OeTrans::NoTrans => &h_mat * e_mat.transpose() + &e_mat * h_mat.transpose(),
        Mb01OeTrans::Trans => h_mat.transpose() * &e_mat + e_mat.transpose() * &h_mat,
    };

    let mut r_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            r_full[(i, j)] = if (uplo == Mb01OeUplo::Upper && i <= j) || (uplo == Mb01OeUplo::Lower && i >= j) {
                r[i + j * ldr]
            } else {
                r[j + i * ldr]
            };
        }
    }
    r_full = alpha * r_full + beta * update;
    for i in 0..n {
        for j in 0..n {
            if (uplo == Mb01OeUplo::Upper && i <= j) || (uplo == Mb01OeUplo::Lower && i >= j) {
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
    fn test_mb01oe_upper_notrans() {
        let n = 2;
        let mut r = vec![0.0, 0.0, 0.0, 0.0];
        let h = [1.0, 0.0, 0.0, 1.0];
        let e = [1.0, 0.0, 0.0, 1.0];
        assert_eq!(
            mb01oe(
                Mb01OeUplo::Upper,
                Mb01OeTrans::NoTrans,
                n,
                0.0,
                1.0,
                &mut r,
                2,
                &h,
                2,
                &e,
                2,
            ),
            0
        );
        assert!((r[0] - 2.0).abs() < 1e-14);
        assert!((r[3] - 2.0).abs() < 1e-14);
    }
}
