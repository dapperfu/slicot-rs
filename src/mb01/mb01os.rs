//! MB01OS — P = H*X or P = X*H (SLICOT MB01OS)
// H upper Hessenberg, X symmetric.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OsUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OsTrans {
    NoTrans, // P = H*X
    Trans,   // P = X*H
}

/// P = H*X (NoTrans) or P = X*H (Trans). Overwrites P.
pub fn mb01os(
    uplo: Mb01OsUplo,
    trans: Mb01OsTrans,
    n: usize,
    h: &[f64],
    ldh: usize,
    x: &[f64],
    ldx: usize,
    p: &mut [f64],
    ldp: usize,
) -> i32 {
    if !matches!(uplo, Mb01OsUplo::Upper | Mb01OsUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01OsTrans::NoTrans | Mb01OsTrans::Trans) {
        return -2;
    }
    if ldh < n.max(1) || ldx < n.max(1) || ldp < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }

    let h_mat = DMatrix::from_fn(n, n, |i, j| if i <= j + 1 { h[i + j * ldh] } else { 0.0 });
    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            x_full[(i, j)] = if (uplo == Mb01OsUplo::Upper && i <= j) || (uplo == Mb01OsUplo::Lower && i >= j) {
                x[i + j * ldx]
            } else {
                x[j + i * ldx]
            };
        }
    }
    let result = match trans {
        Mb01OsTrans::NoTrans => &h_mat * &x_full,
        Mb01OsTrans::Trans => &x_full * &h_mat,
    };
    for j in 0..n {
        for i in 0..n {
            p[i + j * ldp] = result[(i, j)];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01os_upper_notrans() {
        let n = 2;
        let h = [1.0, 0.0, 0.0, 1.0];
        let x = [1.0, 0.0, 0.0, 1.0];
        let mut p = vec![0.0; 4];
        assert_eq!(mb01os(Mb01OsUplo::Upper, Mb01OsTrans::NoTrans, n, &h, 2, &x, 2, &mut p, 2), 0);
        assert!((p[0] - 1.0).abs() < 1e-14);
        assert!((p[3] - 1.0).abs() < 1e-14);
    }
}
