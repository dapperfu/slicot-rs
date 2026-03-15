//! MB01OO — P = op(H)*X*op(E)' (SLICOT MB01OO)
// H upper Hessenberg, X symmetric, E upper triangular.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OoUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01OoTrans {
    NoTrans, // P = H*X*E'
    Trans,   // P' = E'*X*H, so P = (E'*X*H)' = H'*X*E
}

/// Computes P = op(H)*X*op(E)' (NoTrans) or P' = E'*X*H (Trans). Overwrites P.
pub fn mb01oo(
    uplo: Mb01OoUplo,
    trans: Mb01OoTrans,
    n: usize,
    h: &[f64],
    ldh: usize,
    x: &[f64],
    ldx: usize,
    e: &[f64],
    lde: usize,
    p: &mut [f64],
    ldp: usize,
) -> i32 {
    if !matches!(uplo, Mb01OoUplo::Upper | Mb01OoUplo::Lower) {
        return -1;
    }
    if !matches!(trans, Mb01OoTrans::NoTrans | Mb01OoTrans::Trans) {
        return -2;
    }
    if ldh < n.max(1) || ldx < n.max(1) || lde < n.max(1) || ldp < n.max(1) {
        return -7;
    }
    if n == 0 {
        return 0;
    }

    let h_mat = DMatrix::from_fn(n, n, |i, j| if i <= j + 1 { h[i + j * ldh] } else { 0.0 });
    let e_mat = DMatrix::from_fn(n, n, |i, j| if i <= j { e[i + j * lde] } else { 0.0 });
    let mut x_full = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            x_full[(i, j)] = if (uplo == Mb01OoUplo::Upper && i <= j) || (uplo == Mb01OoUplo::Lower && i >= j) {
                x[i + j * ldx]
            } else {
                x[j + i * ldx]
            };
        }
    }
    let result = match trans {
        Mb01OoTrans::NoTrans => &h_mat * &x_full * e_mat.transpose(),
        Mb01OoTrans::Trans => (e_mat.transpose() * &x_full * &h_mat).transpose(),
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
    fn test_mb01oo_upper_notrans() {
        let n = 2;
        let h = [1.0, 0.0, 0.0, 1.0];
        let x = [1.0, 0.0, 0.0, 1.0];
        let e = [1.0, 0.0, 0.0, 1.0];
        let mut p = vec![0.0; 4];
        assert_eq!(
            mb01oo(Mb01OoUplo::Upper, Mb01OoTrans::NoTrans, n, &h, 2, &x, 2, &e, 2, &mut p, 2),
            0
        );
        assert!((p[0] - 1.0).abs() < 1e-14);
        assert!((p[3] - 1.0).abs() < 1e-14);
    }
}
