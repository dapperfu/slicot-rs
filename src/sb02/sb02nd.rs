//! SB02ND — Optimal state feedback matrix (SLICOT). Stub.

use nalgebra::DMatrix;

/// Stub: validates dimensions, returns 0.
pub fn sb02nd(
    _dico: char,
    _fact: char,
    _uplo: char,
    _jobl: char,
    n: usize,
    m: usize,
    a: &DMatrix<f64>,
    b: &mut DMatrix<f64>,
    _r: &mut DMatrix<f64>,
    _l: Option<&DMatrix<f64>>,
    x: &DMatrix<f64>,
    f: &mut DMatrix<f64>,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -8;
    }
    if b.nrows() != n || b.ncols() != m {
        return -10;
    }
    if x.nrows() != n || x.ncols() != n {
        return -16;
    }
    if f.nrows() != m || f.ncols() != n {
        return -18;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02nd() {
        let a = DMatrix::zeros(2, 2);
        let mut b = DMatrix::zeros(2, 1);
        let mut r = DMatrix::zeros(1, 1);
        let x = DMatrix::zeros(2, 2);
        let mut f = DMatrix::zeros(1, 2);
        assert_eq!(sb02nd('C', 'N', 'U', 'Z', 2, 1, &a, &mut b, &mut r, None, &x, &mut f), 0);
    }
}
