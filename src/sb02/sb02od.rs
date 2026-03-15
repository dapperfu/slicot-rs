//! SB02OD — Solution of continuous/discrete-time algebraic Riccati equation (general) (SLICOT). Stub.

use nalgebra::DMatrix;

/// Stub: validates N, M, returns 0.
pub fn sb02od(
    _dico: char,
    _jobb: char,
    _fact: &str,
    _uplo: char,
    _jobl: char,
    _sort: char,
    n: usize,
    m: usize,
    _p: usize,
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    _c: Option<&DMatrix<f64>>,
    _r: &DMatrix<f64>,
    _l: Option<&DMatrix<f64>>,
    x: &mut DMatrix<f64>,
    _tol: f64,
) -> i32 {
    if a.nrows() != n || a.ncols() != n {
        return -10;
    }
    if b.nrows() != n || b.ncols() != m {
        return -12;
    }
    if x.nrows() != n || x.ncols() != n {
        return -18;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02od() {
        let a = DMatrix::zeros(1, 1);
        let b = DMatrix::zeros(1, 1);
        let r = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb02od('C', 'B', "Both", 'U', 'Z', 'S', 1, 1, 1, &a, &b, None, &r, None, &mut x, 0.0), 0);
    }
}
