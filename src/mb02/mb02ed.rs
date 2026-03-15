//! MB02ED — Solve T*X = B or X*T = B with positive definite block Toeplitz T (SLICOT).
//!
//! Forms full T from first block row/column, Cholesky, then solve.

use nalgebra::DMatrix;

/// First block row ('R') => solve X*T = B; first block column ('C') => solve T*X = B.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum TypeT {
    Row,
    Col,
}

/// Solves T*X = B (TYPET='C') or X*T = B (TYPET='R'). T from first block row/column in t.
/// B overwritten with solution X. Returns 0, or 1 if T not positive definite.
pub fn mb02ed(
    typet: TypeT,
    k: usize,
    n: usize,
    nrhs: usize,
    t: &[f64],
    ldt: usize,
    b: &mut [f64],
    ldb: usize,
    _dwork: &mut [f64],
) -> i32 {
    if k == 0 || n == 0 || nrhs == 0 {
        return 0;
    }
    let nk = n * k;
    let mut t_full = DMatrix::<f64>::zeros(nk, nk);
    if typet == TypeT::Row {
        for i in 0..n {
            for j in 0..n {
                let blk = if j >= i { j - i } else { i - j };
                for pi in 0..k {
                    for pj in 0..k {
                        t_full[(i * k + pi, j * k + pj)] = if j >= i {
                            t[pi + (blk * k + pj) * ldt]
                        } else {
                            t[pj + (blk * k + pi) * ldt]
                        };
                    }
                }
            }
        }
    } else {
        for i in 0..n {
            for j in 0..n {
                let blk = if j >= i { j - i } else { i - j };
                for pi in 0..k {
                    for pj in 0..k {
                        t_full[(i * k + pi, j * k + pj)] = if j >= i {
                            t[(blk * k + pi) + pj * ldt]
                        } else {
                            t[(blk * k + pj) + pi * ldt]
                        };
                    }
                }
            }
        }
    }
    let ch = match t_full.cholesky() {
        Some(c) => c,
        None => return 1,
    };
    if typet == TypeT::Col {
        let b_mat = DMatrix::<f64>::from_fn(nk, nrhs, |i, j| b[i + j * ldb]);
        let sol = ch.solve(&b_mat);
        for i in 0..nk {
            for j in 0..nrhs {
                b[i + j * ldb] = sol[(i, j)];
            }
        }
    } else {
        let b_mat = DMatrix::<f64>::from_fn(nk, nrhs, |col, row| b[row + col * ldb]);
        let sol = ch.solve(&b_mat);
        for row in 0..nrhs {
            for col in 0..nk {
                b[row + col * ldb] = sol[(col, row)];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb02ed_trivial() {
        let t = vec![];
        let mut b = vec![0.0];
        let mut dwork = vec![0.0];
        assert_eq!(mb02ed(TypeT::Col, 0, 0, 0, &t, 0, &mut b, 0, &mut dwork), 0);
    }

    #[test]
    fn test_mb02ed_1x1_solve() {
        let k = 1;
        let n = 2;
        let t = vec![2.0, 0.5, 0.5, 1.0];
        let mut b = vec![1.0, 1.0, 1.0, 1.0];
        let mut dwork = vec![0.0; 8];
        assert_eq!(
            mb02ed(TypeT::Col, k, n, 2, &t, 1, &mut b, 2, &mut dwork),
            0
        );
    }
}
