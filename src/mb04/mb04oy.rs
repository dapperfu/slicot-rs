//! MB04OY — Apply Householder reflector from the left to C = [A; B].
//!
//! H = I - tau*u*u', u = (1; v). Updates [A; B] := H*[A; B].
//! A is 1×N, B is M×N. In-line code for order (M+1) < 11.

use crate::mb04::blas::{daxpy, dcopy, dger, dgemv};

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Applies real elementary reflector H from the left to C = [A; B].
///
/// * `m` - number of rows of B (A has one row)
/// * `n` - number of columns of A and B
/// * `v` - vector v in u = (1; v), length m, contiguous
/// * `tau` - scalar factor of H; if zero, returns immediately
/// * `a` - 1×N row (lda >= 1), updated in place
/// * `b` - M×N matrix (ldb >= m), updated in place
/// * `dwork` - workspace length N (not referenced if order < 11)
#[allow(clippy::too_many_arguments)]
pub fn mb04oy(
    m: usize,
    n: usize,
    v: &[f64],
    tau: f64,
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    dwork: &mut [f64],
) {
    if tau == ZERO {
        return;
    }
    let order = m + 1;
    match order {
        1 => {
            let t1 = ONE - tau;
            for j in 0..n {
                a[j * lda] *= t1;
            }
        }
        2 => {
            let v1 = v.get(0).copied().unwrap_or(ZERO);
            let t1 = tau * v1;
            for j in 0..n {
                let sum = a[j * lda] + v1 * b[j * ldb];
                a[j * lda] -= sum * tau;
                b[j * ldb] -= sum * t1;
            }
        }
        3..=10 => {
            let mut vs = [ZERO; 9];
            let mut ts = [ZERO; 9];
            for i in 0..(order - 1) {
                vs[i] = v.get(i).copied().unwrap_or(ZERO);
                ts[i] = tau * vs[i];
            }
            for j in 0..n {
                let mut sum = a[j * lda];
                for i in 0..(order - 1) {
                    sum += vs[i] * b[i + j * ldb];
                }
                a[j * lda] -= sum * tau;
                for i in 0..(order - 1) {
                    b[i + j * ldb] -= sum * ts[i];
                }
            }
        }
        _ => {
            // w := C'*u = A' + B'*v, C := C - tau*u*w'
            dcopy(n, a, lda, dwork, 1);
            dgemv(true, m, n, ONE, b, ldb, v, 1, ONE, dwork, 1);
            daxpy(n, -tau, dwork, 1, a, lda);
            dger(m, n, -tau, v, 1, dwork, 1, b, ldb);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb04oy_1x1_householder() {
        let mut a = vec![3.0];
        let mut b = vec![0.0; 0];
        let v: Vec<f64> = vec![];
        let mut dwork = vec![0.0; 1];
        mb04oy(0, 1, &v, 0.5, &mut a, 1, &mut b, 1, &mut dwork);
        assert!((a[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_mb04oy_2x2_householder() {
        let mut a = vec![1.0, 0.0];
        let mut b = vec![0.0, 1.0];
        let v = vec![1.0];
        let mut dwork = vec![0.0; 2];
        mb04oy(1, 2, &v, 0.5, &mut a, 1, &mut b, 1, &mut dwork);
        // u = [1, 1]', w' = a + v'*b = [1,0] + 1*[0,1] = [1,1]
        // a_new = a - tau*w' = [1,0] - 0.5*[1,1] = [0.5, -0.5]
        // b_new = b - tau*v*w' = [0,1] - 0.5*[1,1] = [-0.5, 0.5]
        assert!((a[0] - 0.5).abs() < 1e-10);
        assert!((a[1] - (-0.5)).abs() < 1e-10);
        assert!((b[0] - (-0.5)).abs() < 1e-10);
        assert!((b[1] - 0.5).abs() < 1e-10);
    }
}
