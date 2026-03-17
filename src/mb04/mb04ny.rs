//! MB04NY — Apply Householder reflector from the right to C = [A  B].
//!
//! H = I - tau*u*u', u = (1; v). Updates [A B] := [A B]*H.
//! A is M×1, B is M×N. In-line code for order (N+1) < 11.

use crate::mb04::blas::{daxpy, dcopy, dger, dgemv};

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Applies real elementary reflector H from the right to C = [A  B].
///
/// * `m` - number of rows of A and B
/// * `n` - number of columns of B (A has one column)
/// * `v` - vector v in u = (1; v), length n, stride `incv`
/// * `tau` - scalar factor of H; if zero, returns immediately
/// * `a` - M×1 column (lda >= m), updated in place
/// * `b` - M×N matrix (ldb >= m), updated in place
/// * `dwork` - workspace length M (not referenced if order < 11)
#[allow(clippy::too_many_arguments)]
pub fn mb04ny(
    m: usize,
    n: usize,
    v: &[f64],
    incv: i32,
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
    let incv_usize = incv.unsigned_abs() as usize;
    let order = n + 1;
    match order {
        1 => {
            // 1×1 Householder: A := (1-tau)*A
            let t1 = ONE - tau;
            for j in 0..m {
                a[j] *= t1;
            }
        }
        2 => {
            let iv = if incv < 0 {
                ((1 - n as i32) * incv + 1).max(0) as usize
            } else {
                0
            };
            let v1 = v.get(iv).copied().unwrap_or(ZERO);
            let t1 = tau * v1;
            for j in 0..m {
                let sum = a[j] + v1 * b[j];
                a[j] -= sum * tau;
                b[j] -= sum * t1;
            }
        }
        3..=10 => {
            let iv_start = if incv < 0 {
                ((1 - n as i32) * incv + 1).max(0) as usize
            } else {
                0
            };
            let mut vs = [ZERO; 9];
            let mut ts = [ZERO; 9];
            for i in 0..(order - 1) {
                let idx = iv_start + i * incv_usize;
                vs[i] = v.get(idx).copied().unwrap_or(ZERO);
                ts[i] = tau * vs[i];
            }
            for j in 0..m {
                let mut sum = a[j];
                for i in 0..(order - 1) {
                    sum += vs[i] * b[j + i * ldb];
                }
                a[j] -= sum * tau;
                for i in 0..(order - 1) {
                    b[j + i * ldb] -= sum * ts[i];
                }
            }
        }
        _ => {
            // General: w := C*u = A + B*v, C := C - tau*w*u'
            dcopy(m, a, 1, dwork, 1);
            dgemv(false, m, n, ONE, b, ldb, v, incv_usize, ONE, dwork, 1);
            daxpy(m, -tau, dwork, 1, a, 1);
            dger(m, n, -tau, dwork, 1, v, incv_usize, b, ldb);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb04ny_1x1_householder() {
        let mut a = vec![3.0];
        let mut b = vec![0.0];
        let v: Vec<f64> = vec![];
        let mut dwork = vec![0.0; 1];
        mb04ny(1, 0, &v, 1, 0.5, &mut a, 1, &mut b, 1, &mut dwork);
        // (1-tau)*a = 0.5*3 = 1.5
        assert!((a[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_mb04ny_2x2_householder() {
        let mut a = vec![1.0, 0.0];
        let mut b = vec![0.0, 1.0];
        let v = vec![1.0];
        let mut dwork = vec![0.0; 2];
        mb04ny(2, 1, &v, 1, 0.5, &mut a, 2, &mut b, 2, &mut dwork);
        // u = [1, 1]', w = [1, 0] + [0, 1]*1 = [1, 1], sum for row0 = 1, row1 = 1
        // a_new = a - tau*w = [1,0] - 0.5*[1,1] = [0.5, -0.5]
        // b_new = b - tau*w*v' = [0,1] - 0.5*[1,1]*[1] = [0,1] - [0.5,0.5] = [-0.5, 0.5]
        assert!((a[0] - 0.5).abs() < 1e-10);
        assert!((a[1] - (-0.5)).abs() < 1e-10);
        assert!((b[0] - (-0.5)).abs() < 1e-10);
        assert!((b[1] - 0.5).abs() < 1e-10);
    }
}
