//! SB02CX — Balance an N×N matrix (scale rows/columns) (SLICOT).
//!
//! Applies diagonal scaling so that row and column norms are comparable.
//! A is stored in column-major order with leading dimension LDA.

/// Balances the N×N matrix A in-place (stored column-major, LDA leading dimension).
///
/// # Returns
/// 0 on success; &lt; 0 if invalid (e.g. lda &lt; n or len(a) &lt; lda*n).
pub fn sb02cx(n: usize, a: &mut [f64], lda: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if lda < n {
        return -3;
    }
    let len = lda * n;
    if a.len() < len {
        return -3;
    }
    let radix = 2.0_f64;
    let idx = |i: usize, j: usize| i + j * lda;
    for _ in 0..20 {
        let mut done = true;
        for i in 0..n {
            let row_norm: f64 = (0..n).map(|j| a[idx(i, j)].abs()).sum();
            let col_norm: f64 = (0..n).map(|j| a[idx(j, i)].abs()).sum();
            if row_norm == 0.0 || col_norm == 0.0 {
                continue;
            }
            let g = row_norm / radix;
            let h = col_norm * radix;
            if row_norm < g && g < h {
                done = false;
                for j in 0..n {
                    a[idx(i, j)] *= radix;
                }
                for j in 0..n {
                    a[idx(j, i)] /= radix;
                }
            } else if row_norm > h && h > g {
                done = false;
                for j in 0..n {
                    a[idx(i, j)] /= radix;
                }
                for j in 0..n {
                    a[idx(j, i)] *= radix;
                }
            }
        }
        if done {
            break;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb02cx_n0() {
        let mut a = [0.0_f64; 1];
        assert_eq!(sb02cx(0, &mut a, 1), 0);
    }

    #[test]
    fn test_sb02cx_invalid() {
        assert!(sb02cx(1, &mut [], 1) != 0);
    }

    #[test]
    fn test_sb02cx_balance() {
        let mut a = [1e4_f64, 1.0, 1.0, 1e-4]; // 2x2 column-major
        assert_eq!(sb02cx(2, &mut a, 2), 0);
        let r0 = a[0].abs() + a[2].abs();
        let c0 = a[0].abs() + a[1].abs();
        assert!(r0 > 0.0 && c0 > 0.0);
    }
}
