//! MB01TD — B := A*B where A and B are upper quasi-triangular (SLICOT MB01TD)
// Simplified: upper triangular A and B; result overwrites B.

/// Computes B := A*B. A and B are N×N upper triangular (elements below diagonal not referenced).
pub fn mb01td(
    n: usize,
    a: &[f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    dwork: &mut [f64],
) -> i32 {
    if lda < n.max(1) || ldb < n.max(1) {
        return -3;
    }
    if dwork.len() < n.saturating_sub(1) && n > 1 {
        return -6;
    }
    if n == 0 {
        return 0;
    }
    if n == 1 {
        b[0] = a[0] * b[0];
        return 0;
    }
    // Column j of result = A * B(:,j). Compute from j=0 to n-1 using temp for column.
    for j in 0..n {
        let mut col: Vec<f64> = (0..n).map(|_| 0.0).collect();
        for i in 0..=j {
            for k in i..=j {
                col[i] += a[i + k * lda] * b[k + j * ldb];
            }
        }
        for i in 0..=j {
            b[i + j * ldb] = col[i];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01td_identity() {
        let n = 2;
        let a = [1.0, 0.0, 0.0, 1.0];
        let mut b = [1.0, 0.0, 2.0, 3.0];
        let mut dwork = vec![0.0; n];
        assert_eq!(mb01td(n, &a, 2, &mut b, 2, &mut dwork), 0);
        assert!((b[0] - 1.0).abs() < 1e-14);
        assert!((b[2] - 2.0).abs() < 1e-14);
        assert!((b[3] - 3.0).abs() < 1e-14);
    }
}
