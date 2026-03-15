//! MB02CU — In-place Cholesky (upper) of s.p.d. matrix (SLICOT MB02). Dense fallback.

use nalgebra::DMatrix;

/// Overwrites upper triangle of A with Cholesky R (A = R'*R). Returns 0, 1 if not s.p.d., -1 if invalid.
pub fn mb02cu(n: usize, a: &mut DMatrix<f64>, _lda: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() < n || a.ncols() < n {
        return -1;
    }
    let view = a.view((0, 0), (n, n));
    let sub = view.into_owned();
    let ch = match sub.cholesky() {
        Some(c) => c,
        None => return 1,
    };
    let l = ch.l();
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = if j >= i { l[(j, i)] } else { 0.0 };
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb02cu_trivial() {
        let mut a = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb02cu(0, &mut a, 0), 0);
    }
}
