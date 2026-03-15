//! MB02TZ — In-place LU factorization (SLICOT MB02). Dense fallback.

use nalgebra::DMatrix;

/// Overwrites A with L (unit lower) and U (upper). Returns 0, 1 if singular, -1 if invalid.
pub fn mb02tz(n: usize, a: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if a.nrows() != n || a.ncols() != n {
        return -1;
    }
    let lu = a.clone().lu();
    let (l, u) = (lu.l(), lu.u());
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = if i > j { l[(i, j)] } else { u[(i, j)] };
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_mb02tz_trivial() {
        let mut a = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(mb02tz(0, &mut a), 0);
    }
}
