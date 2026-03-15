//! MC01RD — Product of two real polynomials (SLICOT MC01RD)
//
// C(x) = A(x)*B(x). Coefficients in increasing powers.

/// DA, DB degrees; A, B coefficients (increasing); C has length DA+DB+1.
pub fn mc01rd(da: i32, db: i32, a: &[f64], b: &[f64], c: &mut [f64]) -> i32 {
    let (da, db) = (da as usize, db as usize);
    if da.saturating_add(1) > a.len() || db.saturating_add(1) > b.len() {
        return -3;
    }
    if da + db + 1 > c.len() {
        return -5;
    }
    for i in 0..=da + db {
        c[i] = 0.0;
    }
    for i in 0..=da {
        for j in 0..=db {
            c[i + j] += a[i] * b[j];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01rd() {
        let a = [1.0, 1.0];
        let b = [1.0, 1.0];
        let mut c = [0.0; 3];
        assert_eq!(mc01rd(1, 1, &a, &b, &mut c), 0);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }
}
