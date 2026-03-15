//! MC01XD — Composition of real polynomials (SLICOT MC01XD)
//
// C(x) = A(B(x)). Coefficients in increasing powers.

/// DA degree of A, DB degree of B. C has degree DA*DB, length DA*DB+1.
pub fn mc01xd(da: i32, db: i32, a: &[f64], b: &[f64], c: &mut [f64], dwork: &mut [f64]) -> i32 {
    let (da, db) = (da as usize, db as usize);
    if da + 1 > a.len() || db + 1 > b.len() {
        return -3;
    }
    let dc = da * db;
    if dc + 1 > c.len() || dc + 1 > dwork.len() {
        return -5;
    }
    for i in 0..=dc {
        c[i] = 0.0;
    }
    c[0] = a[0];
    if da == 0 {
        return 0;
    }
    let mut pow: Vec<f64> = vec![0.0; dc + 1];
    pow[0] = 1.0;
    let mut pow_deg = 0_usize;
    for i in 1..=da {
        for j in 0..=pow_deg + db {
            let mut s = 0.0;
            for k in 0..=pow_deg.min(j) {
                if j - k <= db {
                    s += pow[k] * b[j - k];
                }
            }
            dwork[j] = s;
        }
        pow_deg += db;
        for j in 0..=pow_deg {
            pow[j] = dwork[j];
        }
        for j in 0..=pow_deg.min(dc) {
            c[j] += a[i] * pow[j];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01xd_linear() {
        let a = [0.0, 1.0];
        let b = [1.0, 1.0];
        let mut c = [0.0; 2];
        let mut dwork = [0.0; 2];
        assert_eq!(mc01xd(1, 1, &a, &b, &mut c, &mut dwork), 0);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
    }
}
