//! MC01QD — Quotient and remainder of polynomial division (SLICOT MC01QD)
//
// A(x) = B(x)*Q(x) + R(x), R degree < degree(B). Coefficients in increasing powers.

/// On exit: RQ[0..db] = R (if da >= db), RQ[db..da+1] = Q. DB may be reduced if leading B zeros.
pub fn mc01qd(
    da: i32,
    db: &mut i32,
    a: &[f64],
    b: &[f64],
    rq: &mut [f64],
    iwarn: &mut i32,
) -> i32 {
    let mut da = da as i32;
    let mut db_val = *db as usize;
    *iwarn = 0;
    if da < -1 {
        return -1;
    }
    if *db < 0 {
        return -2;
    }
    let db_in = *db as usize;
    if db_in.saturating_add(1) > b.len() {
        return -4;
    }
    if da >= 0 && da as usize + 1 > a.len() {
        return -3;
    }
    if (da as usize).saturating_add(1) > rq.len() {
        return -5;
    }
    // Reduce DB if leading coefficients of B are zero
    while db_val > 0 && b[db_val].abs() < 1e-15 {
        db_val -= 1;
        *iwarn += 1;
    }
    *db = db_val as i32;
    if da < *db {
        if da >= 0 {
            rq[..=da as usize].copy_from_slice(&a[..=da as usize]);
        }
        return 0;
    }
    if db_val == 0 && b[0].abs() < 1e-15 {
        return 1;
    }
    let mut work: Vec<f64> = a[..=da as usize].to_vec();
    let b_lead = b[db_val];
    let dq = (da - *db) as usize;
    for k in (0..=dq).rev() {
        let qk = work[db_val + k] / b_lead;
        for j in 0..=db_val {
            work[k + j] -= qk * b[j];
        }
        rq[db_val + k] = qk;
    }
    for i in 0..db_val {
        rq[i] = work[i];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01qd_example() {
        // A = 2+2x -x^2 +2x^3 +x^4 (da=4), B = 1 -x +x^2 (db=2)
        let a = [2.0, 2.0, -1.0, 2.0, 1.0];
        let b = [1.0, -1.0, 1.0];
        let mut db = 2i32;
        let mut rq = [0.0; 5];
        let mut iwarn = 0;
        assert_eq!(mc01qd(4, &mut db, &a, &b, &mut rq, &mut iwarn), 0);
        // R = 1+0x, Q = 1+3x+x^2
        assert!((rq[0] - 1.0).abs() < 1e-6);
        assert!((rq[1] - 0.0).abs() < 1e-6);
        assert!((rq[2] - 1.0).abs() < 1e-6);
        assert!((rq[3] - 3.0).abs() < 1e-6);
        assert!((rq[4] - 1.0).abs() < 1e-6);
    }
}
