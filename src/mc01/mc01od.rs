//! MC01OD — Coefficients of complex polynomial from zeros (SLICOT MC01OD)
//
// P(x) = (x-r1)(x-r2)...(x-rK). Output REP, IMP in increasing powers of x.

/// Zeros in REZ, IMZ. REP, IMP hold real and imaginary parts of coefficients (increasing powers).
pub fn mc01od(
    k: i32,
    rez: &[f64],
    imz: &[f64],
    rep: &mut [f64],
    imp: &mut [f64],
    dwork: &mut [f64],
) -> i32 {
    let k = k as usize;
    if k > rez.len() || k > imz.len() {
        return -2;
    }
    if k + 1 > rep.len() || k + 1 > imp.len() {
        return -4;
    }
    if k > 0 && 2 * k + 2 > dwork.len() {
        return -6;
    }
    if k == 0 {
        rep[0] = 1.0;
        imp[0] = 0.0;
        return 0;
    }
    rep[0] = -rez[0];
    imp[0] = -imz[0];
    rep[1] = 1.0;
    imp[1] = 0.0;
    let mut deg = 1;
    for i in 1..k {
        let (zr, zi) = (rez[i], imz[i]);
        dwork[0] = rep[0];
        dwork[1] = imp[0];
        for j in 1..=deg {
            dwork[2 * j] = rep[j];
            dwork[2 * j + 1] = imp[j];
        }
        rep[0] = -(zr * dwork[0] - zi * dwork[1]);
        imp[0] = -(zr * dwork[1] + zi * dwork[0]);
        for j in 1..=deg {
            rep[j] = dwork[2 * (j - 1)] - (zr * dwork[2 * j] - zi * dwork[2 * j + 1]);
            imp[j] = dwork[2 * (j - 1) + 1] - (zr * dwork[2 * j + 1] + zi * dwork[2 * j]);
        }
        deg += 1;
        rep[deg] = dwork[2 * (deg - 1)];
        imp[deg] = dwork[2 * (deg - 1) + 1];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01od_one_zero() {
        let mut rep = [0.0; 2];
        let mut imp = [0.0; 2];
        let mut dwork = [0.0; 4];
        assert_eq!(mc01od(1, &[1.0], &[0.0], &mut rep, &mut imp, &mut dwork), 0);
        assert!((rep[0] - (-1.0)).abs() < 1e-10);
        assert!((rep[1] - 1.0).abs() < 1e-10);
    }
}
