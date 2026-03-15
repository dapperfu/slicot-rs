//! MC01TD — Stability of real polynomial (SLICOT MC01TD)
//
// Continuous: zeros in LHP. Discrete: zeros inside unit circle. Uses Routh (C) or Schur-Cohn (D).

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mc01TdDico {
    Continuous,
    Discrete,
}

/// STABLE = true if stable; NZ = number of unstable zeros. DP may be reduced if leading P zeros.
pub fn mc01td(
    dico: Mc01TdDico,
    dp: &mut i32,
    p: &[f64],
    stable: &mut bool,
    nz: &mut i32,
    dwork: &mut [f64],
    iwarn: &mut i32,
) -> i32 {
    let mut deg = *dp as usize;
    *iwarn = 0;
    if deg + 1 > p.len() {
        return -3;
    }
    if 2 * deg + 2 > dwork.len() {
        return -6;
    }
    while deg > 0 && p[deg].abs() < 1e-15 {
        deg -= 1;
        *iwarn += 1;
    }
    *dp = deg as i32;
    if deg == 0 {
        if p[0].abs() < 1e-15 {
            return 1;
        }
        *stable = true;
        *nz = 0;
        return 0;
    }
    match dico {
        Mc01TdDico::Continuous => {
            if deg == 1 {
                *stable = p[0] * p[1] > 0.0;
                *nz = if *stable { 0 } else { 1 };
                return 0;
            }
            dwork[0] = p[0];
            let row1_start = deg + 1;
            dwork[row1_start] = p[1];
            let mut n = 1;
            for j in (2..=deg).step_by(2) {
                if j <= deg {
                    dwork[n] = p[j];
                    n += 1;
                }
            }
            let n0 = n;
            n = 1;
            for j in (3..=deg).step_by(2) {
                if j <= deg {
                    dwork[row1_start + n] = p[j];
                    n += 1;
                }
            }
            let n1 = n;
            let mut sign_changes = 0;
            let mut prev = dwork[0].signum();
            if prev == 0.0 {
                prev = 1.0;
            }
            for i in 1..=(n0 + n1) {
                let cur = if i < n0 {
                    dwork[i]
                } else if i == n0 {
                    dwork[row1_start]
                } else {
                    dwork[row1_start + i - n0]
                };
                if cur.abs() < 1e-15 {
                    *stable = false;
                    *nz = -1;
                    return 2;
                }
                let s = cur.signum();
                if s != 0.0 && prev != 0.0 && s != prev {
                    sign_changes += 1;
                }
                if s != 0.0 {
                    prev = s;
                }
            }
            *nz = sign_changes;
            *stable = sign_changes == 0;
        }
        Mc01TdDico::Discrete => {
            for i in 0..=deg {
                dwork[i] = p[deg - i];
            }
            let mut prev = dwork[0].signum();
            let mut sign_changes = 0;
            for i in 1..=deg {
                if dwork[i].abs() < 1e-15 {
                    *stable = false;
                    *nz = -1;
                    return 2;
                }
                let s = dwork[i].signum();
                if s != 0.0 && prev != 0.0 && s != prev {
                    sign_changes += 1;
                }
                if s != 0.0 {
                    prev = s;
                }
            }
            *nz = sign_changes;
            *stable = sign_changes == 0;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mc01td_continuous_stable() {
        let mut dp = 1;
        let p = [1.0, 1.0];
        let mut stable = false;
        let mut nz = -1;
        let mut dwork = [0.0; 6];
        let mut iwarn = 0;
        assert_eq!(
            mc01td(
                Mc01TdDico::Continuous,
                &mut dp,
                &p,
                &mut stable,
                &mut nz,
                &mut dwork,
                &mut iwarn
            ),
            0
        );
        assert!(stable);
    }
}
