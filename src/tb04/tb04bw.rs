//! TB04BW — Sum of a rational matrix and a real matrix (SLICOT TB04BW)
//!
//! Computes G + D where G is P-by-M rational (num/den) and D is real.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Order {
    Increasing,
    Decreasing,
}

/// On exit GN/IGN contain the numerator polynomials of G + D; GD/IGD unchanged.
///
/// # Returns
/// 0 success; < 0 invalid argument.
pub fn tb04bw(
    order: Order,
    p: usize,
    m: usize,
    md: usize,
    ign: &mut [i32],
    ldign: usize,
    igd: &[i32],
    ldigd: usize,
    gn: &mut [f64],
    gd: &[f64],
    d: &[f64],
    ldd: usize,
) -> i32 {
    if p == 0 || m == 0 {
        return 0;
    }
    if ldign < p || ldigd < p || ldd < p {
        return -5;
    }
    let pm_md = p * m * md;
    if gn.len() < pm_md || gd.len() < pm_md {
        return -9;
    }

    for j in 0..m {
        for i in 0..p {
            let ij = (j * p + i) * md;
            let deg_d_ij = igd[j * ldigd + i] as usize;
            let d_ij = d.get(j * ldd + i).copied().unwrap_or(0.0);
            if deg_d_ij == 0 {
                if gn.len() > ij {
                    gn[ij] = gn.get(ij).copied().unwrap_or(0.0) + d_ij;
                }
                ign[j * ldign + i] = 0i32;
                continue;
            }
            if order == Order::Increasing {
                for k in 0..=deg_d_ij.min(md.saturating_sub(1)) {
                    let gd_k = gd.get(ij + k).copied().unwrap_or(0.0);
                    let idx = ij + k;
                    if idx < gn.len() {
                        gn[idx] = gn.get(idx).copied().unwrap_or(0.0) + d_ij * gd_k;
                    }
                }
            } else {
                for k in 0..=deg_d_ij.min(md.saturating_sub(1)) {
                    let gd_k = gd.get(ij + k).copied().unwrap_or(0.0);
                    let idx = ij + k;
                    if idx < gn.len() {
                        gn[idx] = gn.get(idx).copied().unwrap_or(0.0) + d_ij * gd_k;
                    }
                }
            }
            let deg_n = ign[j * ldign + i] as usize;
            ign[j * ldign + i] = deg_n.max(deg_d_ij) as i32;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb04bw_add_zero() {
        let mut ign = [1i32];
        let igd = [1i32];
        let mut gn = [1.0, 1.0];
        let gd = [1.0, 1.0];
        let d = [0.0];
        let info = tb04bw(
            Order::Increasing,
            1,
            1,
            2,
            &mut ign,
            1,
            &igd,
            1,
            &mut gn,
            &gd,
            &d,
            1,
        );
        assert_eq!(info, 0);
        assert!((gn[0] - 1.0).abs() < 1e-10);
        assert!((gn[1] - 1.0).abs() < 1e-10);
    }
}
