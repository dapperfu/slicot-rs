//! TB04BV — Strictly proper part of a proper transfer function matrix (SLICOT TB04BV)
//!
//! Separates the strictly proper part G0 from the constant part D of G.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Order {
    Increasing,
    Decreasing,
}

/// Separates strictly proper part from D. On exit GN/IGN hold G0; D holds the constant part.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1 not proper; 2 null denominator.
pub fn tb04bv(
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
    d: &mut [f64],
    ldd: usize,
    tol: f64,
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
    if d.len() < p * ldd {
        return -11;
    }

    for j in 0..m {
        for i in 0..p {
            let ij = (j * p + i) * md;
            let deg_n = ign[j * ldign + i] as usize;
            let deg_d = igd[j * ldigd + i] as usize;
            if deg_d == 0 && gd.get(ij).map(|&c| c == 0.0).unwrap_or(true) {
                return 2;
            }
            let (lead_n, lead_d) = if order == Order::Decreasing {
                (
                    gn.get(ij).copied().unwrap_or(0.0),
                    gd.get(ij).copied().unwrap_or(0.0),
                )
            } else {
                (
                    gn.get(ij + deg_n).copied().unwrap_or(0.0),
                    gd.get(ij + deg_d).copied().unwrap_or(0.0),
                )
            };
            let tol_use = if tol > 0.0 {
                tol
            } else {
                let norm_n: f64 = gn[ij..ij + md.min(gn.len() - ij)]
                    .iter()
                    .map(|x| x.abs())
                    .fold(0.0_f64, |a, b| a.max(b));
                (deg_n as f64) * f64::EPSILON * norm_n.max(1.0)
            };
            if deg_n < deg_d {
                if d.len() > j * ldd + i {
                    d[j * ldd + i] = 0.0;
                }
                // gn[ij..=ij+deg_n] unchanged when deg_n < deg_d
                for k in (deg_n + 1)..md {
                    if ij + k < gn.len() {
                        gn[ij + k] = 0.0;
                    }
                }
                ign[j * ldign + i] = deg_n as i32;
            } else if deg_n == deg_d {
                let dc = if lead_d.abs() > tol_use {
                    lead_n / lead_d
                } else {
                    return 2;
                };
                if d.len() > j * ldd + i {
                    d[j * ldd + i] = dc;
                }
                for k in 0..=deg_n {
                    let idx = ij + k;
                    if idx < gn.len() && idx < gd.len() {
                        gn[idx] = gn[idx] - dc * gd[idx];
                    }
                }
                let mut new_deg = deg_n as i32;
                while new_deg > 0 {
                    let ck = if order == Order::Decreasing {
                        gn.get(ij + (deg_n - new_deg as usize))
                    } else {
                        gn.get(ij + new_deg as usize)
                    };
                    if ck.map(|&c| c.abs() > tol_use).unwrap_or(false) {
                        break;
                    }
                    new_deg -= 1;
                }
                ign[j * ldign + i] = new_deg;
                for k in (new_deg as usize + 1)..md {
                    if ij + k < gn.len() {
                        gn[ij + k] = 0.0;
                    }
                }
            } else {
                return 1;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb04bv_proper() {
        let mut ign = [2i32, 0i32];
        let igd = [2i32, 1i32];
        let mut gn = vec![7.0, 5.0, 1.0, 1.0, 0.0, 0.0];
        let gd = [6.0, 5.0, 1.0, 2.0, 1.0, 0.0];
        let mut d = [0.0, 0.0, 0.0, 0.0];
        let info = tb04bv(
            Order::Increasing,
            2,
            1,
            3,
            &mut ign,
            2,
            &igd,
            2,
            &mut gn,
            &gd,
            &mut d,
            2,
            0.0,
        );
        assert_eq!(info, 0);
        assert!((d[0] - 1.0).abs() < 1e-10);
    }
}
