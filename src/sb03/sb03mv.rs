//! SB03MV — Solve 2×2 discrete Lyapunov: op(T)'*X*op(T) - X = scale*B (symmetric X).
//! T and B are 2×2; X overwrites B. Uses 3×3 equivalent system with complete pivoting.

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const FOUR: f64 = 4.0;

/// Solves op(T)'*X*op(T) - X = scale*B for symmetric 2×2 X.
/// `ltran`: false = op(T)=T, true = op(T)=T'.
/// `lupper`: true = use upper triangle of B/X, false = lower.
/// Returns (scale, xnorm, info). info=1 if T has almost reciprocal eigenvalues (perturbed).
#[allow(clippy::too_many_arguments)]
pub fn sb03mv(
    ltran: bool,
    lupper: bool,
    t: &[f64],
    ldt: usize,
    b: &[f64],
    ldb: usize,
    scale: &mut f64,
    x: &mut [f64],
    ldx: usize,
    xnorm: &mut f64,
    info: &mut i32,
) {
    *info = 0;
    let eps = f64::EPSILON;
    let smlnum = (f64::MIN_POSITIVE / eps).min(1e-10);

    let t11 = t[0];
    let t12 = t[if ltran { 1 } else { ldt }];
    let t21 = t[if ltran { ldt } else { 1 }];
    let t22 = t[1 + ldt];

    let smin = (t11.abs())
        .max(t12.abs())
        .max(t21.abs())
        .max(t22.abs())
        .mul_add(eps, smlnum);

    // Build 3×3 system for discrete: T'*X*T - X = B, vec(X)=[x11,x12,x22]
    let mut t9 = [[0.0f64; 3]; 3];
    t9[0][0] = t11 * t11 - ONE;
    t9[1][1] = t11 * t22 + t12 * t21 - ONE;
    t9[2][2] = t22 * t22 - ONE;
    if ltran {
        t9[0][1] = t11 * t12 + t11 * t12;
        t9[0][2] = t12 * t12;
        t9[1][0] = t11 * t21;
        t9[1][2] = t12 * t22;
        t9[2][0] = t21 * t21;
        t9[2][1] = t21 * t22 + t21 * t22;
    } else {
        t9[0][1] = t11 * t21 + t11 * t21;
        t9[0][2] = t21 * t21;
        t9[1][0] = t11 * t12;
        t9[1][2] = t21 * t22;
        t9[2][0] = t12 * t12;
        t9[2][1] = t12 * t22 + t12 * t22;
    }

    let mut btmp = [
        b[0],
        if lupper { b[ldb] } else { b[1] }, // B(1,2) upper vs B(2,1) lower
        b[1 + ldb],
    ];

    let mut jpiv = [0usize, 1, 2];
    for i in 0..2 {
        let mut xmax = 0.0f64;
        let mut ipsv = i;
        let mut jpsv = i;
        for ip in i..3 {
            for jp in i..3 {
                let a = t9[ip][jp].abs();
                if a >= xmax {
                    xmax = a;
                    ipsv = ip;
                    jpsv = jp;
                }
            }
        }
        if ipsv != i {
            t9.swap(ipsv, i);
            btmp.swap(ipsv, i);
        }
        if jpsv != i {
            for r in 0..3 {
                t9[r].swap(jpsv, i);
            }
        }
        jpiv[i] = jpsv;
        if t9[i][i].abs() < smin {
            *info = 1;
            t9[i][i] = smin;
        }
        for j in (i + 1)..3 {
            let mult = t9[j][i] / t9[i][i];
            t9[j][i] = mult;
            btmp[j] -= mult * btmp[i];
            for k in (i + 1)..3 {
                t9[j][k] -= mult * t9[i][k];
            }
        }
    }
    if t9[2][2].abs() < smin {
        *info = 1;
        t9[2][2] = smin;
    }
    *scale = ONE;
    let bmax = btmp[0].abs().max(btmp[1].abs()).max(btmp[2].abs());
    if (FOUR * smlnum) * bmax > t9[0][0].abs().min(t9[1][1].abs()).min(t9[2][2].abs()) {
        *scale = (ONE / FOUR) / bmax;
        btmp[0] *= *scale;
        btmp[1] *= *scale;
        btmp[2] *= *scale;
    }
    let mut tmp = [0.0f64; 3];
    for &k in [2, 1, 0].iter() {
        let temp = ONE / t9[k][k];
        tmp[k] = btmp[k] * temp;
        for j in (k + 1)..3 {
            tmp[k] -= (temp * t9[k][j]) * tmp[j];
        }
    }
    for i in 0..2 {
        let idx = 2 - i;
        if jpiv[idx] != idx {
            tmp.swap(idx, jpiv[idx]);
        }
    }
    x[0] = tmp[0];
    if lupper {
        x[1] = tmp[1];
    } else {
        x[ldx] = tmp[1];
    }
    x[1 + ldx] = tmp[2];
    *xnorm = (tmp[0].abs() + tmp[1].abs()).max(tmp[1].abs() + tmp[2].abs());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03mv_1x1_equiv() {
        // T = 0.5*I => T'*X*T - X = -0.75*X = B, X = -B/0.75
        let t = [0.5, 0.0, 0.0, 0.5];
        let b = [1.0, 0.0, 0.0, 1.0];
        let mut x = [0.0; 4];
        let mut scale = 0.0;
        let mut xnorm = 0.0;
        let mut info = 0;
        sb03mv(false, true, &t, 2, &b, 2, &mut scale, &mut x, 2, &mut xnorm, &mut info);
        assert_eq!(info, 0);
        let expect = -1.0 / 0.75;
        assert!((x[0] - expect).abs() < 1e-10);
        assert!((x[3] - expect).abs() < 1e-10);
    }
}
