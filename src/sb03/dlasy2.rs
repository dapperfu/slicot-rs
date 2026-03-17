//! DLASY2 — Solve continuous Sylvester op(TL)*X + ISGN*X*op(TR) = SCALE*B.
//! N1, N2 in {1, 2}. Used by SB03OR.


const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Solves op(TL)*X + ISGN*X*op(TR) = SCALE*B. No transpose in this implementation.
/// Returns INFO (0 = ok, 1 = perturbed).
#[allow(clippy::too_many_arguments)]
pub fn dlasy2(
    _ltranl: bool,
    _ltranr: bool,
    isgn: i32,
    n1: usize,
    n2: usize,
    tl: &[f64],
    ldtl: usize,
    tr: &[f64],
    ldtr: usize,
    b: &[f64],
    ldb: usize,
    scale: &mut f64,
    x: &mut [f64],
    ldx: usize,
    xnorm: &mut f64,
) -> i32 {
    let mut info = 0_i32;
    *scale = ONE;
    if n1 == 0 || n2 == 0 {
        *xnorm = ZERO;
        return info;
    }
    let sgn = isgn as f64;
    let eps = f64::EPSILON;
    let smlnum = (ONE / eps).min(f64::MIN_POSITIVE * (ONE / f64::EPSILON));

    if n1 == 1 && n2 == 1 {
        let den = tl[0] + sgn * tr[0];
        let abden = den.abs();
        if abden <= smlnum {
            *scale = ONE / b[0].abs().max(smlnum);
            x[0] = (b[0] * *scale) / smlnum;
            info = 1;
        } else {
            *scale = ONE;
            x[0] = (b[0] * *scale) / den;
        }
        *xnorm = x[0].abs();
        return info;
    }

    if n1 == 1 && n2 == 2 {
        let tr11 = tr[0];
        let tr12 = tr[1];
        let tr21 = if ldtr >= 2 && tr.len() > ldtr {
            tr[ldtr]
        } else {
            ZERO
        };
        let tr22 = if ldtr >= 2 && tr.len() > ldtr + 1 {
            tr[ldtr + 1]
        } else {
            ZERO
        };
        let den = tl[0] + sgn * tr11;
        let den2 = tl[0] + sgn * tr22;
        let b1 = b[0];
        let b2 = if ldb >= 1 && b.len() > 1 {
            b[1]
        } else {
            ZERO
        };
        let abden = den.abs().max(den2.abs());
        if abden <= smlnum {
            *scale = ONE / b1.abs().max(b2.abs()).max(smlnum);
            x[0] = (b1 * *scale) / smlnum;
            if ldx >= 1 && x.len() > 1 {
                x[1] = (b2 * *scale) / smlnum;
            }
            info = 1;
        } else {
            *scale = ONE;
            x[0] = (b1 - sgn * tr21 * (if ldx >= 1 && x.len() > 1 {
                (b2 - sgn * tr12 * b1 / den) / den2
            } else {
                ZERO
            })) / den;
            if ldx >= 1 && x.len() > 1 {
                x[1] = (b2 - sgn * tr12 * x[0]) / den2;
            }
        }
        *xnorm = x[0].abs() + (if x.len() > 1 { x[1].abs() } else { ZERO });
        return info;
    }

    if n1 == 2 && n2 == 1 {
        let tl11 = tl[0];
        let tl12 = tl[1];
        let tl21 = if ldtl >= 2 && tl.len() > ldtl {
            tl[ldtl]
        } else {
            ZERO
        };
        let tl22 = if ldtl >= 2 && tl.len() > ldtl + 1 {
            tl[ldtl + 1]
        } else {
            ZERO
        };
        let den = tl11 + sgn * tr[0];
        let den2 = tl22 + sgn * tr[0];
        let b1 = b[0];
        let b2 = if ldb >= 1 && b.len() > ldb {
            b[ldb]
        } else {
            ZERO
        };
        let abden = den.abs().max(den2.abs());
        if abden <= smlnum {
            *scale = ONE / b1.abs().max(b2.abs()).max(smlnum);
            x[0] = (b1 * *scale) / smlnum;
            if x.len() > ldx {
                x[ldx] = (b2 * *scale) / smlnum;
            }
            info = 1;
        } else {
            *scale = ONE;
            x[0] = (b1 - tl12 * (if x.len() > ldx {
                (b2 - tl21 * b1 / den) / den2
            } else {
                ZERO
            })) / den;
            if x.len() > ldx {
                x[ldx] = (b2 - tl21 * x[0]) / den2;
            }
        }
        *xnorm = x[0].abs() + (if x.len() > ldx { x[ldx].abs() } else { ZERO });
        return info;
    }

    // 2×2: form (I⊗TL + (TR')⊗I)*vec(X) = vec(B), solve 4×4
    let tl11 = tl[0];
    let tl12 = tl[1];
    let tl21 = if ldtl >= 2 { tl[ldtl] } else { ZERO };
    let tl22 = if ldtl >= 2 { tl[ldtl + 1] } else { ZERO };
    let tr11 = tr[0];
    let tr12 = tr[1];
    let tr21 = if ldtr >= 2 && tr.len() > ldtr { tr[ldtr] } else { ZERO };
    let tr22 = if ldtr >= 2 && tr.len() > ldtr + 1 {
        tr[ldtr + 1]
    } else {
        ZERO
    };
    let mut mat = [
        tl11 + sgn * tr11,
        tl21,
        sgn * tr21,
        ZERO,
        tl12,
        tl22 + sgn * tr11,
        ZERO,
        sgn * tr21,
        sgn * tr12,
        ZERO,
        tl11 + sgn * tr22,
        tl21,
        ZERO,
        sgn * tr12,
        tl12,
        tl22 + sgn * tr22,
    ];
    let mut rhs = [b[0], b[ldb], b[1], b[ldb + 1]];
    let smin = mat
        .iter()
        .map(|v| v.abs())
        .fold(f64::MAX, f64::min)
        .max(smlnum)
        * eps;
    for i in 0..4 {
        if mat[i * 4 + i].abs() <= smin {
            mat[i * 4 + i] = smin;
            info = 1;
        }
    }
    for k in 0..3 {
        let mut imax = k;
        let mut vmax = mat[k * 4 + k].abs();
        for i in (k + 1)..4 {
            if mat[i * 4 + k].abs() > vmax {
                vmax = mat[i * 4 + k].abs();
                imax = i;
            }
        }
        if imax != k {
            for j in 0..4 {
                mat.swap(k * 4 + j, imax * 4 + j);
            }
            rhs.swap(k, imax);
        }
        let pivot = mat[k * 4 + k];
        for i in (k + 1)..4 {
            let t = mat[i * 4 + k] / pivot;
            mat[i * 4 + k] = t;
            for j in (k + 1)..4 {
                mat[i * 4 + j] -= t * mat[k * 4 + j];
            }
            rhs[i] -= t * rhs[k];
        }
    }
    for i in (0..4).rev() {
        let mut sum = rhs[i];
        for j in (i + 1)..4 {
            sum -= mat[i * 4 + j] * rhs[j];
        }
        rhs[i] = sum / mat[i * 4 + i];
    }
    x[0] = rhs[0];
    x[ldx] = rhs[1];
    if x.len() > 1 {
        x[1] = rhs[2];
    }
    if x.len() > ldx + 1 {
        x[ldx + 1] = rhs[3];
    }
    *xnorm = rhs[0].abs()
        + rhs[1].abs()
        + rhs[2].abs()
        + rhs[3].abs();
    info
}
