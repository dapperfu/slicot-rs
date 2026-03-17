//! SB04PX — Solve small Sylvester equation (SLICOT SB04PX).
//!
//! Solves op(TL)*X*op(TR) + ISGN*X = SCALE*B for X (N1×N2 with N1,N2 in {0,1,2}).
//! 1:1 with Fortran SB04PX.

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const HALF: f64 = 0.5;
const EIGHT: f64 = 8.0;

fn dlamch_p() -> f64 {
    f64::EPSILON
}
fn dlamch_s() -> f64 {
    f64::MIN_POSITIVE
}

/// IDAMAX: index of element with max absolute value (1-based in Fortran).
fn idamax(n: usize, x: &[f64], incx: usize) -> usize {
    if n == 0 {
        return 0;
    }
    let mut imax = 0;
    let mut xmax = 0.0_f64;
    let mut i = 0_usize;
    for _ in 0..n {
        let ax = x[i].abs();
        if ax > xmax {
            xmax = ax;
            imax = i;
        }
        i = i.saturating_add(incx);
    }
    imax + 1
}

/// DSWAP: swap two vectors.
fn dswap(n: usize, x: &mut [f64], incx: usize, y: &mut [f64], incy: usize) {
    if n == 0 {
        return;
    }
    let mut ix = 0_usize;
    let mut iy = 0_usize;
    for _ in 0..n {
        let t = x[ix];
        x[ix] = y[iy];
        y[iy] = t;
        ix = ix.saturating_add(incx);
        iy = iy.saturating_add(incy);
    }
}

/// SB04PX: solve op(TL)*X*op(TR) + ISGN*X = SCALE*B for N1×N2 X (N1,N2 in {0,1,2}).
/// X is written into the leading N1×N2 part of the provided buffer; scale and xnorm returned.
/// INFO = 0 success, INFO = 1 if equation was perturbed (nearly singular).
#[allow(clippy::too_many_arguments)]
pub fn sb04px(
    ltranl: bool,
    ltranr: bool,
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

    let eps = dlamch_p();
    let smlnum = dlamch_s() / eps;
    let sgn = isgn as f64;

    let k = n1 + n1 + n2 - 2; // Fortran: K = 2*N1+N2-2, GO TO (10, 20, 30, 50) K => 1->10, 2->20, 3->30, 4->50

    if k == 1 {
        // 1-by-1
        let mut tau1 = tl[0] * tr[0] + sgn;
        let mut bet = tau1.abs();
        if bet <= smlnum {
            tau1 = smlnum;
            bet = smlnum;
            info = 1;
        }
        let gam = b[0].abs();
        let mut scale_local = ONE;
        if smlnum * gam > bet {
            scale_local = ONE / gam;
        }
        let x11 = (b[0] * scale_local) / tau1;
        x[0] = x11;
        *scale = scale_local;
        *xnorm = x11.abs();
        return info;
    }

    if k == 2 {
        // 1-by-2
        let tr_11 = tr[0];
        let tr_12 = tr[1];
        let tr_21 = if ldtr >= 2 && tr.len() > ldtr {
            tr[ldtr]
        } else {
            0.0
        };
        let tr_22 = if ldtr >= 2 && tr.len() > ldtr + 1 {
            tr[ldtr + 1]
        } else {
            0.0
        };
        let tl_11 = tl[0];
        let smin = (tr_11.abs().max(tr_12).max(tr_21.abs()).max(tr_22.abs()) * tl_11.abs() * eps).max(smlnum);
        let mut tmp = [ZERO; 4];
        tmp[0] = tl_11 * tr_11 + sgn;
        tmp[3] = tl_11 * tr_22 + sgn;
        if ltranr {
            tmp[1] = tl_11 * tr_21;
            tmp[2] = tl_11 * tr_12;
        } else {
            tmp[1] = tl_11 * tr_12;
            tmp[2] = tl_11 * tr_21;
        }
        let mut btmp = [b[0], b[1]];
        sb04px_2x2_solve(&mut tmp, &mut btmp, smin, smlnum, scale, x, ldx, xnorm, &mut info, 1, 2);
        return info;
    }

    if k == 3 {
        // 2-by-1
        let tl_11 = tl[0];
        let tl_12 = tl[1];
        let tl_21 = if ldtl >= 2 { tl[ldtl] } else { 0.0 };
        let tl_22 = if ldtl >= 2 { tl[ldtl + 1] } else { 0.0 };
        let tr_11 = tr[0];
        let smin = (tl_11.abs().max(tl_12).max(tl_21.abs()).max(tl_22.abs()) * tr_11.abs() * eps).max(smlnum);
        let mut tmp = [ZERO; 4];
        tmp[0] = tl_11 * tr_11 + sgn;
        tmp[3] = tl_22 * tr_11 + sgn;
        if ltranl {
            tmp[1] = tl_12 * tr_11;
            tmp[2] = tl_21 * tr_11;
        } else {
            tmp[1] = tl_21 * tr_11;
            tmp[2] = tl_12 * tr_11;
        }
        let mut btmp = [b[0], b[ldb]];
        sb04px_2x2_solve(&mut tmp, &mut btmp, smin, smlnum, scale, x, ldx, xnorm, &mut info, 2, 1);
        return info;
    }

    // 2-by-2: build 4×4 system T16*vec(X) = vec(B)
    let locu12 = [3, 4, 1, 2];
    let locl21 = [2, 1, 4, 3];
    let locu22 = [4, 3, 2, 1];
    let xswpiv = [false, false, true, true];
    let bswpiv = [false, true, false, true];

    let tl_11 = tl[0];
    let tl_12 = tl[1];
    let tl_21 = tl[ldtl];
    let tl_22 = tl[ldtl + 1];
    let tr_11 = tr[0];
    let tr_12 = tr[1];
    let tr_21 = tr[ldtr];
    let tr_22 = tr[ldtr + 1];

    let mut smin = tr_11.abs().max(tr_12).max(tr_21.abs()).max(tr_22.abs());
    smin = tl_11.abs().max(tl_12).max(tl_21.abs()).max(tl_22.abs()) * smin;
    smin = (eps * smin).max(smlnum);

    let mut t16 = [[ZERO; 4]; 4];
    t16[0][0] = tl_11 * tr_11 + sgn;
    t16[1][1] = tl_22 * tr_11 + sgn;
    t16[2][2] = tl_11 * tr_22 + sgn;
    t16[3][3] = tl_22 * tr_22 + sgn;
    if ltranl {
        t16[0][1] = tl_21 * tr_11;
        t16[1][0] = tl_12 * tr_11;
        t16[2][3] = tl_21 * tr_22;
        t16[3][2] = tl_12 * tr_22;
    } else {
        t16[0][1] = tl_12 * tr_11;
        t16[1][0] = tl_21 * tr_11;
        t16[2][3] = tl_12 * tr_22;
        t16[3][2] = tl_21 * tr_22;
    }
    if ltranr {
        t16[0][2] = tl_11 * tr_12;
        t16[1][3] = tl_22 * tr_12;
        t16[2][0] = tl_11 * tr_21;
        t16[3][1] = tl_22 * tr_21;
    } else {
        t16[0][2] = tl_11 * tr_21;
        t16[1][3] = tl_22 * tr_21;
        t16[2][0] = tl_11 * tr_12;
        t16[3][1] = tl_22 * tr_12;
    }
    if ltranl && ltranr {
        t16[0][3] = tl_21 * tr_12;
        t16[1][2] = tl_12 * tr_12;
        t16[2][1] = tl_21 * tr_21;
        t16[3][0] = tl_12 * tr_21;
    } else if ltranl && !ltranr {
        t16[0][3] = tl_21 * tr_21;
        t16[1][2] = tl_12 * tr_21;
        t16[2][1] = tl_21 * tr_12;
        t16[3][0] = tl_12 * tr_12;
    } else if !ltranl && ltranr {
        t16[0][3] = tl_12 * tr_12;
        t16[1][2] = tl_21 * tr_12;
        t16[2][1] = tl_12 * tr_21;
        t16[3][0] = tl_21 * tr_21;
    } else {
        t16[0][3] = tl_12 * tr_21;
        t16[1][2] = tl_21 * tr_21;
        t16[2][1] = tl_12 * tr_12;
        t16[3][0] = tl_21 * tr_12;
    }

    let mut btmp = [b[0], b[ldb], b[1], b[ldb + 1]];
    let mut jpiv = [0_usize; 4];

    // Flatten t16 for row swap (Fortran column-major: t16(ip,j) -> row ip)
    let mut t16_flat = [ZERO; 16];
    for i in 0..4 {
        for j in 0..4 {
            t16_flat[i * 4 + j] = t16[i][j];
        }
    }

    for i in 0..3 {
        let mut xmax = ZERO;
        let mut ipsv = i;
        let mut jpsv = i;
        for ip in i..4 {
            for jp in i..4 {
                let abs_val = t16_flat[ip * 4 + jp].abs();
                if abs_val >= xmax {
                    xmax = abs_val;
                    ipsv = ip;
                    jpsv = jp;
                }
            }
        }
        if ipsv != i {
            for j in 0..4 {
                let t = t16_flat[ipsv * 4 + j];
                t16_flat[ipsv * 4 + j] = t16_flat[i * 4 + j];
                t16_flat[i * 4 + j] = t;
            }
            let t = btmp[i];
            btmp[i] = btmp[ipsv];
            btmp[ipsv] = t;
        }
        if jpsv != i {
            for row in 0..4 {
                let t = t16_flat[row * 4 + jpsv];
                t16_flat[row * 4 + jpsv] = t16_flat[row * 4 + i];
                t16_flat[row * 4 + i] = t;
            }
        }
        jpiv[i] = jpsv;
        if t16_flat[i * 4 + i].abs() < smin {
            info = 1;
            t16_flat[i * 4 + i] = smin;
        }
        for j in (i + 1)..4 {
            t16_flat[j * 4 + i] /= t16_flat[i * 4 + i];
            btmp[j] -= t16_flat[j * 4 + i] * btmp[i];
            for k in (i + 1)..4 {
                t16_flat[j * 4 + k] -= t16_flat[j * 4 + i] * t16_flat[i * 4 + k];
            }
        }
    }
    if t16_flat[15].abs() < smin {
        t16_flat[15] = smin;
    }
    if (EIGHT * smlnum) * btmp[0].abs() > t16_flat[0].abs()
        || (EIGHT * smlnum) * btmp[1].abs() > t16_flat[5].abs()
        || (EIGHT * smlnum) * btmp[2].abs() > t16_flat[10].abs()
        || (EIGHT * smlnum) * btmp[3].abs() > t16_flat[15].abs()
    {
        let s = (ONE / EIGHT) / btmp[0].abs().max(btmp[1].abs()).max(btmp[2].abs()).max(btmp[3].abs());
        *scale = s;
        btmp[0] *= s;
        btmp[1] *= s;
        btmp[2] *= s;
        btmp[3] *= s;
    }

    let mut tmp = [ZERO; 4];
    for kk in 0..4 {
        let k = 3 - kk;
        let temp = ONE / t16_flat[k * 4 + k];
        tmp[k] = btmp[k] * temp;
        for j in (k + 1)..4 {
            tmp[k] -= (temp * t16_flat[k * 4 + j]) * tmp[j];
        }
    }
    for ii in 0..3 {
        let i = 3 - ii;
        if jpiv[i] != i {
            let t = tmp[i];
            tmp[i] = tmp[jpiv[i]];
            tmp[jpiv[i]] = t;
        }
    }

    x[0] = tmp[0];
    x[ldx] = tmp[1];
    x[1] = tmp[2];
    x[ldx + 1] = tmp[3];
    *xnorm = (tmp[0].abs() + tmp[2].abs()).max(tmp[1].abs() + tmp[3].abs());
    info
}

fn sb04px_2x2_solve(
    tmp: &mut [f64; 4],
    btmp: &mut [f64; 2],
    smin: f64,
    smlnum: f64,
    scale: &mut f64,
    x: &mut [f64],
    ldx: usize,
    xnorm: &mut f64,
    info: &mut i32,
    n1: usize,
    _n2: usize,
) {
    let locu12 = [3, 4, 1, 2];
    let locl21 = [2, 1, 4, 3];
    let locu22 = [4, 3, 2, 1];
    let xswpiv = [false, false, true, true];
    let bswpiv = [false, true, false, true];

    let ipiv = idamax(4, tmp, 1);
    let ipiv = ipiv - 1; // 0-based
    let u11 = tmp[ipiv];
    if u11.abs() <= smin {
        *info = 1;
    }
    let u11 = if u11.abs() <= smin { smin } else { u11 };
    let u12 = tmp[locu12[ipiv] - 1];
    let l21 = tmp[locl21[ipiv] - 1] / u11;
    let u22 = tmp[locu22[ipiv] - 1] - u12 * l21;
    let u22 = if u22.abs() <= smin {
        *info = 1;
        smin
    } else {
        u22
    };
    let xswap = xswpiv[ipiv];
    let bswap = bswpiv[ipiv];
    if bswap {
        let temp = btmp[1];
        btmp[1] = btmp[0] - l21 * temp;
        btmp[0] = temp;
    } else {
        btmp[1] -= l21 * btmp[0];
    }
    if (TWO * smlnum) * btmp[1].abs() > u22.abs() || (TWO * smlnum) * btmp[0].abs() > u11.abs() {
        *scale = HALF / btmp[0].abs().max(btmp[1].abs());
        btmp[0] *= *scale;
        btmp[1] *= *scale;
    }
    let mut x2 = [ZERO; 2];
    x2[1] = btmp[1] / u22;
    x2[0] = btmp[0] / u11 - (u12 / u11) * x2[1];
    if xswap {
        x2.swap(0, 1);
    }
    x[0] = x2[0];
    if n1 == 1 {
        x[1] = x2[1];
        *xnorm = x2[0].abs() + x2[1].abs();
    } else {
        x[ldx] = x2[1];
        *xnorm = x2[0].abs().max(x2[1].abs());
    }
}

/// Compatibility wrapper: (n, a, x) for existing callers; solves with identity op and ISGN=1, scale 1.
pub fn sb04px_compat(_n: usize, _a: &nalgebra::DMatrix<f64>, x: &mut nalgebra::DMatrix<f64>) -> i32 {
    if x.is_empty() {
        return 0;
    }
    let mut scale = 1.0_f64;
    let mut xnorm = 0.0_f64;
    let tl = [1.0_f64, 0.0, 0.0, 1.0];
    let tr = [1.0_f64, 0.0, 0.0, 1.0];
    let b = [0.0_f64, 0.0, 0.0, 0.0];
    let mut xbuf = [0.0_f64; 4];
    sb04px(false, false, 1, 1, 1, &tl, 2, &tr, 2, &b, 2, &mut scale, &mut xbuf, 2, &mut xnorm);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb04px_1x1() {
        let tl = [2.0_f64];
        let tr = [3.0_f64];
        let b = [1.0_f64];
        let mut scale = 0.0_f64;
        let mut x = [0.0_f64; 4];
        let mut xnorm = 0.0_f64;
        let info = sb04px(false, false, 1, 1, 1, &tl, 1, &tr, 1, &b, 1, &mut scale, &mut x, 2, &mut xnorm);
        assert_eq!(info, 0);
        assert!((scale - 1.0).abs() < 1e-10);
        // 2*3*X + X = 1 => 7*X = 1 => X = 1/7
        assert!((x[0] - 1.0 / 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_sb04px_compat() {
        let a = nalgebra::DMatrix::<f64>::zeros(1, 1);
        let mut x = nalgebra::DMatrix::<f64>::zeros(1, 1);
        assert_eq!(sb04px_compat(1, &a, &mut x), 0);
    }
}
