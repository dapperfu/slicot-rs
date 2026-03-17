//! SB03OY — Solve 2×2 Lyapunov equation for Cholesky factor (SLICOT).
//!
//! Full implementation: continuous or discrete, op(K)=K or K', returns Cholesky factor U
//! and auxiliary matrices B, A. Uses DLANV2, SB03OV, DLAPY2, DLAPY3, DLAMCH, DLABAD.

use crate::sb03::sb03ov::{dlapy3, sb03ov};

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;

/// DLAMCH('P') = epsilon, DLAMCH('S') = safe min
fn dlamch_p() -> f64 {
    f64::EPSILON
}
fn dlamch_s() -> f64 {
    f64::MIN_POSITIVE
}

/// DLABAD: rescale SMLNUM and BIGNUM to avoid underflow/overflow
fn dlabad(smlnum: &mut f64, bignum: &mut f64) {
    let eps = dlamch_p();
    if *smlnum <= ZERO {
        return;
    }
    *smlnum = (*smlnum / eps).min(*bignum);
    *bignum = ONE / (*smlnum);
}

#[inline]
fn dlapy2(x: f64, y: f64) -> f64 {
    x.hypot(y)
}

/// DLANV2: 2×2 real Schur form and eigenvalues (LAPACK).
/// Overwrites a,b,c,d with Schur form; returns (rt1r,rt1i), (rt2r,rt2i), (cs,sn).
/// For complex conjugate pair, rt1i > 0.
fn dlanv2(
    a: &mut f64,
    b: &mut f64,
    c: &mut f64,
    d: &mut f64,
    rt1r: &mut f64,
    rt1i: &mut f64,
    rt2r: &mut f64,
    rt2i: &mut f64,
    cs: &mut f64,
    sn: &mut f64,
) {
    let mut aa = *a;
    let mut bb = *b;
    let mut cc = *c;
    let mut dd = *d;
    let tr = aa + dd;
    let det = aa * dd - bb * cc;
    let discr = tr * tr - 4.0 * det;
    if discr >= ZERO {
        let sqrt_d = discr.sqrt();
        let (l1, l2) = if tr >= ZERO {
            ((tr + sqrt_d) / 2.0, (tr - sqrt_d) / 2.0)
        } else {
            ((tr - sqrt_d) / 2.0, (tr + sqrt_d) / 2.0)
        };
        *rt1r = l1;
        *rt1i = ZERO;
        *rt2r = l2;
        *rt2i = ZERO;
        if cc != ZERO {
            *cs = (l1 - dd) / cc;
            *sn = -bb / cc;
        } else if bb != ZERO {
            *cs = -cc / bb;
            *sn = (l1 - aa) / bb;
        } else {
            *cs = ONE;
            *sn = ZERO;
        }
        *a = l1;
        *b = ZERO;
        *c = ZERO;
        *d = l2;
    } else {
        let sr = tr / 2.0;
        let si = (-discr).sqrt() / 2.0;
        *rt1r = sr;
        *rt1i = si;
        *rt2r = sr;
        *rt2i = -si;
        let r = (aa - sr).hypot(bb);
        if r > ZERO {
            *cs = (aa - sr) / r;
            *sn = bb / r;
        } else {
            *cs = ONE;
            *sn = ZERO;
        }
        *a = sr;
        *b = r;
        *c = -cc / r;
        *d = sr;
    }
}

/// SB03OY full: solve 2×2 Lyapunov for Cholesky factor U.
///
/// DISCR: false = continuous (S'*X + X*S = -ISGN*scale^2*R'*R), true = discrete.
/// LTRANS: op(K)=K' if true, op(K)=K if false.
/// ISGN: 1 or -1.
/// S (2×2), R (2×2 upper tri) input; on exit R = U (Cholesky of X), S = B, A (2×2 upper tri) output.
///
/// Returns: 0 = success, 1 = (nearly) singular warning, 2 = stability check failed, 4 = real eigenvalues.
pub fn sb03oy_full(
    discr: bool,
    ltrans: bool,
    isgn: i32,
    s: &mut [f64],
    lds: usize,
    r: &mut [f64],
    ldr: usize,
    a: &mut [f64],
    lda: usize,
    scale: &mut f64,
) -> i32 {
    let sgn = isgn as f64;
    let mut s11 = s[0];
    let mut s12 = s[1];
    let mut s21 = s[lds];
    let mut s22 = s[lds + 1];

    let mut eps = dlamch_p();
    let mut smlnum = dlamch_s();
    let mut bignum = ONE / smlnum;
    dlabad(&mut smlnum, &mut bignum);
    smlnum = smlnum * FOUR / eps;
    bignum = ONE / smlnum;

    let mut smin = smlnum;
    *scale = ONE;

    let mut tempr = 0.0_f64;
    let mut tempi = 0.0_f64;
    let mut e1 = 0.0_f64;
    let mut e2 = 0.0_f64;
    let mut csp = [0.0_f64; 2];
    let mut csq = [0.0_f64; 2];
    dlanv2(
        &mut s11,
        &mut s12,
        &mut s21,
        &mut s22,
        &mut tempr,
        &mut tempi,
        &mut e1,
        &mut e2,
        &mut csp[0],
        &mut csq[0],
    );
    csp[1] = ZERO;
    csq[1] = ZERO;

    let mut info = 0;
    if tempi == ZERO {
        return 4;
    }
    let absb = dlapy2(e1, e2);
    if discr {
        if sgn * (absb - ONE) >= ZERO {
            return 2;
        }
    } else if sgn * e1 >= ZERO {
        return 2;
    }

    let mut temp = [s[0] - e1, e2];
    if ltrans {
        temp[1] = -e2;
    }
    let mut snq = 0.0_f64;
    sb03ov(&mut temp, s[lds], smlnum, &mut csq, &mut snq);

    let t1 = csq[0] * s[1] - snq * s[0];
    let t2 = csq[1] * s[1];
    let tempr_q = csq[0] * s[lds + 1] - snq * s[lds];
    let tempi_q = csq[1] * s[lds + 1];
    let t = [
        csq[0] * t1 - csq[1] * t2 + snq * tempr_q,
        csq[0] * t2 + csq[1] * t1 + snq * tempi_q,
    ];

    let (p1, p2, p3r, p3i, snp) = if ltrans {
        let mut temp_p = [csq[0] * r[ldr + 1] - snq * r[1], -csq[1] * r[ldr + 1]];
        let mut snp = 0.0_f64;
        sb03ov(&mut temp_p, -snq * r[0], smlnum, &mut csp, &mut snp);
        let p1 = temp_p[0];
        let temp2_1 = csq[0] * r[1] + snq * r[ldr + 1];
        let temp2_2 = -csq[1] * r[1];
        let p2_1 = csp[0] * temp2_1 - csp[1] * temp2_2 + snp * (csq[0] * r[0]);
        let p2_2 = -csp[0] * temp2_2 - csp[1] * temp2_1 - snp * (-csq[1] * r[0]);
        let p3r = csp[0] * (csq[0] * r[0]) + csp[1] * (-csq[1] * r[0]) - snp * temp2_1;
        let p3i = csp[0] * (-csq[1] * r[0]) - csp[1] * (csq[0] * r[0]) - snp * temp2_2;
        (p1, [p2_1, p2_2], p3r, p3i, snp)
    } else {
        let mut temp_p = [csq[0] * r[0] + snq * r[1], csq[1] * r[0]];
        let mut snp = 0.0_f64;
        sb03ov(&mut temp_p, snq * r[ldr + 1], smlnum, &mut csp, &mut snp);
        let p1 = temp_p[0];
        let temp2_1 = csq[0] * r[1] - snq * r[0];
        let temp2_2 = csq[1] * r[1];
        let p2_1 = csp[0] * temp2_1 - csp[1] * temp2_2 + snp * (csq[0] * r[ldr + 1]);
        let p2_2 = csp[0] * temp2_2 + csp[1] * temp2_1 + snp * (csq[1] * r[ldr + 1]);
        let p3r = csp[0] * (csq[0] * r[ldr + 1]) + csp[1] * (csq[1] * r[ldr + 1]) - snp * temp2_1;
        let p3i = csp[1] * (csq[0] * r[ldr + 1]) - csp[0] * (csq[1] * r[ldr + 1]) + snp * temp2_2;
        (p1, [p2_1, p2_2], p3r, p3i, snp)
    };

    let (p3, dp1, dp2) = if p3i == ZERO {
        (p3r.abs(), p3r.signum(), ZERO)
    } else {
        let p3_abs = dlapy2(p3r, p3i);
        (p3_abs, p3r / p3_abs, -p3i / p3_abs)
    };

    let alpha = if discr {
        ((ONE - absb) * (ONE + absb)).abs().sqrt()
    } else {
        (TWO * e1).abs().sqrt()
    };

    let mut scaloc = ONE;
    let mut alpha_use = alpha;
    if alpha < smin {
        alpha_use = smin;
        info = 1;
    }
    let abst = p1.abs();
    if alpha_use < ONE && abst > ONE && abst > bignum * alpha_use {
        scaloc = ONE / abst;
    }
    let (mut p1_s, mut p2_s, mut p3_s) = (p1, p2, p3);
    if scaloc != ONE {
        p1_s *= scaloc;
        p2_s[0] *= scaloc;
        p2_s[1] *= scaloc;
        p3_s *= scaloc;
        *scale *= scaloc;
    }
    let mut v1 = p1_s / alpha_use;

    let (mut v2, mut y) = if discr {
        let g1 = (ONE - e1) * (ONE + e1) + e2 * e2;
        let g2 = -TWO * e1 * e2;
        let mut absg = dlapy2(g1, g2);
        if absg < smin {
            absg = smin;
            info = 1;
        }
        let mut temp1 = [
            sgn * alpha_use * p2_s[0] + v1 * (e1 * t[0] - e2 * t[1]),
            sgn * alpha_use * p2_s[1] + v1 * (e1 * t[1] + e2 * t[0]),
        ];
        let mut abst = temp1[0].abs().max(temp1[1].abs());
        if absg < ONE && abst > ONE && abst > bignum * absg {
            let sc = ONE / abst;
            v1 *= sc;
            temp1[0] *= sc;
            temp1[1] *= sc;
            p1_s *= sc;
            p2_s[0] *= sc;
            p2_s[1] *= sc;
            p3_s *= sc;
            *scale *= sc;
        }
        temp1[0] /= absg;
        temp1[1] /= absg;
        let mut v2_ = [
            g1 * temp1[0] + g2 * temp1[1],
            g1 * temp1[1] - g2 * temp1[0],
        ];
        abst = v2_[0].abs().max(v2_[1].abs());
        if absg < ONE && abst > ONE && abst > bignum * absg {
            let sc = ONE / abst;
            v1 *= sc;
            v2_[0] *= sc;
            v2_[1] *= sc;
            p1_s *= sc;
            p2_s[0] *= sc;
            p2_s[1] *= sc;
            p3_s *= sc;
            *scale *= sc;
        }
        v2_[0] /= absg;
        v2_[1] /= absg;
        let mut temp2 = [p1_s * t[0] - TWO * e2 * p2_s[1], p1_s * t[1] + TWO * e2 * p2_s[0]];
        abst = temp2[0].abs().max(temp2[1].abs());
        if absg < ONE && abst > ONE && abst > bignum * absg {
            let sc = ONE / abst;
            temp2[0] *= sc;
            temp2[1] *= sc;
            v1 *= sc;
            v2_[0] *= sc;
            v2_[1] *= sc;
            p3_s *= sc;
            *scale *= sc;
        }
        temp2[0] /= absg;
        temp2[1] /= absg;
        let mut y_ = [
            -(g1 * temp2[0] + g2 * temp2[1]),
            -(g1 * temp2[1] - g2 * temp2[0]),
        ];
        abst = y_[0].abs().max(y_[1].abs());
        if absg < ONE && abst > ONE && abst > bignum * absg {
            let sc = ONE / abst;
            y_[0] *= sc;
            y_[1] *= sc;
            v1 *= sc;
            v2_[0] *= sc;
            v2_[1] *= sc;
            p3_s *= sc;
            *scale *= sc;
        }
        (v2_, [y_[0] / absg, y_[1] / absg])
    } else {
        let mut absb_use = absb;
        if absb_use < smin {
            absb_use = smin;
            info = 1;
        }
        let mut temp1 = [
            sgn * alpha_use * p2_s[0] + v1 * t[0],
            sgn * alpha_use * p2_s[1] + v1 * t[1],
        ];
        let mut abst = temp1[0].abs().max(temp1[1].abs());
        if absb_use < ONE && abst > ONE && abst > bignum * absb_use {
            let sc = ONE / abst;
            v1 *= sc;
            temp1[0] *= sc;
            temp1[1] *= sc;
            p2_s[0] *= sc;
            p2_s[1] *= sc;
            p3_s *= sc;
            *scale *= sc;
        }
        temp1[0] /= TWO * absb_use;
        temp1[1] /= TWO * absb_use;
        let mut v2_ = [
            -(e1 * temp1[0] + e2 * temp1[1]),
            -(e1 * temp1[1] - e2 * temp1[0]),
        ];
        abst = v2_[0].abs().max(v2_[1].abs());
        if absb_use < ONE && abst > ONE && abst > bignum * absb_use {
            let sc = ONE / abst;
            v1 *= sc;
            v2_[0] *= sc;
            v2_[1] *= sc;
            p2_s[0] *= sc;
            p2_s[1] *= sc;
            p3_s *= sc;
            *scale *= sc;
        }
        v2_[0] /= absb_use;
        v2_[1] /= absb_use;
        let y_ = [p2_s[0] - alpha_use * v2_[0], p2_s[1] - alpha_use * v2_[1]];
        (v2_, y_)
    };

    let mut v3 = dlapy3(p3_s, y[0], y[1]);
    scaloc = ONE;
    if alpha_use < ONE && v3 > ONE && v3 > bignum * alpha_use {
        scaloc = ONE / v3;
    }
    if scaloc != ONE {
        v1 *= scaloc;
        v2[0] *= scaloc;
        v2[1] *= scaloc;
        v3 *= scaloc;
        p3_s *= scaloc;
        *scale *= scaloc;
    }
    v3 /= alpha_use;

    let mut cst = [0.0_f64; 2];
    let mut snt = 0.0_f64;
    let (r11, r12, r21, r22, dt1, dt2) = if ltrans {
        let mut x11 = [csq[0] * v3, csq[1] * v3];
        let x21_1 = snq * v3;
        let mut x12 = [
            csq[0] * v2[0] + csq[1] * v2[1] - snq * v1,
            -csq[0] * v2[1] + csq[1] * v2[0],
        ];
        let mut x22 = [csq[0] * v1 + snq * v2[0], -csq[1] * v1 - snq * v2[1]];
        x22[1] = -x22[1];
        sb03ov(&mut x22, x21_1, smlnum, &mut cst, &mut snt);
        r[ldr + 1] = x22[0];
        r[1] = cst[0] * x12[0] - cst[1] * x12[1] + snt * x11[0];
        let tempr = cst[0] * x11[0] + cst[1] * x11[1] - snt * x12[0];
        let tempi = cst[0] * x11[1] - cst[1] * x11[0] - snt * x12[1];
        let (r00, d1, d2) = if tempi == ZERO {
            (tempr.abs(), tempr.signum(), ZERO)
        } else {
            let rr = dlapy2(tempr, tempi);
            (rr, tempr / rr, -tempi / rr)
        };
        r[0] = r00;
        (r00, r[1], 0.0_f64, r[ldr + 1], d1, d2)
    } else {
        let mut x11 = [csq[0] * v1 - snq * v2[0], -csq[1] * v1 + snq * v2[1]];
        let x21_1 = -snq * v3;
        let mut x12 = [
            csq[0] * v2[0] + csq[1] * v2[1] + snq * v1,
            -csq[0] * v2[1] + csq[1] * v2[0],
        ];
        let mut x22 = [csq[0] * v3, csq[1] * v3];
        sb03ov(&mut x11, x21_1, smlnum, &mut cst, &mut snt);
        r[0] = x11[0];
        r[1] = cst[0] * x12[0] + cst[1] * x12[1] + snt * x22[0];
        let tempr = cst[0] * x22[0] - cst[1] * x22[1] - snt * x12[0];
        let tempi = cst[0] * x22[1] + cst[1] * x22[0] - snt * x12[1];
        let (r22_val, d1, d2) = if tempi == ZERO {
            (tempr.abs(), tempr.signum(), ZERO)
        } else {
            let rr = dlapy2(tempr, tempi);
            (rr, tempr / rr, -tempi / rr)
        };
        r[ldr + 1] = r22_val;
        (r[0], r[1], 0.0_f64, r22_val, d1, d2)
    };

    let (delta1, delta2, gamma1, gamma2, eta) = if y[0].abs() < smlnum && y[1].abs() <= smlnum {
        (ZERO, ZERO, ZERO, ZERO, alpha_use)
    } else {
        let d1 = y[0] / v3;
        let d2 = y[1] / v3;
        let (g1, g2) = (-alpha_use * d1, -alpha_use * d2);
        let eta = p3_s / v3;
        let (dd1, dd2) = if discr {
            (e1 * d1 - e2 * d2, e1 * d2 + e2 * d1)
        } else {
            (d1, d2)
        };
        (dd1, dd2, g1, g2, eta)
    };

    if ltrans {
        let x11_1 = cst[0] * e1 + cst[1] * e2;
        let x11_2 = -cst[0] * e2 + cst[1] * e1;
        let x12_1 = sgn * (cst[0] * gamma1 + cst[1] * gamma2) - snt * e1;
        let x12_2 = sgn * (-cst[0] * gamma2 + cst[1] * gamma1) - snt * e2;
        let x22_1 = cst[0] * e1 + cst[1] * e2 + sgn * snt * gamma1;
        let x22_2 = cst[0] * e2 - cst[1] * e1 - sgn * snt * gamma2;

        s[0] = cst[0] * x11_1 + cst[1] * x11_2 - snt * x12_1;
        let tempr_s = cst[0] * (snt * e1) + cst[1] * (-snt * e2) - snt * x22_1;
        let tempi_s = cst[0] * (-snt * e2) - cst[1] * (snt * e1) - snt * x22_2;
        s[lds] = dt1 * tempr_s - dt2 * tempi_s;
        let tempr_s12 = cst[0] * x12_1 - cst[1] * x12_2 + snt * x11_1;
        let tempi_s12 = cst[0] * x12_2 + cst[1] * x12_1 + snt * x11_2;
        s[1] = dt1 * tempr_s12 + dt2 * tempi_s12;
        s[lds + 1] = cst[0] * x22_1 - cst[1] * x22_2 + snt * (snt * e1);

        let tempr = dp1 * eta;
        let tempi = -dp2 * eta;
        a[0] = dt1 * (cst[0] * (csp[0] * tempr - csp[1] * tempi + csp[0] * delta1 - csp[1] * delta2) - cst[1] * (csp[0] * tempi + csp[1] * tempr - csp[1] * delta1 - csp[0] * delta2) - snt * (csp[0] * alpha_use)) + dt2 * (cst[0] * (csp[0] * alpha_use) + cst[1] * (-csp[1] * alpha_use));
        a[1] = dt1 * (cst[0] * (-csp[0] * tempr - csp[1] * tempi + csp[0] * delta1 - csp[1] * delta2) - cst[1] * (-csp[0] * tempi + csp[1] * tempr - csp[1] * delta1 - csp[0] * delta2) - snt * (csp[0] * alpha_use)) + dt2 * (cst[0] * (-csp[1] * alpha_use) - cst[1] * (csp[0] * alpha_use));
        a[lda] = ZERO;
        a[lda + 1] = cst[0] * (csp[0] * alpha_use) + cst[1] * (-csp[1] * alpha_use) + snt * (-csp[0] * tempr - csp[1] * tempi + csp[0] * delta1 - csp[1] * delta2);
    } else {
        s[0] = cst[0] * (cst[0] * e1 + cst[1] * e2) - cst[1] * (cst[0] * e2 - cst[1] * e1) + snt * (sgn * (cst[0] * gamma1 - cst[1] * gamma2) + snt * e1);
        s[1] = dt1 * (cst[0] * (sgn * (cst[0] * gamma1 - cst[1] * gamma2) + snt * e1) - cst[1] * (sgn * (-cst[0] * gamma2 - cst[1] * gamma1) - snt * e2) + snt * (cst[0] * e1 + cst[1] * e2)) + dt2 * (cst[0] * (sgn * (-cst[0] * gamma2 - cst[1] * gamma1) - snt * e2) + cst[1] * (sgn * (cst[0] * gamma1 - cst[1] * gamma2) + snt * e1) - snt * (-cst[0] * e2 + cst[1] * e1));
        s[lds] = dt1 * (cst[0] * (-snt * e1) - cst[1] * (-snt * e2) + snt * (cst[0] * e1 + cst[1] * e2 - sgn * snt * gamma1)) - dt2 * (cst[0] * (-snt * e2) + cst[1] * (-snt * e1) + snt * (-cst[0] * e2 + cst[1] * e1 + sgn * snt * gamma2));
        s[lds + 1] = cst[0] * (cst[0] * e1 + cst[1] * e2 - sgn * snt * gamma1) + cst[1] * (-cst[0] * e2 + cst[1] * e1 + sgn * snt * gamma2) - snt * (-snt * e1);

        let tempr = dp1 * eta;
        let tempi = -dp2 * eta;
        let x11_1 = csp[0] * alpha_use;
        let x11_2 = csp[1] * alpha_use;
        let x21_1 = snp * alpha_use;
        let x12_1 = csp[0] * delta1 + csp[1] * delta2 - snp * tempr;
        let x12_2 = -csp[0] * delta2 + csp[1] * delta1 - snp * tempi;
        let x22_1 = csp[0] * tempr + csp[1] * tempi + snp * delta1;
        let x22_2 = csp[0] * tempi - csp[1] * tempr - snp * delta2;
        a[0] = cst[0] * x11_1 - cst[1] * x11_2 + snt * x12_1;
        a[lda] = ZERO;
        a[1] = cst[0] * x12_1 + cst[1] * x12_2 - snt * x11_1;
        let tempr_a = cst[0] * x22_1 + cst[1] * x22_2 - snt * x21_1;
        let tempi_a = cst[0] * x22_2 - cst[1] * x22_1;
        a[lda + 1] = dt1 * tempr_a + dt2 * tempi_a;
    }

    if *scale != ONE {
        a[0] *= *scale;
        a[1] *= *scale;
        a[lda + 1] *= *scale;
    }

    info
}

/// Compatibility wrapper: (n, a, x). Calls full 2×2 solver when n=2; otherwise returns 1.
pub fn sb03oy(n: usize, a: &nalgebra::DMatrix<f64>, x: &mut nalgebra::DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    if n == 2 {
        let mut s = [a[(0, 0)], a[(0, 1)], a[(1, 0)], a[(1, 1)]];
        let mut r = [1.0_f64, 0.0, 0.0, 1.0];
        let mut a_out = [0.0_f64; 4];
        let mut scale = 1.0_f64;
        let info = sb03oy_full(false, false, 1, &mut s, 2, &mut r, 2, &mut a_out, 2, &mut scale);
        if info == 0 {
            let ldr = 2;
            x[(0, 0)] = r[0] * r[0];
            x[(0, 1)] = r[0] * r[1];
            x[(1, 0)] = r[1] * r[0];
            x[(1, 1)] = r[1] * r[1] + r[ldr + 1] * r[ldr + 1];
        }
        return info;
    }
    1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03oy_cont_2x2() {
        let mut s = [-1.0, 0.5, -0.5, -1.0];
        let mut r = [1.0, 0.0, 0.0, 1.0];
        let mut a = [0.0; 4];
        let mut scale = 1.0;
        let info = sb03oy_full(false, false, 1, &mut s, 2, &mut r, 2, &mut a, 2, &mut scale);
        assert_eq!(info, 0);
    }
}
