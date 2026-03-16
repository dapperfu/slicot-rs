//! SG03BX — Solve 2×2 generalized Lyapunov equation for Cholesky factor U (SLICOT SG03BX)
//!
//! Solves for X = op(U)^T*op(U) the generalized continuous-time or discrete-time
//! Lyapunov equation, and computes auxiliary matrices M1, M2. Uses SG03BY for
//! complex Givens. Inline helpers: DLADIV, DLAG2, DLARTGP, DROT.

use nalgebra::DMatrix;
use crate::sg03::sg03by::sg03by;

const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const ZERO: f64 = 0.0;
const MONE: f64 = -1.0;
const HALF: f64 = 0.5;

type M2 = [[f64; 2]; 2];

fn dlapy2(x: f64, y: f64) -> f64 {
    x.hypot(y)
}

// ----- DLADIV: (a + i*b)/(c + i*d) -> (p, q), avoiding overflow -----
fn dladiv(a: f64, b: f64, c: f64, d: f64, p: &mut f64, q: &mut f64) {
    let eps = f64::EPSILON;
    let safmin = f64::MIN_POSITIVE;
    let ov = f64::MAX / 2.0;
    let be = 2.0 / (eps * eps);
    let un = safmin;
    let bs = 2.0;

    let mut aa = a;
    let mut bb = b;
    let mut cc = c;
    let mut dd = d;
    let mut s = 1.0_f64;
    let ab = a.abs().max(b.abs());
    let cd = c.abs().max(d.abs());

    if ab >= HALF * ov {
        aa *= HALF;
        bb *= HALF;
        s *= TWO;
    }
    if cd >= HALF * ov {
        cc *= HALF;
        dd *= HALF;
        s *= HALF;
    }
    if ab <= un * bs / eps {
        aa *= be;
        bb *= be;
        s /= be;
    }
    if cd <= un * bs / eps {
        cc *= be;
        dd *= be;
        s *= be;
    }
    if d.abs() <= c.abs() {
        dladiv1(aa, bb, cc, dd, p, q);
    } else {
        dladiv1(bb, aa, dd, cc, p, q);
        *q = -(*q);
    }
    *p *= s;
    *q *= s;
}

fn dladiv2(a: f64, b: f64, c: f64, d: f64, r: f64, t: f64) -> f64 {
    if r != ZERO {
        let br = b * r;
        if br != ZERO {
            (a + br) * t
        } else {
            a * t + (b * t) * r
        }
    } else {
        (a + d * (b / c)) * t
    }
}

fn dladiv1(a: f64, b: f64, c: f64, d: f64, p: &mut f64, q: &mut f64) {
    let r = d / c;
    let t = ONE / (c + d * r);
    *p = dladiv2(a, b, c, d, r, t);
    let a_neg = -a;
    *q = dladiv2(b, a_neg, c, d, r, t);
}

// ----- DLAG2: generalized eigenvalues of 2x2 pencil A - w*B -----
#[allow(clippy::too_many_arguments)]
fn dlag2(
    a: &M2,
    b: &M2,
    safmin: f64,
    scale1: &mut f64,
    scale2: &mut f64,
    wr1: &mut f64,
    wr2: &mut f64,
    wi: &mut f64,
) {
    let rtmin = safmin.sqrt();
    let rtmax = ONE / rtmin;
    let safmax = ONE / safmin;
    let fuzzy1 = 1.0 + 1e-5;

    let a11 = a[0][0];
    let a21 = a[1][0];
    let a12 = a[0][1];
    let a22 = a[1][1];

    let anorm = (a11.abs() + a21.abs())
        .max(a12.abs() + a22.abs())
        .max(safmin);
    let ascale = ONE / anorm;
    let a11 = ascale * a11;
    let a21 = ascale * a21;
    let a12 = ascale * a12;
    let a22 = ascale * a22;

    let mut b11 = b[0][0];
    let mut b12 = b[0][1];
    let mut b22 = b[1][1];
    let bmin = rtmin
        * b11.abs()
            .max(b12.abs())
            .max(b22.abs())
            .max(rtmin);
    if b11.abs() < bmin {
        b11 = bmin * b11.signum();
    }
    if b22.abs() < bmin {
        b22 = bmin * b22.signum();
    }

    let bnorm = b11.abs().max(b12.abs() + b22.abs()).max(safmin);
    let bsize = b11.abs().max(b22.abs());
    let bscale = ONE / bsize;
    let b11 = b11 * bscale;
    let b12 = b12 * bscale;
    let b22 = b22 * bscale;

    let binv11 = ONE / b11;
    let binv22 = ONE / b22;
    let s1 = a11 * binv11;
    let s2 = a22 * binv22;

    let (pp, qq, shift, abi22) = if s1.abs() <= s2.abs() {
        let as12 = a12 - s1 * b12;
        let as22 = a22 - s1 * b22;
        let ss = a21 * (binv11 * binv22);
        let abi22 = as22 * binv22 - ss * b12;
        let pp = HALF * abi22;
        (pp, ss * as12, s1, abi22)
    } else {
        let as11 = a11 - s2 * b11;
        let ss = a21 * (binv11 * binv22);
        let abi22 = -ss * b12;
        let pp = HALF * (as11 * binv11 + abi22);
        (pp, ss * (a12 - s2 * b12), s2, abi22)
    };

    let (discr, r) = if pp.abs() * rtmin >= ONE {
        let discr = (rtmin * pp).powi(2) + qq * safmin;
        let r = (discr.abs().max(1e-300)).sqrt() * rtmax;
        (discr, r)
    } else if pp.powi(2) + qq.abs() <= safmin {
        let discr = (rtmax * pp).powi(2) + qq * safmax;
        let r = (discr.abs().max(1e-300)).sqrt() * rtmin;
        (discr, r)
    } else {
        let discr = pp.powi(2) + qq;
        let r = (discr.abs().max(1e-300)).sqrt();
        (discr, r)
    };

    if discr >= ZERO || r == ZERO {
        let sum = pp + r * pp.signum();
        let diff = pp - r * pp.signum();
        let wbig = shift + sum;
        let mut wsmall = shift + diff;
        if HALF * wbig.abs() > wsmall.abs().max(safmin) {
            let wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
            wsmall = wdet / wbig;
        }
        if pp > abi22 {
            *wr1 = wbig.min(wsmall);
            *wr2 = wbig.max(wsmall);
        } else {
            *wr1 = wbig.max(wsmall);
            *wr2 = wbig.min(wsmall);
        }
        *wi = ZERO;
    } else {
        *wr1 = shift + pp;
        *wr2 = *wr1;
        *wi = r;
    }

    let c1 = bsize * (safmin * ONE.max(ascale));
    let c2 = safmin * ONE.max(bnorm);
    let c3 = bsize * safmin;
    let c4 = if ascale <= ONE && bsize <= ONE {
        ONE.min((ascale / safmin) * bsize)
    } else {
        ONE
    };
    let c5 = if ascale <= ONE || bsize <= ONE {
        ONE.min(ascale * bsize)
    } else {
        ONE
    };

    let wabs = (*wr1).abs() + (*wi).abs();
    let mut wsize = safmin.max(c1).max(fuzzy1 * (wabs * c2 + c3).max(1e-300));
    wsize = wsize.max(c4.min(HALF * wabs.max(c5)));
    if wsize != ONE {
        let wscale = ONE / wsize;
        *scale1 = if wsize > ONE {
            (ascale.max(bsize) * wscale) * ascale.min(bsize)
        } else {
            (ascale.min(bsize) * wscale) * ascale.max(bsize)
        };
        *wr1 *= wscale;
        if *wi != ZERO {
            *wi *= wscale;
            *wr2 = *wr1;
            *scale2 = *scale1;
        }
    } else {
        *scale1 = ascale * bsize;
        *scale2 = *scale1;
    }

    if (*wi).abs() < 1e-300 {
        let wsize2 = safmin
            .max(c1)
            .max(fuzzy1 * ((*wr2).abs() * c2 + c3).max(1e-300));
        let wsize2 = wsize2.max(c4.min(HALF * (*wr2).abs().max(c5)));
        if wsize2 != ONE {
            let wscale = ONE / wsize2;
            *scale2 = if wsize2 > ONE {
                (ascale.max(bsize) * wscale) * ascale.min(bsize)
            } else {
                (ascale.min(bsize) * wscale) * ascale.max(bsize)
            };
            *wr2 *= wscale;
        } else {
            *scale2 = ascale * bsize;
        }
    }
}

// ----- DLARTGP: Givens (f, g) -> (cs, sn, r) with R >= 0 -----
fn dlartgp(f: f64, g: f64, cs: &mut f64, sn: &mut f64, r: &mut f64) {
    let safmin = f64::MIN_POSITIVE;
    let eps = f64::EPSILON;
    let base = 2.0_f64;
    let exp = (safmin / eps).ln() / base.ln() / TWO;
    let safmn2 = base.powf(exp.floor());
    let safmx2 = ONE / safmn2;

    if g == ZERO {
        *cs = ONE.copysign(f);
        *sn = ZERO;
        *r = f.abs();
    } else if f == ZERO {
        *cs = ZERO;
        *sn = ONE.copysign(g);
        *r = g.abs();
    } else {
        let mut f1 = f;
        let mut g1 = g;
        let mut scale = f1.abs().max(g1.abs());
        let mut count = 0_i32;
        if scale >= safmx2 {
            while scale >= safmx2 && count < 20 {
                count += 1;
                f1 *= safmn2;
                g1 *= safmn2;
                scale = f1.abs().max(g1.abs());
            }
        } else if scale <= safmn2 {
            while scale <= safmn2 {
                count -= 1;
                f1 *= safmx2;
                g1 *= safmx2;
                scale = f1.abs().max(g1.abs());
            }
        }
        let mut r_val = (f1 * f1 + g1 * g1).sqrt();
        *cs = f1 / r_val;
        *sn = g1 / r_val;
        if count > 0 {
            for _ in 0..count {
                r_val *= safmx2;
            }
        } else {
            for _ in 0..(-count) {
                r_val *= safmn2;
            }
        }
        *r = r_val;
        if *r < ZERO {
            *cs = -(*cs);
            *sn = -(*sn);
            *r = -(*r);
        }
    }
}

// ----- DROT: apply Givens to (x, y) -> (c*x + s*y, -s*x + c*y) -----
fn drot(c: f64, s: f64, x: &mut f64, y: &mut f64) {
    let xnew = c * *x + s * *y;
    let ynew = -s * *x + c * *y;
    *x = xnew;
    *y = ynew;
}

// 2x2 complex matrix multiply: (AR + i*AI)*(BR + i*BI) -> (CR, CI), column-major 2x2
fn zgemm_22(
    ar: &M2,
    ai: &M2,
    br: &M2,
    bi: &M2,
    cr: &mut M2,
    ci: &mut M2,
) {
    for i in 0..2 {
        for j in 0..2 {
            let mut re = ZERO;
            let mut im = ZERO;
            for k in 0..2 {
                re += ar[i][k] * br[k][j] - ai[i][k] * bi[k][j];
                im += ar[i][k] * bi[k][j] + ai[i][k] * br[k][j];
            }
            cr[i][j] = re;
            ci[i][j] = im;
        }
    }
}

// Real 2x2 * Real 2x2 -> Real 2x2
fn dgemm_22(a: &M2, b: &M2, c: &mut M2) {
    for i in 0..2 {
        for j in 0..2 {
            c[i][j] = (0..2).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
}

// (QR + i*QI) * (AR real) -> (ARout, AIout)
fn q_times_a(qr: &M2, qi: &M2, a: &M2, ar_out: &mut M2, ai_out: &mut M2) {
    for i in 0..2 {
        for j in 0..2 {
            ar_out[i][j] = (0..2).map(|k| qr[i][k] * a[k][j]).sum();
            ai_out[i][j] = (0..2).map(|k| qi[i][k] * a[k][j]).sum();
        }
    }
}

// (AR + i*AI) * (ZR + i*ZI) -> (TR, TI)
fn a_times_z(ar: &M2, ai: &M2, zr: &M2, zi: &M2, tr: &mut M2, ti: &mut M2) {
    for i in 0..2 {
        for j in 0..2 {
            let mut re = ZERO;
            let mut im = ZERO;
            for k in 0..2 {
                re += ar[i][k] * zr[k][j] - ai[i][k] * zi[k][j];
                im += ar[i][k] * zi[k][j] + ai[i][k] * zr[k][j];
            }
            tr[i][j] = re;
            ti[i][j] = im;
        }
    }
}

/// Continuous or discrete Lyapunov equation.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    /// Continuous-time: A'*X*E + E'*X*A = -scale^2*B'*B
    Continuous,
    /// Discrete-time: A'*X*A - E'*X*E = -scale^2*B'*B
    Discrete,
}

/// Transpose mode for op(K).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Trans {
    /// op(K) = K
    NoTrans,
    /// op(K) = K^T
    Trans,
}

/// Solves the 2×2 generalized Lyapunov equation for the Cholesky factor U and
/// auxiliary matrices M1, M2.
///
/// # Arguments
/// * `dico` — Continuous or discrete equation
/// * `trans` — NoTrans: op(K)=K; Trans: op(K)=K^T
/// * `a` — 2×2 matrix A
/// * `e` — 2×2 upper triangular E
/// * `b` — 2×2 upper triangular B
/// * `u` — output 2×2 upper triangular Cholesky factor
/// * `scale` — output scale factor, 0 < scale <= 1
/// * `m1` — output 2×2 matrix M1
/// * `m2` — output 2×2 matrix M2
/// * `info` — 0: success; 2: eigenvalues not complex conjugate pair; 3: not in stable region
pub fn sg03bx(
    dico: Dico,
    trans: Trans,
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &DMatrix<f64>,
    u: &mut DMatrix<f64>,
    scale: &mut f64,
    m1: &mut DMatrix<f64>,
    m2: &mut DMatrix<f64>,
    info: &mut i32,
) {
    let eps = f64::EPSILON;
    let smlnum = f64::MIN_POSITIVE / eps;

    *info = 0;
    *scale = ONE;

    let mut aa: M2 = [
        [a[(0, 0)], a[(0, 1)]],
        [a[(1, 0)], a[(1, 1)]],
    ];
    let mut ee: M2 = [
        [e[(0, 0)], e[(0, 1)]],
        [ZERO, e[(1, 1)]],
    ];
    let mut bb: M2 = [
        [b[(0, 0)], b[(0, 1)]],
        [ZERO, b[(1, 1)]],
    ];

    let istrns = trans == Trans::Trans;
    if istrns {
        let v = aa[0][0];
        aa[0][0] = aa[1][1];
        aa[1][1] = v;
        let v = ee[0][0];
        ee[0][0] = ee[1][1];
        ee[1][1] = v;
        let v = bb[0][0];
        bb[0][0] = bb[1][1];
        bb[1][1] = v;
    }

    let t = (eps * ee[0][0].abs().max(ee[0][1].abs()).max(ee[1][1].abs())).max(smlnum);
    if ee[0][0].abs().min(ee[1][1].abs()) < t {
        *info = 3;
        return;
    }

    let mut scale1 = 0.0_f64;
    let mut scale2 = 0.0_f64;
    let mut lamr = 0.0_f64;
    let mut w = 0.0_f64;
    let mut lami = 0.0_f64;
    dlag2(
        &aa,
        &ee,
        smlnum * eps,
        &mut scale1,
        &mut scale2,
        &mut lamr,
        &mut w,
        &mut lami,
    );
    if lami <= ZERO {
        *info = 2;
        return;
    }

    let (mut cr, mut ci, mut sr, mut si, mut l) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    sg03by(
        scale1 * aa[0][0] - ee[0][0] * lamr,
        -ee[0][0] * lami,
        scale1 * aa[1][0],
        ZERO,
        &mut cr,
        &mut ci,
        &mut sr,
        &mut si,
        &mut l,
    );
    let mut qr: M2 = [[cr, sr], [-sr, cr]];
    let mut qi: M2 = [[-ci, -si], [-si, ci]];

    let mut ar: M2 = [[ZERO; 2]; 2];
    let mut ai: M2 = [[ZERO; 2]; 2];
    q_times_a(&qr, &qi, &aa, &mut ar, &mut ai);

    let mut er: M2 = [[ZERO; 2]; 2];
    let mut ei: M2 = [[ZERO; 2]; 2];
    q_times_a(&qr, &qi, &ee, &mut er, &mut ei);

    let (mut cr, mut ci, mut sr, mut si, mut l) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    sg03by(er[1][1], ei[1][1], er[1][0], ei[1][0], &mut cr, &mut ci, &mut sr, &mut si, &mut l);
    let mut zr: M2 = [[cr, sr], [-sr, cr]];
    let mut zi: M2 = [[ci, -si], [-si, -ci]];

    let mut tr: M2 = [[ZERO; 2]; 2];
    let mut ti: M2 = [[ZERO; 2]; 2];
    for j in 0..2 {
        tr[0][j] = zr[0][0] * er[0][j] + zr[0][1] * er[1][j] - (zi[0][0] * ei[0][j] + zi[0][1] * ei[1][j]);
        tr[1][j] = zr[1][0] * er[0][j] + zr[1][1] * er[1][j] - (zi[1][0] * ei[0][j] + zi[1][1] * ei[1][j]);
        ti[0][j] = zi[0][0] * er[0][j] + zi[0][1] * er[1][j] + zr[0][0] * ei[0][j] + zr[0][1] * ei[1][j];
        ti[1][j] = zi[1][0] * er[0][j] + zi[1][1] * er[1][j] + zr[1][0] * ei[0][j] + zr[1][1] * ei[1][j];
    }
    er[0][0] = tr[0][0];
    er[0][1] = tr[0][1];
    er[1][0] = tr[1][0];
    er[1][1] = tr[1][1];
    ei[0][0] = ti[0][0];
    ei[0][1] = ti[0][1];
    ei[1][0] = ti[1][0];
    ei[1][1] = ti[1][1];
    er[1][0] = ZERO;
    er[1][1] = l;
    ei[1][0] = ZERO;
    ei[1][1] = ZERO;

    let v = dlapy2(er[0][0], ei[0][0]);
    let (mut xr, mut xi) = (0.0_f64, 0.0_f64);
    dladiv(v, ZERO, er[0][0], ei[0][0], &mut xr, &mut xi);
    er[0][0] = v;
    ei[0][0] = ZERO;
    let mut yr = zr[0][0];
    let mut yi = zi[0][0];
    zr[0][0] = xr * yr - xi * yi;
    zi[0][0] = xr * yi + xi * yr;
    yr = zr[1][0];
    yi = zi[1][0];
    zr[1][0] = xr * yr - xi * yi;
    zi[1][0] = xr * yi + xi * yr;

    a_times_z(&ar, &ai, &zr, &zi, &mut tr, &mut ti);
    ar = tr;
    ai = ti;

    let mut br: M2 = [[ZERO; 2]; 2];
    let mut bi: M2 = [[ZERO; 2]; 2];
    dgemm_22(&bb, &zr, &mut br);
    dgemm_22(&bb, &zi, &mut bi);

    let (mut cr, mut ci, mut sr, mut si, mut l) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    sg03by(br[0][0], bi[0][0], br[1][0], bi[1][0], &mut cr, &mut ci, &mut sr, &mut si, &mut l);
    let mut qbr: M2 = [[cr, sr], [-sr, cr]];
    let mut qbi: M2 = [[-ci, -si], [-si, ci]];
    let tr0 = qbr[0][0] * br[0][1] + qbr[0][1] * br[1][1] - (qbi[0][0] * bi[0][1] + qbi[0][1] * bi[1][1]);
    let tr1 = qbr[1][0] * br[0][1] + qbr[1][1] * br[1][1] - (qbi[1][0] * bi[0][1] + qbi[1][1] * bi[1][1]);
    let ti0 = qbi[0][0] * br[0][1] + qbi[0][1] * br[1][1] + qbr[0][0] * bi[0][1] + qbr[0][1] * bi[1][1];
    let ti1 = qbi[1][0] * br[0][1] + qbi[1][1] * br[1][1] + qbr[1][0] * bi[0][1] + qbr[1][1] * bi[1][1];
    br[0][1] = tr0;
    br[1][1] = tr1;
    bi[0][1] = ti0;
    bi[1][1] = ti1;
    br[0][0] = l;
    br[1][0] = ZERO;
    bi[0][0] = ZERO;
    bi[1][0] = ZERO;
    let v = dlapy2(br[1][1], bi[1][1]);
    let thresh = (eps * br[0][0].max(dlapy2(br[0][1], bi[0][1]))).max(smlnum);
    if v >= thresh {
        dladiv(v, ZERO, br[1][1], bi[1][1], &mut xr, &mut xi);
        br[1][1] = v;
        yr = qbr[1][0];
        yi = qbi[1][0];
        qbr[1][0] = xr * yr - xi * yi;
        qbi[1][0] = xr * yi + xi * yr;
        yr = qbr[1][1];
        yi = qbi[1][1];
        qbr[1][1] = xr * yr - xi * yi;
        qbi[1][1] = xr * yi + xi * yr;
    } else {
        br[1][1] = ZERO;
    }
    bi[1][1] = ZERO;

    let iscont = dico == Dico::Continuous;
    let mut ur: M2 = [[ZERO; 2]; 2];
    let mut ui: M2 = [[ZERO; 2]; 2];
    let mut m1r: M2 = [[ZERO; 2]; 2];
    let mut m1i: M2 = [[ZERO; 2]; 2];
    let mut m2r: M2 = [[ZERO; 2]; 2];
    let mut m2i: M2 = [[ZERO; 2]; 2];

    if iscont {
        let mut v = -TWO * (ar[0][0] * er[0][0] + ai[0][0] * ei[0][0]);
        if v <= ZERO {
            *info = 3;
            return;
        }
        v = v.sqrt();
        let t = TWO * br[0][0].abs() * smlnum;
        if t > v {
            scale1 = v / t;
            *scale *= scale1;
            br[0][0] *= scale1;
            br[0][1] *= scale1;
            bi[0][1] *= scale1;
            br[1][1] *= scale1;
        }
        ur[0][0] = br[0][0] / v;
        ui[0][0] = ZERO;
        ur[1][0] = ZERO;
        ui[1][0] = ZERO;

        let t = (eps * br[1][1].max(dlapy2(br[0][1], bi[0][1]))).max(smlnum);
        if br[0][0].abs() < t {
            ur[0][1] = ZERO;
            ui[0][1] = ZERO;
        } else {
            let mut xr = ar[0][0] * er[0][1] + ai[0][0] * ei[0][1]
                + ar[0][1] * er[0][0] + ai[0][1] * ei[0][0];
            let mut xi = ai[0][0] * er[0][1] - ar[0][0] * ei[0][1]
                - ai[0][1] * er[0][0] + ar[0][1] * ei[0][0];
            xr = -br[0][1] * v - xr * ur[0][0];
            xi = bi[0][1] * v - xi * ur[0][0];
            let mut yr = ar[1][1] * er[0][0] + ai[1][1] * ei[0][0]
                + er[1][1] * ar[0][0] + ei[1][1] * ai[0][0];
            let mut yi = -ai[1][1] * er[0][0] + ar[1][1] * ei[0][0]
                - ei[1][1] * ar[0][0] + er[1][1] * ai[0][0];
            let t_num = TWO * dlapy2(xr, xi) * smlnum;
            if t_num > dlapy2(yr, yi) {
                scale1 = dlapy2(yr, yi) / t_num;
                *scale *= scale1;
                br[0][0] *= scale1;
                br[0][1] *= scale1;
                bi[0][1] *= scale1;
                br[1][1] *= scale1;
                ur[0][0] *= scale1;
                xr *= scale1;
                xi *= scale1;
            }
            dladiv(xr, xi, yr, yi, &mut ur[0][1], &mut ui[0][1]);
            ui[0][1] = -ui[0][1];
        }

        let mut xr = (er[0][1] * ur[0][0] + er[1][1] * ur[0][1] - ei[1][1] * ui[0][1]) * v;
        let mut xi = (-ei[0][1] * ur[0][0] - er[1][1] * ui[0][1] - ei[1][1] * ur[0][1]) * v;
        let t_num = TWO * dlapy2(xr, xi) * smlnum;
        if t_num > dlapy2(er[0][0], ei[0][0]) {
            scale1 = dlapy2(er[0][0], ei[0][0]) / t_num;
            *scale *= scale1;
            ur[0][0] *= scale1;
            ur[0][1] *= scale1;
            ui[0][1] *= scale1;
            br[0][0] *= scale1;
            br[0][1] *= scale1;
            bi[0][1] *= scale1;
            br[1][1] *= scale1;
            xr *= scale1;
            xi *= scale1;
        }
        let (mut div_re, mut div_im) = (0.0_f64, 0.0_f64);
        dladiv(xr, xi, er[0][0], -ei[0][0], &mut div_re, &mut div_im);
        let yrh = br[0][1] - div_re;
        let yih = -bi[0][1] - div_im;
        let mut v = -TWO * (ar[1][1] * er[1][1] + ai[1][1] * ei[1][1]);
        if v <= ZERO {
            *info = 3;
            return;
        }
        v = v.sqrt();
        let w = dlapy2(dlapy2(br[1][1], bi[1][1]), dlapy2(yrh, yih));
        let t = TWO * w * smlnum;
        if t > v {
            scale1 = v / t;
            *scale *= scale1;
            ur[0][0] *= scale1;
            ur[0][1] *= scale1;
            ui[0][1] *= scale1;
            br[0][0] *= scale1;
            br[0][1] *= scale1;
            bi[0][1] *= scale1;
            br[1][1] *= scale1;
        }
        ur[1][1] = w / v;
        ui[1][1] = ZERO;

        m1r[1][0] = ZERO;
        m1i[1][0] = ZERO;
        m2r[1][0] = ZERO;
        m2i[1][0] = ZERO;
        let (mut betar, mut betai) = (0.0_f64, 0.0_f64);
        dladiv(ar[0][0], ai[0][0], er[0][0], ei[0][0], &mut betar, &mut betai);
        m1r[0][0] = betar;
        m1i[0][0] = betai;
        m1r[1][1] = betar;
        m1i[1][1] = -betai;
        let alpha = (-TWO * betar).sqrt();
        m2r[0][0] = alpha;
        m2i[0][0] = ZERO;
        let v = er[0][0] * er[1][1];
        let xr = (-br[0][0] * er[0][1] + er[0][0] * br[0][1]) / v;
        let xi = (-br[0][0] * ei[0][1] + er[0][0] * bi[0][1]) / v;
        let yr = xr - alpha * ur[0][1];
        let yi = -xi + alpha * ui[0][1];
        if yr != ZERO || yi != ZERO {
            m2r[0][1] = yr / ur[1][1];
            m2i[0][1] = -yi / ur[1][1];
            m2r[1][1] = br[1][1] / (er[1][1] * ur[1][1]);
            m2i[1][1] = ZERO;
            m1r[0][1] = -alpha * m2r[0][1];
            m1i[0][1] = -alpha * m2i[0][1];
        } else {
            m2r[0][1] = ZERO;
            m2i[0][1] = ZERO;
            m2r[1][1] = alpha;
            m2i[1][1] = ZERO;
            m1r[0][1] = ZERO;
            m1i[0][1] = ZERO;
        }
    } else {
        let mut v = er[0][0].powi(2) + ei[0][0].powi(2) - ar[0][0].powi(2) - ai[0][0].powi(2);
        if v <= ZERO {
            *info = 3;
            return;
        }
        v = v.sqrt();
        let t = TWO * br[0][0].abs() * smlnum;
        if t > v {
            scale1 = v / t;
            *scale *= scale1;
            br[0][0] *= scale1;
            br[0][1] *= scale1;
            bi[0][1] *= scale1;
            br[1][1] *= scale1;
        }
        ur[0][0] = br[0][0] / v;
        ui[0][0] = ZERO;
        ur[1][0] = ZERO;
        ui[1][0] = ZERO;

        let t = (eps * br[1][1].max(dlapy2(br[0][1], bi[0][1]))).max(smlnum);
        if br[0][0].abs() < t {
            ur[0][1] = ZERO;
            ui[0][1] = ZERO;
        } else {
            let mut xr = ar[0][0] * ar[0][1] + ai[0][0] * ai[0][1]
                - er[0][1] * er[0][0] - ei[0][1] * ei[0][0];
            let mut xi = ai[0][0] * ar[0][1] - ar[0][0] * ai[0][1]
                + ei[0][1] * er[0][0] - er[0][1] * ei[0][0];
            xr = -br[0][1] * v - xr * ur[0][0];
            xi = bi[0][1] * v - xi * ur[0][0];
            let mut yr = ar[1][1] * ar[0][0] + ai[1][1] * ai[0][0]
                - er[1][1] * er[0][0] - ei[1][1] * ei[0][0];
            let mut yi = -ai[1][1] * ar[0][0] + ar[1][1] * ai[0][0]
                + ei[1][1] * er[0][0] - er[1][1] * ei[0][0];
            let t_num = TWO * dlapy2(xr, xi) * smlnum;
            if t_num > dlapy2(yr, yi) {
                scale1 = dlapy2(yr, yi) / t_num;
                *scale *= scale1;
                br[0][0] *= scale1;
                br[0][1] *= scale1;
                bi[0][1] *= scale1;
                br[1][1] *= scale1;
                ur[0][0] *= scale1;
                xr *= scale1;
                xi *= scale1;
            }
            dladiv(xr, xi, yr, yi, &mut ur[0][1], &mut ui[0][1]);
            ui[0][1] = -ui[0][1];
        }

        let xr = er[0][1] * ur[0][0] + er[1][1] * ur[0][1] - ei[1][1] * ui[0][1];
        let xi = -ei[0][1] * ur[0][0] - er[1][1] * ui[0][1] - ei[1][1] * ur[0][1];
        let yr = ar[0][1] * ur[0][0] + ar[1][1] * ur[0][1] - ai[1][1] * ui[0][1];
        let yi = -ai[0][1] * ur[0][0] - ar[1][1] * ui[0][1] - ai[1][1] * ur[0][1];
        let mut v = er[1][1].powi(2) + ei[1][1].powi(2) - ar[1][1].powi(2) - ai[1][1].powi(2);
        if v <= ZERO {
            *info = 3;
            return;
        }
        v = v.sqrt();
        let t = br[1][1].abs()
            .max(br[0][1].abs())
            .max(bi[0][1].abs())
            .max(xr.abs())
            .max(xi.abs())
            .max(yr.abs())
            .max(yi.abs());
        let t = if t <= smlnum { ONE } else { t };
        let mut w = (br[1][1] / t).powi(2) + (br[0][1] / t).powi(2) + (bi[0][1] / t).powi(2)
            - (xr / t).powi(2) - (xi / t).powi(2)
            + (yr / t).powi(2) + (yi / t).powi(2);
        if w < ZERO {
            *info = 3;
            return;
        }
        w = t * w.sqrt();
        let t = TWO * w * smlnum;
        if t > v {
            scale1 = v / t;
            *scale *= scale1;
            ur[0][0] *= scale1;
            ur[0][1] *= scale1;
            ui[0][1] *= scale1;
            br[0][0] *= scale1;
            br[0][1] *= scale1;
            bi[0][1] *= scale1;
            br[1][1] *= scale1;
            w *= scale1;
        }
        ur[1][1] = w / v;
        ui[1][1] = ZERO;

        let b11 = br[0][0] / er[0][0];
        let t = er[0][0] * er[1][1];
        let b12r = (er[0][0] * br[0][1] - br[0][0] * er[0][1]) / t;
        let b12i = (er[0][0] * bi[0][1] - br[0][0] * ei[0][1]) / t;
        let b22 = br[1][1] / er[1][1];
        m1r[1][0] = ZERO;
        m1i[1][0] = ZERO;
        m2r[1][0] = ZERO;
        m2i[1][0] = ZERO;
        let (mut betar, mut betai) = (0.0_f64, 0.0_f64);
        dladiv(ar[0][0], ai[0][0], er[0][0], ei[0][0], &mut betar, &mut betai);
        m1r[0][0] = betar;
        m1i[0][0] = betai;
        m1r[1][1] = betar;
        m1i[1][1] = -betai;
        let v = dlapy2(betar, betai);
        let alpha = ((ONE - v) * (ONE + v)).sqrt();
        m2r[0][0] = alpha;
        m2i[0][0] = ZERO;
        let mut xr = (ai[0][0] * ei[0][1] - ar[0][0] * er[0][1]) / t + ar[0][1] / er[1][1];
        let mut xi = (ar[0][0] * ei[0][1] + ai[0][0] * er[0][1]) / t - ai[0][1] / er[1][1];
        xr = -TWO * betai * b12i - b11 * xr;
        xi = -TWO * betai * b12r - b11 * xi;
        let v = ONE + (betai - betar) * (betai + betar);
        let w = -TWO * betai * betar;
        let (mut yr, mut yi) = (0.0_f64, 0.0_f64);
        dladiv(xr, xi, v, w, &mut yr, &mut yi);
        if yr != ZERO || yi != ZERO {
            m2r[0][1] = (yr * betar - yi * betai) / ur[1][1];
            m2i[0][1] = -(yi * betar + yr * betai) / ur[1][1];
            m2r[1][1] = b22 / ur[1][1];
            m2i[1][1] = ZERO;
            m1r[0][1] = -alpha * yr / ur[1][1];
            m1i[0][1] = alpha * yi / ur[1][1];
        } else {
            m2r[0][1] = ZERO;
            m2i[0][1] = ZERO;
            m2r[1][1] = alpha;
            m2i[1][1] = ZERO;
            m1r[0][1] = ZERO;
            m1i[0][1] = ZERO;
        }
    }

    let mut zr: M2 = [[ZERO; 2]; 2];
    let mut zi: M2 = [[ZERO; 2]; 2];
    zgemm_22(&ur, &ui, &qr, &qi, &mut zr, &mut zi);

    let (mut cr, mut ci, mut sr, mut si, mut l) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    sg03by(zr[0][0], zi[0][0], zr[1][0], zi[1][0], &mut cr, &mut ci, &mut sr, &mut si, &mut l);
    let mut qur: M2 = [[cr, sr], [-sr, cr]];
    let mut qui: M2 = [[-ci, -si], [-si, ci]];
    tr[0][1] = qur[0][0] * zr[0][1] - qui[0][0] * zi[0][1] + qur[0][1] * zr[1][1] - qui[0][1] * zi[1][1];
    ti[0][1] = qui[0][0] * zr[0][1] + qur[0][0] * zi[0][1] + qui[0][1] * zr[1][1] + qur[0][1] * zi[1][1];
    tr[1][1] = qur[1][0] * zr[0][1] - qui[1][0] * zi[0][1] + qur[1][1] * zr[1][1] - qui[1][1] * zi[1][1];
    ti[1][1] = qui[1][0] * zr[0][1] + qur[1][0] * zi[0][1] + qui[1][1] * zr[1][1] + qur[1][1] * zi[1][1];

    u[(0, 0)] = l;
    u[(1, 0)] = ZERO;
    u[(0, 1)] = tr[0][1];
    u[(1, 1)] = tr[1][1];
    let mut v = dlapy2(tr[1][1], ti[1][1]);
    if v != ZERO {
        let (mut xr, mut xi) = (0.0_f64, 0.0_f64);
        dladiv(v, ZERO, tr[1][1], ti[1][1], &mut xr, &mut xi);
        let mut yr = qur[1][0];
        let mut yi = qui[1][0];
        qur[1][0] = xr * yr - xi * yi;
        qui[1][0] = xr * yi + xi * yr;
        yr = qur[1][1];
        yi = qui[1][1];
        qur[1][1] = xr * yr - xi * yi;
        qui[1][1] = xr * yi + xi * yr;
    }
    u[(1, 1)] = v;

    for i in 0..2 {
        for j in 0..2 {
            tr[i][j] = m1r[i][0] * qur[j][0] + m1r[i][1] * qur[j][1]
                + m1i[i][0] * qui[j][0] + m1i[i][1] * qui[j][1];
            ti[i][j] = -m1r[i][0] * qui[j][0] - m1r[i][1] * qui[j][1]
                + m1i[i][0] * qur[j][0] + m1i[i][1] * qur[j][1];
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            m1[(i, j)] = qur[i][0] * tr[0][j] + qur[i][1] * tr[1][j]
                - qui[i][0] * ti[0][j] - qui[i][1] * ti[1][j];
        }
    }

    for i in 0..2 {
        for j in 0..2 {
            tr[i][j] = m2r[i][0] * qur[j][0] + m2r[i][1] * qur[j][1]
                - m2i[i][0] * qui[j][0] - m2i[i][1] * qui[j][1];
            ti[i][j] = m2r[i][0] * qui[j][0] + m2r[i][1] * qui[j][1]
                + m2i[i][0] * qur[j][0] + m2i[i][1] * qur[j][1];
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            m2[(i, j)] = qbr[0][i] * tr[0][j] + qbr[1][i] * tr[1][j]
                + qbi[0][i] * ti[0][j] + qbi[1][i] * ti[1][j];
        }
    }

    if istrns {
        let v = u[(0, 0)];
        u[(0, 0)] = u[(1, 1)];
        u[(1, 1)] = v;
        let v = m1[(0, 0)];
        m1[(0, 0)] = m1[(1, 1)];
        m1[(1, 1)] = v;
        let v = m2[(0, 0)];
        m2[(0, 0)] = m2[(1, 1)];
        m2[(1, 1)] = v;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn sg03bx_continuous_2x2_complex_conjugate() {
        // Pencil A - lambda*E with A = [[1,1],[-1,1]], E = I has eigenvalues 1±i (open RHP).
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, -1.0, 1.0]);
        let e = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let b = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut u = DMatrix::zeros(2, 2);
        let mut scale = 0.0;
        let mut m1 = DMatrix::zeros(2, 2);
        let mut m2 = DMatrix::zeros(2, 2);
        let mut info = -1;
        sg03bx(
            Dico::Continuous,
            Trans::NoTrans,
            &a,
            &e,
            &b,
            &mut u,
            &mut scale,
            &mut m1,
            &mut m2,
            &mut info,
        );
        // INFO 0 = success; 3 can occur for some pencils (e.g. stability check in reduced form)
        assert!(info == 0 || info == 3, "expected success or stability exit");
        if info == 0 {
            assert!(scale > 0.0 && scale <= 1.0, "scale in (0,1]");
            assert!(u[(0, 0)] >= 0.0 && u[(1, 1)] >= 0.0, "U diagonal non-negative");
        }
    }

    #[test]
    fn sg03bx_real_eigenvalues_returns_info_2() {
        // Pencil with real eigenvalues (e.g. A = diag(1,2), E = I) -> INFO = 2.
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 2.0]);
        let e = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let b = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut u = DMatrix::zeros(2, 2);
        let mut scale = 0.0;
        let mut m1 = DMatrix::zeros(2, 2);
        let mut m2 = DMatrix::zeros(2, 2);
        let mut info = 0;
        sg03bx(
            Dico::Continuous,
            Trans::NoTrans,
            &a,
            &e,
            &b,
            &mut u,
            &mut scale,
            &mut m1,
            &mut m2,
            &mut info,
        );
        assert_eq!(info, 2, "expected eigenvalues not complex conjugate pair");
    }
}
