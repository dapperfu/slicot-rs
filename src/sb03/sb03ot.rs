//! SB03OT — Solve Lyapunov for Cholesky factor U (S and R block upper / upper triangular).
//! Continuous: op(S)'*X + X*op(S) = -scale^2*op(R)'*op(R). Discrete: op(S)'*X*op(S) - X = -scale^2*op(R)'*op(R).
//! Full port from SLICOT SB03OT.f.

use crate::mb04::blas::{daxpy, dcopy, dscal, dswap, dtrmm_left, dtrmm_right, dtrmv, dlarfg};
use crate::mb04::mb04nd::{mb04nd, Mb04ndUplo};
use crate::mb04::mb04od::{mb04od, Mb04odUplo};
use crate::sb03::sb03or::sb03or;
use crate::sb03::sb03oy::sb03oy_full;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

fn dlamch_p() -> f64 {
    f64::EPSILON
}
fn dlamch_s() -> f64 {
    f64::MIN_POSITIVE
}
fn dlabad(smlnum: &mut f64, bignum: &mut f64) {
    let eps = dlamch_p();
    if *smlnum <= ZERO {
        return;
    }
    *smlnum = (*smlnum / eps).min(*bignum);
    *bignum = ONE / (*smlnum);
}

/// Max norm of upper Hessenberg matrix (|a(i,j)| for i<=j+1).
fn dlanhs_max(n: usize, s: &[f64], lds: usize) -> f64 {
    let mut v = ZERO;
    for j in 0..n {
        for i in 0..=(j + 1).min(n - 1) {
            let a = s[i + j * lds].abs();
            if a > v {
                v = a;
            }
        }
    }
    v
}

/// Solve for U (Cholesky of X). S block upper Hessenberg (1×1 and 2×2 blocks), R upper triangular.
/// On exit R contains U. DWORK length >= 4*N.
pub fn sb03ot(
    discr: bool,
    ltrans: bool,
    n: usize,
    s: &[f64],
    lds: usize,
    r: &mut [f64],
    ldr: usize,
    scale: &mut f64,
    dwork: &mut [f64],
) -> i32 {
    if lds < n.max(1) || ldr < n.max(1) {
        return -5;
    }
    if dwork.len() < 4 * n {
        return -9;
    }
    *scale = ONE;
    if n == 0 {
        return 0;
    }

    let mut eps = dlamch_p();
    let mut smlnum = dlamch_s();
    let mut bignum = ONE / smlnum;
    dlabad(&mut smlnum, &mut bignum);
    smlnum = smlnum * (n * n) as f64 / eps;
    bignum = ONE / smlnum;

    let snorm = dlanhs_max(n, s, lds);
    let mut smin = smlnum.max(eps * snorm);
    let mut infom = 0_i32;
    let cont = !discr;
    let isgn = 1;

    // 2×2 local arrays (column-major: [0]= (0,0), [1]= (1,0), [2]= (0,1), [3]= (1,1))
    let mut b = [ZERO; 4];
    let mut u = [ZERO; 4];
    let mut a = [ZERO; 4];

    if !ltrans {
        // Forward: op(M) = M
        let mut kount = 0_usize;
        while kount < n {
            let k = kount;
            let tbyt = if kount >= n - 1 {
                false
            } else if s[k + 1 + k * lds] == ZERO {
                false
            } else {
                if k + 2 < n && s[k + 2 + (k + 1) * lds] != ZERO {
                    return 3;
                }
                true
            };

            if tbyt {
                kount += 2;
                // 2×2 block: B, U from S and R
                b[0] = s[k + k * lds];
                b[1] = s[k + 1 + k * lds];
                b[2] = s[k + (k + 1) * lds];
                b[3] = s[k + 1 + (k + 1) * lds];
                u[0] = r[k + k * ldr];
                u[2] = r[k + (k + 1) * ldr];
                u[3] = r[k + 1 + (k + 1) * ldr];

                let info = sb03oy_full(discr, ltrans, isgn, &mut b, 2, &mut u, 2, &mut a, 2, &mut *scale);
                if info > 1 {
                    return info;
                }
                if info != 0 {
                    infom = infom.max(info);
                }
                let scaloc = *scale;
                if scaloc != ONE {
                    for j in 0..n {
                        dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                    }
                }
                r[k + k * ldr] = u[0];
                r[k + (k + 1) * ldr] = u[2];
                r[k + 1 + (k + 1) * ldr] = u[3];

                if kount <= n {
                    let ksize = n - k - 2;
                    if ksize > 0 {
                        let k1 = ksize;
                        let k2 = ksize + k1;
                        let k3 = ksize + k2;
                        dcopy(ksize, &r[k + (k + 2) * ldr..], ldr, &mut dwork[0..], 1);
                        dcopy(ksize, &r[k + 1 + (k + 2) * ldr..], ldr, &mut dwork[k1..], 1);
                        dtrmm_right(false, ksize, 2, -ONE, &a, 2, &mut dwork[0..], ksize);
                        if cont {
                            daxpy(ksize, -r[k + k * ldr], &s[k + (k + 2) * lds..], lds, &mut dwork[0..], 1);
                            daxpy(ksize, -r[k + (k + 1) * ldr], &s[k + 1 + (k + 2) * lds..], lds, &mut dwork[0..], 1);
                            daxpy(ksize, -r[k + 1 + (k + 1) * ldr], &s[k + 1 + (k + 2) * lds..], lds, &mut dwork[k1..], 1);
                        } else {
                            daxpy(ksize, -r[k + k * ldr] * b[0], &s[k + (k + 2) * lds..], lds, &mut dwork[0..], 1);
                            daxpy(ksize, -(r[k + (k + 1) * ldr] * b[0] + r[k + 1 + (k + 1) * ldr] * b[1]), &s[k + 1 + (k + 2) * lds..], lds, &mut dwork[0..], 1);
                            daxpy(ksize, -r[k + k * ldr] * b[2], &s[k + (k + 2) * lds..], lds, &mut dwork[k1..], 1);
                            daxpy(ksize, -(r[k + (k + 1) * ldr] * b[2] + r[k + 1 + (k + 1) * ldr] * b[3]), &s[k + 1 + (k + 2) * lds..], lds, &mut dwork[k1..], 1);
                        }

                        let mut scaloc = ONE;
                        let info = sb03or(
                            discr,
                            ltrans,
                            ksize,
                            2,
                            &s[k + 2 + (k + 2) * lds..],
                            lds,
                            &b,
                            2,
                            &mut dwork[0..],
                            ksize,
                            &mut scaloc,
                        );
                        infom = infom.max(info);
                        if scaloc != ONE {
                            for j in 0..n {
                                dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                            }
                            *scale *= scaloc;
                        }
                        {
                            let (d0, d1) = dwork.split_at_mut(k2);
                            dcopy(2 * ksize, d0, 1, &mut d1[0..2 * ksize], 1);
                        }

                        if cont {
                            dswap(ksize, &mut dwork[0..], 1, &mut r[k + (k + 2) * ldr..], ldr);
                            dswap(ksize, &mut dwork[k1..], 1, &mut r[k + 1 + (k + 2) * ldr..], ldr);
                            let (d0a, d1a) = dwork.split_at_mut(k2);
                            daxpy(ksize, -a[0], &d1a[0..ksize], 1, &mut d0a[0..ksize], 1);
                            daxpy(ksize, -a[2], &d1a[ksize..2 * ksize], 1, &mut d0a[0..ksize], 1);
                            daxpy(ksize, -a[3], &d1a[ksize..2 * ksize], 1, &mut d0a[k1..], 1);
                        } else {
                            dtrmm_left(true, ksize, 2, ONE, &s[k + 2 + (k + 2) * lds..], lds, &mut dwork[0..], ksize);
                            for j in 0..ksize - 1 {
                                if s[k + 2 + j + 1 + (k + 2 + j) * lds] != ZERO {
                                    dwork[j] += s[k + 2 + j + 1 + (k + 2 + j) * lds] * dwork[k2 + j];
                                    dwork[k1 + j] += s[k + 2 + j + 1 + (k + 2 + j) * lds] * dwork[k3 + j];
                                }
                            }
                            daxpy(ksize, r[k + k * ldr], &s[k + (k + 2) * lds..], lds, &mut dwork[0..], 1);
                            daxpy(ksize, r[k + (k + 1) * ldr], &s[k + 1 + (k + 2) * lds..], lds, &mut dwork[0..], 1);
                            daxpy(ksize, r[k + 1 + (k + 1) * ldr], &s[k + 1 + (k + 2) * lds..], lds, &mut dwork[k1..], 1);
                            dswap(ksize, &mut dwork[k2..], 1, &mut r[k + (k + 2) * ldr..], ldr);
                            dswap(ksize, &mut dwork[k3..], 1, &mut r[k + 1 + (k + 2) * ldr..], ldr);

                            let mut tau1 = 0.0_f64;
                            dlarfg(3, &mut a[0], &mut b[0..2], 1, &mut tau1);
                            let v1 = b[0];
                            let t1 = tau1 * v1;
                            let v2 = b[1];
                            let t2 = tau1 * v2;
                            let sum = a[2] + v1 * b[2] + v2 * b[3];
                            b[2] -= sum * t1;
                            b[3] -= sum * t2;
                            let mut tau2 = 0.0_f64;
                            dlarfg(3, &mut a[3], &mut b[2..4], 1, &mut tau2);
                            let v3 = b[2];
                            let t3 = tau2 * v3;
                            let v4 = b[3];
                            let t4 = tau2 * v4;
                            for j in 0..ksize {
                                let sum = dwork[k2 + j] + v1 * dwork[j] + v2 * dwork[k1 + j];
                                let d1 = dwork[j] - sum * t1;
                                let d2 = dwork[k1 + j] - sum * t2;
                                let sum2 = dwork[k3 + j] + v3 * d1 + v4 * d2;
                                dwork[j] = d1 - sum2 * t3;
                                dwork[k1 + j] = d2 - sum2 * t4;
                            }
                        }

                        {
                            let (left, right) = dwork.split_at_mut(k2);
                            dcopy(ksize, &left[0..], 1, &mut right[0..ksize], 1);
                            dcopy(ksize, &left[k1..], 1, &mut right[ksize..2 * ksize], 1);
                        }
                        for j in 0..ksize {
                            dwork[2 * j + 1] = dwork[k3 + j];
                            dwork[2 * j] = dwork[k2 + j];
                        }
                        {
                            let (a_sl, rest) = dwork.split_at_mut(2 * ksize);
                            let (tau_sl, int_sl) = rest.split_at_mut(ksize);
                            let (b_dum, rest2) = int_sl.split_at_mut(1);
                            let (c_dum, work) = rest2.split_at_mut(1);
                            mb04od(
                                Mb04odUplo::Full,
                                ksize,
                                0,
                                2,
                                &mut r[k + 2 + (k + 2) * ldr..],
                                ldr,
                                a_sl,
                                2,
                                b_dum,
                                1,
                                c_dum,
                                1,
                                tau_sl,
                                work,
                            );
                        }
                    }
                }
            } else {
                kount += 1;
                // 1×1 block
                let skk = s[k + k * lds];
                let (temp, ok) = if discr {
                    let absskk = skk.abs();
                    if absskk >= ONE {
                        return 2;
                    }
                    (((ONE - absskk) * (ONE + absskk)).sqrt(), true)
                } else {
                    if skk >= ZERO {
                        return 2;
                    }
                    ((TWO * skk).abs().sqrt(), true)
                };
                let temp = if temp < smin { smin = temp; infom = 1; smin } else { temp };
                let mut scaloc = ONE;
                let dr = r[k + k * ldr].abs();
                if temp < ONE && dr > ONE && dr > bignum * temp {
                    scaloc = ONE / dr;
                }
                let alpha = temp * r[k + k * ldr].signum();
                r[k + k * ldr] /= alpha;
                if scaloc != ONE {
                    for j in 0..n {
                        dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                    }
                    *scale *= scaloc;
                }

                if kount <= n {
                    let ksize = n - k - 1;
                    if ksize > 0 {
                        let k1 = ksize;
                        let k2 = ksize + k1;
                        dcopy(ksize, &r[k + (k + 1) * ldr..], ldr, &mut dwork[0..], 1);
                        dscal(ksize, -alpha, &mut dwork[0..], 1);
                        if cont {
                            daxpy(ksize, -r[k + k * ldr], &s[k + (k + 1) * lds..], lds, &mut dwork[0..], 1);
                        } else {
                            daxpy(ksize, -skk * r[k + k * ldr], &s[k + (k + 1) * lds..], lds, &mut dwork[0..], 1);
                        }
                        let mut scaloc = ONE;
                        let info = sb03or(
                            discr,
                            ltrans,
                            ksize,
                            1,
                            &s[k + 1 + (k + 1) * lds..],
                            lds,
                            &s[k + k * lds..],
                            1,
                            &mut dwork[0..],
                            ksize,
                            &mut scaloc,
                        );
                        infom = infom.max(info);
                        if scaloc != ONE {
                            for j in 0..n {
                                dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                            }
                            *scale *= scaloc;
                        }
                        {
                            let (left, right) = dwork.split_at_mut(k1);
                            dcopy(ksize, &left[0..], 1, &mut right[0..ksize], 1);
                        }
                        dswap(ksize, &mut dwork[0..], 1, &mut r[k + (k + 1) * ldr..], ldr);

                        if cont {
                            let (left, right) = dwork.split_at_mut(k1);
                            daxpy(ksize, -alpha, &right[0..ksize], 1, &mut left[0..ksize], 1);
                        } else {
                            dscal(ksize, -skk, &mut dwork[0..], 1);
                            daxpy(ksize, alpha * r[k + k * ldr], &s[k + (k + 1) * lds..], lds, &mut dwork[0..], 1);
                            for j in 0..ksize - 1 {
                                if s[k + 1 + j + 1 + (k + 1 + j) * lds] != ZERO {
                                    dwork[j] += alpha * s[k + 1 + j + 1 + (k + 1 + j) * lds] * dwork[k1 + j];
                                }
                            }
                            dtrmv(true, true, ksize, &s[k + 1 + (k + 1) * lds..], lds, &mut dwork[k1..], 1);
                            let (left, right) = dwork.split_at_mut(k1);
                            daxpy(ksize, alpha, &right[0..ksize], 1, &mut left[0..ksize], 1);
                        }
                        {
                            let (a_sl, rest) = dwork.split_at_mut(ksize);
                            let (tau_sl, int_sl) = rest.split_at_mut(ksize);
                            let (b_dum, rest2) = int_sl.split_at_mut(1);
                            let (c_dum, work) = rest2.split_at_mut(1);
                            mb04od(
                                Mb04odUplo::Full,
                                ksize,
                                0,
                                1,
                                &mut r[k + 1 + (k + 1) * ldr..],
                                ldr,
                                a_sl,
                                1,
                                b_dum,
                                1,
                                c_dum,
                                1,
                                tau_sl,
                                work,
                            );
                        }
                    }
                }
            }
        }
    } else {
        // Backward: op(M) = M'
        let mut kount = n;
        while kount >= 1 {
            let k = kount - 1;
            let (tbyt, next_kount) = if kount == 1 {
                (false, 0)
            } else if s[k + (k - 1) * lds] == ZERO {
                (false, kount - 1)
            } else {
                let k_prev = k - 1;
                if k_prev > 0 && s[k_prev + (k_prev - 1) * lds] != ZERO {
                    return 3;
                }
                (true, kount.saturating_sub(2))
            };

            if tbyt {
                kount = next_kount;
                let k = k - 1; // block starts at k (0-based)
                b[0] = s[k + k * lds];
                b[1] = s[k + 1 + k * lds];
                b[2] = s[k + (k + 1) * lds];
                b[3] = s[k + 1 + (k + 1) * lds];
                u[0] = r[k + k * ldr];
                u[2] = r[k + (k + 1) * ldr];
                u[3] = r[k + 1 + (k + 1) * ldr];

                let info = sb03oy_full(discr, ltrans, isgn, &mut b, 2, &mut u, 2, &mut a, 2, &mut *scale);
                if info > 1 {
                    return info;
                }
                if info != 0 {
                    infom = infom.max(info);
                }
                let scaloc = *scale;
                if scaloc != ONE {
                    for j in 0..n {
                        dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                    }
                }
                r[k + k * ldr] = u[0];
                r[k + (k + 1) * ldr] = u[2];
                r[k + 1 + (k + 1) * ldr] = u[3];

                if kount >= 1 {
                    let ksize = k;
                    if ksize > 0 {
                        let k1 = ksize;
                        let k2 = ksize + k1;
                        let k3 = ksize + k2;
                        dcopy(ksize, &r[k * ldr..], 1, &mut dwork[0..], 1);
                        dcopy(ksize, &r[(k + 1) * ldr..], 1, &mut dwork[k1..], 1);
                        dtrmm_right(true, ksize, 2, -ONE, &a, 2, &mut dwork[0..], ksize);
                        if cont {
                            daxpy(ksize, -r[k + k * ldr], &s[k * lds..], 1, &mut dwork[0..], 1);
                            daxpy(ksize, -r[k + (k + 1) * ldr], &s[k * lds..], 1, &mut dwork[k1..], 1);
                            daxpy(ksize, -r[k + 1 + (k + 1) * ldr], &s[(k + 1) * lds..], 1, &mut dwork[k1..], 1);
                        } else {
                            daxpy(ksize, -(r[k + k * ldr] * b[0] + r[k + (k + 1) * ldr] * b[2]), &s[k * lds..], 1, &mut dwork[0..], 1);
                            daxpy(ksize, -r[k + 1 + (k + 1) * ldr] * b[2], &s[(k + 1) * lds..], 1, &mut dwork[0..], 1);
                            daxpy(ksize, -(r[k + k * ldr] * b[1] + r[k + (k + 1) * ldr] * b[3]), &s[k * lds..], 1, &mut dwork[k1..], 1);
                            daxpy(ksize, -r[k + 1 + (k + 1) * ldr] * b[3], &s[(k + 1) * lds..], 1, &mut dwork[k1..], 1);
                        }
                        let mut scaloc = ONE;
                        let info = sb03or(
                            discr,
                            ltrans,
                            ksize,
                            2,
                            &s[0..],
                            lds,
                            &b,
                            2,
                            &mut dwork[0..],
                            ksize,
                            &mut scaloc,
                        );
                        infom = infom.max(info);
                        if scaloc != ONE {
                            for j in 0..n {
                                dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                            }
                            *scale *= scaloc;
                        }
                        {
                            let (d0, d1) = dwork.split_at_mut(k2);
                            dcopy(2 * ksize, d0, 1, &mut d1[0..2 * ksize], 1);
                        }

                        if cont {
                            dswap(ksize, &mut dwork[0..], 1, &mut r[k * ldr..], 1);
                            dswap(ksize, &mut dwork[k1..], 1, &mut r[(k + 1) * ldr..], 1);
                            let (d0a, d1a) = dwork.split_at_mut(k2);
                            daxpy(ksize, -a[0], &d1a[0..ksize], 1, &mut d0a[0..ksize], 1);
                            daxpy(ksize, -a[2], &d1a[0..ksize], 1, &mut d0a[k1..], 1);
                            daxpy(ksize, -a[3], &d1a[ksize..2 * ksize], 1, &mut d0a[k1..], 1);
                        } else {
                            dtrmm_left(false, ksize, 2, ONE, &s[0..], lds, &mut dwork[0..], ksize);
                            for j in 1..ksize {
                                if s[j + (j - 1) * lds] != ZERO {
                                    dwork[j] += s[j + (j - 1) * lds] * dwork[k2 + j - 1];
                                    dwork[k1 + j] += s[j + (j - 1) * lds] * dwork[k3 + j - 1];
                                }
                            }
                            daxpy(ksize, r[k + k * ldr], &s[k * lds..], 1, &mut dwork[0..], 1);
                            daxpy(ksize, r[k + (k + 1) * ldr], &s[k * lds..], 1, &mut dwork[k1..], 1);
                            daxpy(ksize, r[k + 1 + (k + 1) * ldr], &s[(k + 1) * lds..], 1, &mut dwork[k1..], 1);
                            dswap(ksize, &mut dwork[k2..], 1, &mut r[k * ldr..], 1);
                            dswap(ksize, &mut dwork[k3..], 1, &mut r[(k + 1) * ldr..], 1);

                            let mut tau1 = 0.0_f64;
                            dlarfg(3, &mut a[3], &mut b[1..=3], 2, &mut tau1);
                            let v1 = b[1];
                            let t1 = tau1 * v1;
                            let v2 = b[3];
                            let t2 = tau1 * v2;
                            let sum = a[2] + v1 * b[0] + v2 * b[2];
                            b[0] -= sum * t1;
                            b[2] -= sum * t2;
                            let mut tau2 = 0.0_f64;
                            dlarfg(3, &mut a[0], &mut b[0..=2], 2, &mut tau2);
                            let v3 = b[0];
                            let t3 = tau2 * v3;
                            let v4 = b[2];
                            let t4 = tau2 * v4;
                            for j in 0..ksize {
                                let sum = dwork[k3 + j] + v1 * dwork[j] + v2 * dwork[k1 + j];
                                let d1 = dwork[j] - sum * t1;
                                let d2 = dwork[k1 + j] - sum * t2;
                                let sum2 = dwork[k2 + j] + v3 * d1 + v4 * d2;
                                dwork[j] = d1 - sum2 * t3;
                                dwork[k1 + j] = d2 - sum2 * t4;
                            }
                        }
                        {
                            let (a_sl, rest) = dwork.split_at_mut(2 * ksize);
                            let (tau_sl, int_sl) = rest.split_at_mut(ksize);
                            let (b_dum, rest2) = int_sl.split_at_mut(1);
                            let (c_dum, work) = rest2.split_at_mut(1);
                            mb04nd(
                                Mb04ndUplo::Full,
                                ksize,
                                0,
                                2,
                                &mut r[0..],
                                ldr,
                                a_sl,
                                ksize,
                                b_dum,
                                1,
                                c_dum,
                                1,
                                tau_sl,
                                work,
                            );
                        }
                    }
                }
            } else {
                kount = next_kount;
                let skk = s[k + k * lds];
                let (temp, _) = if discr {
                    let absskk = skk.abs();
                    if absskk >= ONE {
                        return 2;
                    }
                    (((ONE - absskk) * (ONE + absskk)).sqrt(), true)
                } else {
                    if skk >= ZERO {
                        return 2;
                    }
                    ((TWO * skk).abs().sqrt(), true)
                };
                let temp = if temp < smin { smin = temp; infom = 1; smin } else { temp };
                let mut scaloc = ONE;
                let dr = r[k + k * ldr].abs();
                if temp < ONE && dr > ONE && dr > bignum * temp {
                    scaloc = ONE / dr;
                }
                let alpha = temp * r[k + k * ldr].signum();
                r[k + k * ldr] /= alpha;
                if scaloc != ONE {
                    for j in 0..n {
                        dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                    }
                    *scale *= scaloc;
                }

                if kount >= 1 {
                    let ksize = k;
                    if ksize > 0 {
                        let k1 = ksize;
                        let k2 = ksize + k1;
                        dcopy(ksize, &r[k * ldr..], 1, &mut dwork[0..], 1);
                        dscal(ksize, -alpha, &mut dwork[0..], 1);
                        if cont {
                            daxpy(ksize, -r[k + k * ldr], &s[k * lds..], 1, &mut dwork[0..], 1);
                        } else {
                            daxpy(ksize, -skk * r[k + k * ldr], &s[k * lds..], 1, &mut dwork[0..], 1);
                        }
                        let mut scaloc = ONE;
                        let info = sb03or(
                            discr,
                            ltrans,
                            ksize,
                            1,
                            &s[0..],
                            lds,
                            &s[k + k * lds..],
                            1,
                            &mut dwork[0..],
                            ksize,
                            &mut scaloc,
                        );
                        infom = infom.max(info);
                        if scaloc != ONE {
                            for j in 0..n {
                                dscal(j + 1, scaloc, &mut r[j * ldr..], 1);
                            }
                            *scale *= scaloc;
                        }
                        {
                            let (left, right) = dwork.split_at_mut(k1);
                            dcopy(ksize, &left[0..], 1, &mut right[0..ksize], 1);
                        }
                        dswap(ksize, &mut dwork[0..], 1, &mut r[k * ldr..], 1);

                        if cont {
                            let (left, right) = dwork.split_at_mut(k1);
                            daxpy(ksize, -alpha, &right[0..ksize], 1, &mut left[0..ksize], 1);
                        } else {
                            dscal(ksize, -skk, &mut dwork[0..], 1);
                            daxpy(ksize, alpha * r[k + k * ldr], &s[k * lds..], 1, &mut dwork[0..], 1);
                            for j in 0..ksize - 1 {
                                if s[j + 1 + j * lds] != ZERO {
                                    dwork[j] += alpha * s[j + 1 + j * lds] * dwork[k1 + j];
                                }
                            }
                            dtrmv(true, false, ksize, &s[0..], lds, &mut dwork[k1..], 1);
                            let (left, right) = dwork.split_at_mut(k1);
                            daxpy(ksize, alpha, &right[0..ksize], 1, &mut left[0..ksize], 1);
                        }
                        {
                            let (a_sl, rest) = dwork.split_at_mut(ksize);
                            let (tau_sl, int_sl) = rest.split_at_mut(ksize);
                            let (b_dum, rest2) = int_sl.split_at_mut(1);
                            let (c_dum, work) = rest2.split_at_mut(1);
                            mb04nd(
                                Mb04ndUplo::Full,
                                ksize,
                                0,
                                1,
                                &mut r[0..],
                                ldr,
                                a_sl,
                                ksize,
                                b_dum,
                                1,
                                c_dum,
                                1,
                                tau_sl,
                                work,
                            );
                        }
                    }
                }
            }
        }
    }

    infom
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03ot_1x1_cont() {
        let n = 1;
        let s = [-1.0];
        let lds = 1;
        let mut r = [1.0];
        let ldr = 1;
        let mut scale = 0.0;
        let mut dwork = vec![0.0; 4 * n];
        let info = sb03ot(false, false, n, &s, lds, &mut r, ldr, &mut scale, &mut dwork);
        assert_eq!(info, 0);
        assert!((scale - 1.0).abs() < 1e-10);
        assert!(r[0].abs() > 0.0);
    }
}
