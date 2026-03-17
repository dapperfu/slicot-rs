//! SB03MX — Solve discrete Lyapunov op(A)'*X*op(A) - X = scale*C with A upper quasi-triangular.
//! C is symmetric; solution X overwrites C. Bartels-Stewart with 1×1 and 2×2 blocks.

use super::sb03mv;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Solve TL'*X*TR - X = B for 2×2 X (SB04PX with LTRANL=true, LTRANR=false, ISGN=-1).
/// Column-major tl(ldtl,*), tr(ldtr,*), b(ldb,*), x(ldx,*).
fn solve_discrete_sylvester_2x2(
    tl: &[f64], ldtl: usize,
    tr: &[f64], ldtr: usize,
    b: &[f64], ldb: usize,
    x: &mut [f64], ldx: usize,
    scaloc: &mut f64,
    smin: f64,
    info: &mut i32,
) {
    let tl11 = tl[0]; let tl12 = tl[ldtl]; let tl21 = tl[1]; let tl22 = tl[1 + ldtl];
    let tr11 = tr[0]; let tr12 = tr[ldtr]; let tr21 = tr[1]; let tr22 = tr[1 + ldtr];
    // (TR' kron TL' - I)*vec(X) = vec(B), vec(X)=[x11,x21,x12,x22]
    let mut m = [
        [tl11 * tr11 - ONE, tl21 * tr11, tl11 * tr21, tl21 * tr21],
        [tl12 * tr11, tl22 * tr11 - ONE, tl12 * tr21, tl22 * tr21],
        [tl11 * tr12, tl21 * tr12, tl11 * tr22 - ONE, tl21 * tr22],
        [tl12 * tr12, tl22 * tr12, tl12 * tr22, tl22 * tr22 - ONE],
    ];
    let mut rhs = [b[0], b[1], b[ldb], b[1 + ldb]];
    for i in 0..4 {
        if m[i][i].abs() < smin {
            *info = 1;
            m[i][i] = smin;
        }
        for j in (i + 1)..4 {
            let mult = m[j][i] / m[i][i];
            m[j][i] = mult;
            rhs[j] -= mult * rhs[i];
            for k in (i + 1)..4 {
                m[j][k] -= mult * m[i][k];
            }
        }
    }
    if m[3][3].abs() < smin {
        *info = 1;
        m[3][3] = smin;
    }
    let mut t = [ZERO; 4];
    for k in (0..4).rev() {
        t[k] = rhs[k];
        for j in (k + 1)..4 {
            t[k] -= m[k][j] * t[j];
        }
        t[k] /= m[k][k];
    }
    x[0] = t[0];
    x[1] = t[1];
    x[ldx] = t[2];
    x[1 + ldx] = t[3];
    *scaloc = ONE;
}

/// DLALN2-style: solve (A'*X*B + beta*X = C) for 2x1 X with A 2x2, B 1x1 => (A'*b - beta)*x = c.
/// Here we need A(K,K)'*x*B(L,L) - x = vec => (A(K,K)'*B(L,L) - I)*x = vec. So 2x2 system.
fn solve_2x1_discrete(smin: f64, a: &[f64], lda: usize, b_scalar: f64, vec: &[f64], x: &mut [f64], scaloc: &mut f64, info: &mut i32) {
    let a11 = a[0]; let a12 = a[lda]; let a21 = a[1]; let a22 = a[1 + lda];
    let m11 = a11 * b_scalar - ONE;
    let m12 = a21 * b_scalar;
    let m21 = a12 * b_scalar;
    let m22 = a22 * b_scalar - ONE;
    let det = m11 * m22 - m12 * m21;
    if det.abs() <= smin {
        *info = 1;
        *scaloc = ONE;
        let sd = smin;
        x[0] = (m22 * vec[0] - m12 * vec[1]) / sd;
        x[1] = (-m21 * vec[0] + m11 * vec[1]) / sd;
        return;
    }
    *scaloc = ONE;
    x[0] = (m22 * vec[0] - m12 * vec[1]) / det;
    x[1] = (-m21 * vec[0] + m11 * vec[1]) / det;
}

/// TRANA: 'N' => op(A)=A, 'T' or 'C' => op(A)=A'. A, C column-major. dwork length >= 2*n.
pub fn sb03mx(
    trana: char,
    n: usize,
    a: &[f64],
    lda: usize,
    c: &mut [f64],
    ldc: usize,
    scale: &mut f64,
    dwork: &mut [f64],
    info: &mut i32,
) {
    *info = 0;
    *scale = ONE;
    if n == 0 {
        return;
    }
    let notrana = matches!(trana, 'N' | 'n');
    let eps = f64::EPSILON;
    let smlnum = (f64::MIN_POSITIVE / eps).min(1e-10);
    let bignum = ONE / smlnum;
    let smlnum = smlnum * (n * n) as f64 / eps;
    let bignum = ONE / smlnum;
    let anorm = dlanhs_max(n, a, lda, dwork);
    let smin = smlnum.max(eps * anorm);
    let np1 = n + 1;

    let at = |i: usize, j: usize| -> f64 {
        if notrana {
            a[i + j * lda]
        } else {
            a[j + i * lda]
        }
    };

    if notrana {
        // Solve A'*X*A - X = scale*C. Upper-left order.
        let mut l = 0;
        while l < n {
            let l1 = l;
            let l2 = if l < n - 1 && a[l + 1 + l * lda].abs() > ZERO { l + 1 } else { l };
            l = l2 + 1;

            dwork[l1] = ZERO;
            if l2 > l1 {
                dwork[np1 + l1] = ZERO;
            }
            dsymv_lower(l1, c, ldc, &a[l1 * lda..], 1, &mut dwork[..n]);
            if l2 > l1 {
                dsymv_lower(l1, c, ldc, &a[l2 * lda..], 1, &mut dwork[np1..]);
            }

            let mut k = l1;
            while k < n {
                let k1 = k;
                let k2 = if k < n - 1 && a[k + 1 + k * lda].abs() > ZERO { k + 1 } else { k };
                k = k2 + 1;

                if l2 > l1 {
                    dwork[k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l1 * lda]).sum();
                    if k1 != k2 {
                        dwork[k2] = (0..l1).map(|i| c[k2 + i * ldc] * a[i + l1 * lda]).sum();
                    }
                    dwork[np1 + k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l2 * lda]).sum();
                    if k1 != k2 {
                        dwork[np1 + k2] = (0..l1).map(|i| c[k2 + i * ldc] * a[i + l2 * lda]).sum();
                    }
                }

                if l1 == l2 && k1 == k2 {
                    dwork[k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l1 * lda]).sum();
                    let mut rhs = c[k1 + l1 * ldc]
                        - ((0..=k1).map(|i| a[i + k1 * lda] * dwork[i]).sum::<f64>()
                            + at(l1, l1) * (0..k1).map(|i| a[i + k1 * lda] * c[i + l1 * ldc]).sum::<f64>());
                    let a11 = at(k1, k1) * at(l1, l1) - ONE;
                    let a11 = if a11.abs() <= smin { smin } else { a11 };
                    if a11.abs() <= smin {
                        *info = 1;
                    }
                    let mut scaloc = ONE;
                    if a11.abs() < ONE && rhs.abs() > ONE && rhs.abs() > bignum * a11.abs() {
                        scaloc = ONE / rhs.abs();
                        rhs *= scaloc;
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    let x11 = rhs / a11;
                    c[k1 + l1 * ldc] = x11;
                    if k1 != l1 {
                        c[l1 + k1 * ldc] = x11;
                    }
                } else if l1 == l2 && k1 != k2 {
                    dwork[k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l1 * lda]).sum();
                    dwork[k2] = (0..l1).map(|i| c[k2 + i * ldc] * a[i + l1 * lda]).sum();
                    let mut vec = [
                        c[k1 + l1 * ldc]
                            - ((0..=k2).map(|i| a[i + k1 * lda] * dwork[i]).sum::<f64>()
                                + at(l1, l1) * (0..k1).map(|i| a[i + k1 * lda] * c[i + l1 * ldc]).sum::<f64>()),
                        c[k2 + l1 * ldc]
                            - ((0..=k2).map(|i| a[i + k2 * lda] * dwork[i]).sum::<f64>()
                                + at(l1, l1) * (0..k1).map(|i| a[i + k2 * lda] * c[i + l1 * ldc]).sum::<f64>()),
                    ];
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_discrete(smin, &a[k1 + k1 * lda..], lda, at(l1, l1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k2 + l1 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l1 + k2 * ldc] = x[1];
                } else if l1 != l2 && k1 == k2 {
                    let p11: f64 = (0..k1).map(|i| a[i + k1 * lda] * c[i + l1 * ldc]).sum();
                    let p12: f64 = (0..k1).map(|i| a[i + k1 * lda] * c[i + l2 * ldc]).sum();
                    dwork[k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l1 * lda]).sum();
                    dwork[np1 + k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l2 * lda]).sum();
                    let mut vec = [
                        c[k1 + l1 * ldc]
                            - ((0..=k1).map(|i| a[i + k1 * lda] * dwork[i]).sum::<f64>()
                                + p11 * at(l1, l1) + p12 * at(l2, l1)),
                        c[k1 + l2 * ldc]
                            - ((0..=k1).map(|i| a[i + k1 * lda] * dwork[np1 + i]).sum::<f64>()
                                + p11 * at(l1, l2) + p12 * at(l2, l2)),
                    ];
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_discrete(smin, &a[k1 + k1 * lda..], lda, at(l1, l1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k1 + l2 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l2 + k1 * ldc] = x[1];
                } else {
                    let p11: f64 = (0..k1).map(|i| a[i + k1 * lda] * c[i + l1 * ldc]).sum();
                    let p12: f64 = (0..k1).map(|i| a[i + k1 * lda] * c[i + l2 * ldc]).sum();
                    let p21: f64 = (0..k1).map(|i| a[i + k2 * lda] * c[i + l1 * ldc]).sum();
                    let p22: f64 = (0..k1).map(|i| a[i + k2 * lda] * c[i + l2 * ldc]).sum();
                    dwork[k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l1 * lda]).sum();
                    dwork[k2] = (0..l1).map(|i| c[k2 + i * ldc] * a[i + l1 * lda]).sum();
                    dwork[np1 + k1] = (0..l1).map(|i| c[k1 + i * ldc] * a[i + l2 * lda]).sum();
                    dwork[np1 + k2] = (0..l1).map(|i| c[k2 + i * ldc] * a[i + l2 * lda]).sum();
                    let mut vec = [
                        c[k1 + l1 * ldc]
                            - ((0..=k2).map(|i| a[i + k1 * lda] * dwork[i]).sum::<f64>()
                                + p11 * at(l1, l1) + p12 * at(l2, l1)),
                        c[k1 + l2 * ldc]
                            - ((0..=k2).map(|i| a[i + k1 * lda] * dwork[np1 + i]).sum::<f64>()
                                + p11 * at(l1, l2) + p12 * at(l2, l2)),
                        c[k2 + l1 * ldc]
                            - ((0..=k2).map(|i| a[i + k2 * lda] * dwork[i]).sum::<f64>()
                                + p21 * at(l1, l1) + p22 * at(l2, l1)),
                        c[k2 + l2 * ldc]
                            - ((0..=k2).map(|i| a[i + k2 * lda] * dwork[np1 + i]).sum::<f64>()
                                + p21 * at(l1, l2) + p22 * at(l2, l2)),
                    ];
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut xnorm = ZERO;
                    let mut ierr = 0;
                    if k1 == l1 {
                        let mut t = [a[k1 + k1 * lda], a[k1 + 1 + k1 * lda], a[k1 + (k1 + 1) * lda], a[k1 + 1 + (k1 + 1) * lda]];
                        let mut b = [vec[0], vec[2], vec[1], vec[3]];
                        sb03mv::sb03mv(false, true, &t, 2, &b, 2, &mut scaloc, &mut x, 2, &mut xnorm, &mut ierr);
                        x[1] = x[2];
                    } else {
                        solve_discrete_sylvester_2x2(
                            &a[k1 + k1 * lda..], lda,
                            &a[l1 + l1 * lda..], lda,
                            &[vec[0], vec[2], vec[1], vec[3]], 2,
                            &mut x, 2,
                            &mut scaloc, smin, &mut ierr,
                        );
                    }
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k1 + l2 * ldc] = x[2];
                    c[k2 + l1 * ldc] = x[1];
                    c[k2 + l2 * ldc] = x[3];
                    if k1 != l1 {
                        c[l1 + k1 * ldc] = x[0];
                        c[l2 + k1 * ldc] = x[2];
                        c[l1 + k2 * ldc] = x[1];
                        c[l2 + k2 * ldc] = x[3];
                    }
                }
            }
        }
    } else {
        // Solve A*X*A' - X = scale*C. Bottom-right order.
        let mut l = n;
        while l > 0 {
            l -= 1;
            let l2 = l;
            let l1 = if l > 0 && a[l + (l - 1) * lda].abs() > ZERO { l - 1 } else { l };
            let minl1n = (l1 + 1).min(n);
            let minl2n = (l2 + 1).min(n);

            if l2 < n - 1 {
                dsymv_upper(n - l2 - 1, c, ldc, l2 + 1, &a[l1 + (l2 + 1) * lda..], lda, &mut dwork[l2 + 1..]);
                dsymv_upper(n - l2 - 1, c, ldc, l2 + 1, &a[l2 + (l2 + 1) * lda..], lda, &mut dwork[np1 + l2..]);
            }

            let mut k = l2;
            loop {
                let k2 = k;
                let k1 = if k > 0 && a[k + (k - 1) * lda].abs() > ZERO { k - 1 } else { k };
                let mink1n = (k1 + 1).min(n);
                let mink2n = (k2 + 1).min(n);

                if l1 == l2 && k1 == k2 {
                    dwork[k1] = (minl1n..n).map(|j| c[k1 + j * ldc] * a[l1 + j * lda]).sum();
                    let mut rhs = c[k1 + l1 * ldc]
                        - ((k1..n).map(|j| a[k1 + j * lda] * dwork[j]).sum::<f64>()
                            + (mink1n..n).map(|j| c[j + l1 * ldc] * a[k1 + j * lda]).sum::<f64>() * at(l1, l1));
                    let a11 = at(k1, k1) * at(l1, l1) - ONE;
                    let a11 = if a11.abs() <= smin { smin } else { a11 };
                    if a11.abs() <= smin {
                        *info = 1;
                    }
                    let mut scaloc = ONE;
                    if a11.abs() < ONE && rhs.abs() > ONE && rhs.abs() > bignum * a11.abs() {
                        scaloc = ONE / rhs.abs();
                        rhs *= scaloc;
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    let x11 = rhs / a11;
                    c[k1 + l1 * ldc] = x11;
                    if k1 != l1 {
                        c[l1 + k1 * ldc] = x11;
                    }
                } else if l1 == l2 && k1 != k2 {
                    dwork[k1] = (minl1n..n).map(|j| c[k1 + j * ldc] * a[l1 + j * lda]).sum();
                    dwork[k2] = (minl1n..n).map(|j| c[k2 + j * ldc] * a[l1 + j * lda]).sum();
                    let mut vec = [
                        c[k1 + l1 * ldc]
                            - ((k1..n).map(|j| a[k1 + j * lda] * dwork[j]).sum::<f64>()
                                + (mink2n..n).map(|j| c[j + l1 * ldc] * a[k1 + j * lda]).sum::<f64>() * at(l1, l1)),
                        c[k2 + l1 * ldc]
                            - ((k1..n).map(|j| a[k2 + j * lda] * dwork[j]).sum::<f64>()
                                + (mink2n..n).map(|j| c[j + l1 * ldc] * a[k2 + j * lda]).sum::<f64>() * at(l1, l1)),
                    ];
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_discrete(smin, &a[k1 + k1 * lda..], lda, at(l1, l1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k2 + l1 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l1 + k2 * ldc] = x[1];
                } else if l1 != l2 && k1 == k2 {
                    let p11: f64 = (mink1n..n).map(|j| a[k1 + j * lda] * c[j + l1 * ldc]).sum();
                    let p12: f64 = (mink1n..n).map(|j| a[k1 + j * lda] * c[j + l2 * ldc]).sum();
                    dwork[k1] = (minl2n..n).map(|j| c[k1 + j * ldc] * a[l1 + j * lda]).sum();
                    dwork[np1 + k1] = (minl2n..n).map(|j| c[k1 + j * ldc] * a[l2 + j * lda]).sum();
                    let mut vec = [
                        c[k1 + l1 * ldc]
                            - ((k1..n).map(|j| a[k1 + j * lda] * dwork[j]).sum::<f64>()
                                + p11 * at(l1, l1) + p12 * at(l1, l2)),
                        c[k1 + l2 * ldc]
                            - ((k1..n).map(|j| a[k1 + j * lda] * dwork[np1 + j]).sum::<f64>()
                                + p11 * at(l2, l1) + p12 * at(l2, l2)),
                    ];
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_discrete(smin, &a[k1 + k1 * lda..], lda, at(l1, l1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k1 + l2 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l2 + k1 * ldc] = x[1];
                } else {
                    let p11: f64 = (mink2n..n).map(|j| a[k1 + j * lda] * c[j + l1 * ldc]).sum();
                    let p12: f64 = (mink2n..n).map(|j| a[k1 + j * lda] * c[j + l2 * ldc]).sum();
                    let p21: f64 = (mink2n..n).map(|j| a[k2 + j * lda] * c[j + l1 * ldc]).sum();
                    let p22: f64 = (mink2n..n).map(|j| a[k2 + j * lda] * c[j + l2 * ldc]).sum();
                    dwork[k1] = (minl2n..n).map(|j| c[k1 + j * ldc] * a[l1 + j * lda]).sum();
                    dwork[k2] = (minl2n..n).map(|j| c[k2 + j * ldc] * a[l1 + j * lda]).sum();
                    dwork[np1 + k1] = (minl2n..n).map(|j| c[k1 + j * ldc] * a[l2 + j * lda]).sum();
                    dwork[np1 + k2] = (minl2n..n).map(|j| c[k2 + j * ldc] * a[l2 + j * lda]).sum();
                    let mut vec = [
                        c[k1 + l1 * ldc]
                            - ((k1..n).map(|j| a[k1 + j * lda] * dwork[j]).sum::<f64>()
                                + p11 * at(l1, l1) + p12 * at(l1, l2)),
                        c[k1 + l2 * ldc]
                            - ((k1..n).map(|j| a[k1 + j * lda] * dwork[np1 + j]).sum::<f64>()
                                + p11 * at(l2, l1) + p12 * at(l2, l2)),
                        c[k2 + l1 * ldc]
                            - ((k1..n).map(|j| a[k2 + j * lda] * dwork[j]).sum::<f64>()
                                + p21 * at(l1, l1) + p22 * at(l1, l2)),
                        c[k2 + l2 * ldc]
                            - ((k1..n).map(|j| a[k2 + j * lda] * dwork[np1 + j]).sum::<f64>()
                                + p21 * at(l2, l1) + p22 * at(l2, l2)),
                    ];
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut xnorm = ZERO;
                    let mut ierr = 0;
                    if k1 == l1 {
                        let mut t = [a[k1 + k1 * lda], a[k1 + 1 + k1 * lda], a[k1 + (k1 + 1) * lda], a[k1 + 1 + (k1 + 1) * lda]];
                        let mut b = [vec[0], vec[2], vec[1], vec[3]];
                        sb03mv::sb03mv(true, true, &t, 2, &b, 2, &mut scaloc, &mut x, 2, &mut xnorm, &mut ierr);
                        x[1] = x[2];
                    } else {
                        solve_discrete_sylvester_2x2(
                            &a[l1 + l1 * lda..], lda,
                            &a[k1 + k1 * lda..], lda,
                            &[vec[0], vec[2], vec[1], vec[3]], 2,
                            &mut x, 2,
                            &mut scaloc, smin, &mut ierr,
                        );
                    }
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        for i in 0..n {
                            dwork[i] *= scaloc;
                            dwork[np1 + i] *= scaloc;
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k1 + l2 * ldc] = x[2];
                    c[k2 + l1 * ldc] = x[1];
                    c[k2 + l2 * ldc] = x[3];
                    if k1 != l1 {
                        c[l1 + k1 * ldc] = x[0];
                        c[l2 + k1 * ldc] = x[2];
                        c[l1 + k2 * ldc] = x[1];
                        c[l2 + k2 * ldc] = x[3];
                    }
                }

                if k1 == 0 {
                    break;
                }
                k = k1 - 1;
            }
        }
    }
}

fn dlanhs_max(n: usize, a: &[f64], lda: usize, dwork: &mut [f64]) -> f64 {
    let mut anorm = 0.0;
    for i in 0..n {
        for j in 0..=i.min(n - 1) {
            anorm += a[i + j * lda].abs();
        }
    }
    anorm
}

/// DSYMV Lower: y = A*x where A is lower part of C (column-major), x is first column of a_vec
fn dsymv_lower(n: usize, c: &[f64], ldc: usize, a_vec: &[f64], _inc: usize, y: &mut [f64]) {
    for i in 0..n {
        y[i] = (0..=i).map(|j| c[i + j * ldc] * a_vec[j]).sum();
    }
}

/// DSYMV Upper: y(off:) = C(off:,off:)*x
fn dsymv_upper(size: usize, c: &[f64], ldc: usize, off: usize, a_row: &[f64], lda: usize, y: &mut [f64]) {
    for i in 0..size {
        let row = off + i;
        y[i] = (0..size).map(|j| c[row + (off + j) * ldc] * a_row[j * lda]).sum();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03mx_1x1() {
        let n = 1;
        let mut a = [0.5];
        let mut c = [1.0];
        let mut scale = 0.0;
        let mut dwork = [0.0; 4];
        let mut info = 0;
        sb03mx('N', n, &a, 1, &mut c, 1, &mut scale, &mut dwork, &mut info);
        assert_eq!(info, 0);
        assert!((c[0] - (-4.0 / 3.0)).abs() < 1e-10); // 0.25*X - X = 1 => X = -4/3
    }
}
