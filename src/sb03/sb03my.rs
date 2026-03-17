//! SB03MY — Solve continuous Lyapunov op(A)'*X + X*op(A) = scale*C with A upper quasi-triangular.
//! C is symmetric; solution X overwrites C. Bartels-Stewart with 1×1 and 2×2 blocks.

use super::sb03mw;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Solve (A' - beta*I)*x = b for 2×1 x (LAPACK DLALN2 style, trans=true).
fn solve_2x1_trans(smin: f64, a: &[f64], lda: usize, beta: f64, b: &[f64], x: &mut [f64], scaloc: &mut f64, info: &mut i32) {
    let a11 = a[0];
    let a12 = a[lda];
    let a21 = a[1];
    let a22 = a[1 + lda];
    let m = [[a11 - beta, a21], [a12, a22 - beta]];
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    if det.abs() <= smin {
        *info = 1;
        *scaloc = ONE;
        let scale_det = smin;
        x[0] = (m[1][1] * b[0] - m[0][1] * b[1]) / scale_det;
        x[1] = (-m[1][0] * b[0] + m[0][0] * b[1]) / scale_det;
        return;
    }
    *scaloc = ONE;
    x[0] = (m[1][1] * b[0] - m[0][1] * b[1]) / det;
    x[1] = (-m[1][0] * b[0] + m[0][0] * b[1]) / det;
}

/// Solve 2×2 Sylvester A'*X + X*B = C (op(A)=A', op(B)=B). 4×4 system in vec(X).
fn solve_sylvester_2x2(
    a: &[f64], lda: usize,
    b: &[f64], ldb: usize,
    c: &[f64], ldc: usize,
    x: &mut [f64], ldx: usize,
    scaloc: &mut f64,
    smin: f64,
    info: &mut i32,
) {
    let a11 = a[0]; let a12 = a[lda]; let a21 = a[1]; let a22 = a[1 + lda];
    let b11 = b[0]; let b12 = b[ldb]; let b21 = b[1]; let b22 = b[1 + ldb];
    // vec(X) = [x11, x21, x12, x22]. (A'*X + X*B)_{ij} = sum_k A'_{ik} X_{kj} + sum_k X_{ik} B_{kj}
    // Row 1: (A'*X + X*B)_{11} = a11*x11 + a21*x21 + x11*b11 + x12*b21 = c11
    // => (a11+b11)*x11 + a21*x21 + b21*x12 = c11
    // Row 2: a12*x11 + a22*x21 + x21*b11 + x22*b21 = c21 => a12*x11 + (a22+b11)*x21 + b21*x22 = c21
    // Row 3: a11*x12 + a21*x22 + x11*b12 + x12*b22 = c12 => b12*x11 + (a11+b22)*x12 + a21*x22 = c12
    // Row 4: a12*x12 + a22*x22 + x21*b12 + x22*b22 = c22 => b12*x21 + a12*x12 + (a22+b22)*x22 = c22
    let mut m = [
        [a11 + b11, a21, b21, ZERO],
        [a12, a22 + b11, ZERO, b21],
        [b12, ZERO, a11 + b22, a21],
        [ZERO, b12, a12, a22 + b22],
    ];
    let mut rhs = [c[0], c[1], c[ldc], c[1 + ldc]];
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

/// TRANA: 'N' => op(A)=A, 'T' or 'C' => op(A)=A'. A, C column-major LDA, LDC.
pub fn sb03my(
    trana: char,
    n: usize,
    a: &[f64],
    lda: usize,
    c: &mut [f64],
    ldc: usize,
    scale: &mut f64,
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
    let mut dwork = [0.0f64];
    let anorm = dlanhs_max(n, a, lda, &mut dwork);
    let smin = smlnum.max(eps * anorm);

    let at = |i: usize, j: usize| -> f64 {
        if notrana {
            a[i + j * lda]
        } else {
            a[j + i * lda]
        }
    };

    if notrana {
        // Solve A'*X + X*A = scale*C. Upper-left to bottom-right.
        let mut lnext = 1;
        let mut l = 0;
        while l < n {
            let l1 = l;
            let l2 = if l < n - 1 && a[l + 1 + l * lda].abs() > ZERO { l + 1 } else { l };
            l = l2 + 1;

            let mut knext = l1;
            let mut k = l1;
            while k < n {
                let k1 = k;
                let k2 = if k < n - 1 && a[k + 1 + k * lda].abs() > ZERO { k + 1 } else { k };
                k = k2 + 1;

                if l1 == l2 && k1 == k2 {
                    let mut rhs = c[k1 + l1 * ldc];
                    for i in 0..k1 {
                        rhs -= at(i, k1) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(i, l1);
                    }
                    let a11 = at(k1, k1) + at(l1, l1);
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
                        *scale *= scaloc;
                    }
                    let x11 = rhs / a11;
                    c[k1 + l1 * ldc] = x11;
                    if k1 != l1 {
                        c[l1 + k1 * ldc] = x11;
                    }
                } else if l1 == l2 && k1 != k2 {
                    let mut vec = [ZERO; 2];
                    vec[0] = c[k1 + l1 * ldc];
                    vec[1] = c[k2 + l1 * ldc];
                    for i in 0..k1 {
                        vec[0] -= at(i, k1) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(i, l1);
                        vec[1] -= at(i, k2) * c[i + l1 * ldc] + c[k2 + i * ldc] * at(i, l1);
                    }
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_trans(smin, &a[k1 + k1 * lda..], lda, -at(l1, l1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k2 + l1 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l1 + k2 * ldc] = x[1];
                } else if l1 != l2 && k1 == k2 {
                    let mut vec = [ZERO; 2];
                    vec[0] = c[k1 + l1 * ldc];
                    vec[1] = c[k1 + l2 * ldc];
                    for i in 0..k1 {
                        vec[0] -= at(i, k1) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(i, l1);
                        vec[1] -= at(i, k1) * c[i + l2 * ldc] + c[k1 + i * ldc] * at(i, l2);
                    }
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_trans(smin, &a[l1 + l1 * lda..], lda, -at(k1, k1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k1 + l2 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l2 + k1 * ldc] = x[1];
                } else {
                    let mut vec = [[ZERO; 2]; 2];
                    vec[0][0] = c[k1 + l1 * ldc];
                    vec[0][1] = c[k1 + l2 * ldc];
                    vec[1][0] = c[k2 + l1 * ldc];
                    vec[1][1] = c[k2 + l2 * ldc];
                    for i in 0..k1 {
                        vec[0][0] -= at(i, k1) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(i, l1);
                        vec[0][1] -= at(i, k1) * c[i + l2 * ldc] + c[k1 + i * ldc] * at(i, l2);
                        vec[1][0] -= at(i, k2) * c[i + l1 * ldc] + c[k2 + i * ldc] * at(i, l1);
                        vec[1][1] -= at(i, k2) * c[i + l2 * ldc] + c[k2 + i * ldc] * at(i, l2);
                    }
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut xnorm = ZERO;
                    let mut ierr = 0;
                    if k1 == l1 {
                        let mut t = [a[k1 + k1 * lda], a[k1 + 1 + k1 * lda], a[k1 + (k1 + 1) * lda], a[k1 + 1 + (k1 + 1) * lda]];
                        let mut b = [vec[0][0], vec[1][0], vec[0][1], vec[1][1]];
                        sb03mw::sb03mw(false, true, &t, 2, &b, 2, &mut scaloc, &mut x, 2, &mut xnorm, &mut ierr);
                        x[1] = x[2]; // upper: X(2,1)=X(1,2)
                    } else {
                        solve_sylvester_2x2(
                            &a[k1 + k1 * lda..], lda,
                            &a[l1 + l1 * lda..], lda,
                            &[vec[0][0], vec[1][0], vec[0][1], vec[1][1]], 2,
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
        // Solve A*X + X*A' = scale*C. Bottom-right to upper-left.
        let mut lnext = n;
        let mut l = n;
        while l > 0 {
            l -= 1;
            let l2 = l;
            let l1 = if l > 0 && a[l + (l - 1) * lda].abs() > ZERO { l - 1 } else { l };

            let mut k = l2;
            while k >= l1 && k < n {
                let k2 = k;
                let k1 = if k > 0 && a[k + (k - 1) * lda].abs() > ZERO { k - 1 } else { k };
                if k1 > 0 {
                    k = k1 - 1;
                } else {
                    k = 0;
                    if k1 != 0 {
                        break;
                    }
                }

                if l1 == l2 && k1 == k2 {
                    let mut rhs = c[k1 + l1 * ldc];
                    for i in (k1 + 1)..n {
                        rhs -= at(k1, i) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(l1, i);
                    }
                    let a11 = at(k1, k1) + at(l1, l1);
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
                        *scale *= scaloc;
                    }
                    let x11 = rhs / a11;
                    c[k1 + l1 * ldc] = x11;
                    if k1 != l1 {
                        c[l1 + k1 * ldc] = x11;
                    }
                } else if l1 == l2 && k1 != k2 {
                    let mut vec = [ZERO; 2];
                    vec[0] = c[k1 + l1 * ldc];
                    vec[1] = c[k2 + l1 * ldc];
                    for i in (k2 + 1)..n {
                        vec[0] -= at(k1, i) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(l1, i);
                        vec[1] -= at(k2, i) * c[i + l1 * ldc] + c[k2 + i * ldc] * at(l1, i);
                    }
                    // (A - a_ll*I)'*x = vec => (A' - a_ll*I)*x = vec
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_trans(smin, &a[k1 + k1 * lda..], lda, -at(l1, l1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k2 + l1 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l1 + k2 * ldc] = x[1];
                } else if l1 != l2 && k1 == k2 {
                    let mut vec = [ZERO; 2];
                    vec[0] = c[k1 + l1 * ldc];
                    vec[1] = c[k1 + l2 * ldc];
                    for i in (k1 + 1)..n {
                        vec[0] -= at(k1, i) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(l1, i);
                        vec[1] -= at(k1, i) * c[i + l2 * ldc] + c[k1 + i * ldc] * at(l2, i);
                    }
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut ierr = 0;
                    solve_2x1_trans(smin, &a[l1 + l1 * lda..], lda, -at(k1, k1), &vec, &mut x, &mut scaloc, &mut ierr);
                    if ierr != 0 {
                        *info = 1;
                    }
                    if scaloc != ONE {
                        for j in 0..n {
                            for i in 0..n {
                                c[i + j * ldc] *= scaloc;
                            }
                        }
                        *scale *= scaloc;
                    }
                    c[k1 + l1 * ldc] = x[0];
                    c[k1 + l2 * ldc] = x[1];
                    c[l1 + k1 * ldc] = x[0];
                    c[l2 + k1 * ldc] = x[1];
                } else {
                    let mut vec = [[ZERO; 2]; 2];
                    vec[0][0] = c[k1 + l1 * ldc];
                    vec[0][1] = c[k1 + l2 * ldc];
                    vec[1][0] = c[k2 + l1 * ldc];
                    vec[1][1] = c[k2 + l2 * ldc];
                    for i in (k2 + 1)..n {
                        vec[0][0] -= at(k1, i) * c[i + l1 * ldc] + c[k1 + i * ldc] * at(l1, i);
                        vec[0][1] -= at(k1, i) * c[i + l2 * ldc] + c[k1 + i * ldc] * at(l2, i);
                        vec[1][0] -= at(k2, i) * c[i + l1 * ldc] + c[k2 + i * ldc] * at(l1, i);
                        vec[1][1] -= at(k2, i) * c[i + l2 * ldc] + c[k2 + i * ldc] * at(l2, i);
                    }
                    let mut x = [ZERO; 4];
                    let mut scaloc = ONE;
                    let mut xnorm = ZERO;
                    let mut ierr = 0;
                    if k1 == l1 {
                        let mut t = [a[k1 + k1 * lda], a[k1 + 1 + k1 * lda], a[k1 + (k1 + 1) * lda], a[k1 + 1 + (k1 + 1) * lda]];
                        let mut b = [vec[0][0], vec[1][0], vec[0][1], vec[1][1]];
                        sb03mw::sb03mw(true, true, &t, 2, &b, 2, &mut scaloc, &mut x, 2, &mut xnorm, &mut ierr);
                        x[1] = x[2];
                    } else {
                        solve_sylvester_2x2(
                            &a[l1 + l1 * lda..], lda,
                            &a[k1 + k1 * lda..], lda,
                            &[vec[0][0], vec[1][0], vec[0][1], vec[1][1]], 2,
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
    }
}

fn dlanhs_max(n: usize, a: &[f64], lda: usize, _dwork: &mut [f64]) -> f64 {
    let mut anorm = 0.0;
    for i in 0..n {
        for j in 0..=(i.min(n - 1)) {
            anorm += a[i + j * lda].abs();
        }
    }
    anorm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03my_1x1() {
        let n = 1;
        let mut a = [2.0];
        let mut c = [1.0];
        let mut scale = 0.0;
        let mut info = 0;
        sb03my('N', n, &a, 1, &mut c, 1, &mut scale, &mut info);
        assert_eq!(info, 0);
        assert!((c[0] - 0.25).abs() < 1e-10); // 2*x + 2*x = 1 => x = 0.25
    }
}
