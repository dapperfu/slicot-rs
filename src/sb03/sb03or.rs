//! SB03OR — Solve Sylvester equation for block upper triangular S and small A (M=1 or 2).
//! Continuous: op(S)'*X + X*op(A) = scale*C. Discrete: op(S)'*X*op(A) - X = scale*C.

use crate::mb04::blas::{ddot, dscal};
use crate::sb03::dlasy2::dlasy2;
use crate::sb04::sb04px::sb04px;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Solves Sylvester equation. DISCR: discrete (use SB04PX), else continuous (use DLASY2).
/// LTRANS: use transpose of S and A in the equation.
/// S is N×N block upper Hessenberg (1×1 and 2×2 blocks), A is M×M (M=1 or 2).
/// C is overwritten with solution X. Returns INFO (0 = ok, 1 = perturbed).
#[allow(clippy::too_many_arguments)]
pub fn sb03or(
    discr: bool,
    ltrans: bool,
    n: usize,
    m: usize,
    s: &[f64],
    lds: usize,
    a: &[f64],
    lda: usize,
    c: &mut [f64],
    ldc: usize,
    scale: &mut f64,
) -> i32 {
    if n == 0 {
        return 0;
    }
    if !(m == 1 || m == 2) {
        return -4; // invalid M
    }
    *scale = ONE;
    let mut infom = 0_i32;
    let isgn = 1;
    let tbyt = m == 2;

    let mut at = [ZERO; 4];
    at[0] = a[0];
    if tbyt {
        at[1] = a[lda];
        at[2] = a[1];
        at[3] = a[lda + 1];
    }

    let mut vec = [ZERO; 4];
    let mut x = [ZERO; 4];

    if ltrans {
        let mut lnext = n;
        for l in (0..n).rev() {
            if l >= lnext {
                continue;
            }
            let mut l1 = l;
            let mut l2 = l;
            if l > 0 {
                if s.get(l + (l - 1) * lds).copied().unwrap_or(ZERO) != ZERO {
                    l1 = l - 1;
                }
                lnext = l1.saturating_sub(1);
            }
            let dl = l2 - l1 + 1;
            let l2p1 = (l2 + 1).min(n);

            if discr {
                let g11 = -ddot(
                    n - l2,
                    &s[l1 + l2p1 * lds..],
                    lds,
                    &c[l2p1..],
                    1,
                );
                if tbyt {
                    let g12 = -ddot(
                        n - l2,
                        &s[l1 + l2p1 * lds..],
                        lds,
                        &c[l2p1 + ldc..],
                        1,
                    );
                    vec[0] = c[l1] + g11 * at[0] + g12 * at[2];
                    vec[1] = c[l1 + ldc] + g11 * at[1] + g12 * at[3];
                } else {
                    vec[0] = c[l1] + g11 * at[0];
                }
                if dl != 1 {
                    let g21 = -ddot(
                        n - l2,
                        &s[l2 + l2p1 * lds..],
                        lds,
                        &c[l2p1..],
                        1,
                    );
                    if tbyt {
                        let g22 = -ddot(
                            n - l2,
                            &s[l2 + l2p1 * lds..],
                            lds,
                            &c[l2p1 + ldc..],
                            1,
                        );
                        vec[2] = c[l2] + g21 * at[0] + g22 * at[2];
                        vec[3] = c[l2 + ldc] + g21 * at[1] + g22 * at[3];
                    } else {
                        vec[2] = c[l2] + g21 * at[0];
                    }
                }
                let mut scaloc = ONE;
                let mut xnorm = ZERO;
                let info = sb04px(
                    false,
                    false,
                    -isgn,
                    dl,
                    m,
                    &s[l1 + l1 * lds..],
                    lds,
                    &at,
                    2,
                    &vec,
                    2,
                    &mut scaloc,
                    &mut x,
                    2,
                    &mut xnorm,
                );
                infom = infom.max(info);
                if scaloc != ONE {
                    for j in 0..m {
                        dscal(n, scaloc, &mut c[j * ldc..], 1);
                    }
                    *scale *= scaloc;
                }
            } else {
                vec[0] = c[l1]
                    - ddot(n - l2, &s[l1 + l2p1 * lds..], lds, &c[l2p1..], 1);
                if tbyt {
                    vec[1] = c[l1 + ldc]
                        - ddot(n - l2, &s[l1 + l2p1 * lds..], lds, &c[l2p1 + ldc..], 1);
                }
                if dl != 1 {
                    vec[2] = c[l2]
                        - ddot(n - l2, &s[l2 + l2p1 * lds..], lds, &c[l2p1..], 1);
                    if tbyt {
                        vec[3] = c[l2 + ldc]
                            - ddot(
                                n - l2,
                                &s[l2 + l2p1 * lds..],
                                lds,
                                &c[l2p1 + ldc..],
                                1,
                            );
                    }
                }
                let mut scaloc = ONE;
                let mut xnorm = ZERO;
                let info = dlasy2(
                    false,
                    false,
                    isgn,
                    dl,
                    m,
                    &s[l1 + l1 * lds..],
                    lds,
                    &at,
                    2,
                    &vec,
                    2,
                    &mut scaloc,
                    &mut x,
                    2,
                    &mut xnorm,
                );
                infom = infom.max(info);
                if scaloc != ONE {
                    for j in 0..m {
                        dscal(n, scaloc, &mut c[j * ldc..], 1);
                    }
                    *scale *= scaloc;
                }
            }
            c[l1] = x[0];
            if tbyt {
                c[l1 + ldc] = x[1];
            }
            if dl != 1 {
                c[l2] = x[2];
                if tbyt {
                    c[l2 + ldc] = x[3];
                }
            }
        }
    } else {
        let mut lnext = 0;
        for l in 0..n {
            if l < lnext {
                continue;
            }
            let mut l1 = l;
            let mut l2 = l;
            if l < n - 1 {
                if s.get(l + 1 + l * lds).copied().unwrap_or(ZERO) != ZERO {
                    l2 = l + 1;
                }
                lnext = l2 + 1;
            }
            let dl = l2 - l1 + 1;

            if discr {
                let g11 = -ddot(l1, &c[..], 1, &s[l1 * lds..], 1);
                if tbyt {
                    let g21 = -ddot(l1, &c[ldc..], 1, &s[l1 * lds..], 1);
                    vec[0] = c[l1] + at[0] * g11 + at[1] * g21;
                    vec[2] = c[l1 + ldc] + at[2] * g11 + at[3] * g21;
                } else {
                    vec[0] = c[l1] + at[0] * g11;
                }
                if dl != 1 {
                    let g12 = -ddot(l1, &c[..], 1, &s[l2 * lds..], 1);
                    if tbyt {
                        let g22 = -ddot(l1, &c[ldc..], 1, &s[l2 * lds..], 1);
                        vec[1] = c[l2] + at[0] * g12 + at[1] * g22;
                        vec[3] = c[l2 + ldc] + at[2] * g12 + at[3] * g22;
                    } else {
                        vec[1] = c[l2] + at[0] * g12;
                    }
                }
                let mut scaloc = ONE;
                let mut xnorm = ZERO;
                let info = sb04px(
                    false,
                    false,
                    -isgn,
                    m,
                    dl,
                    &at,
                    2,
                    &s[l1 + l1 * lds..],
                    lds,
                    &vec,
                    2,
                    &mut scaloc,
                    &mut x,
                    2,
                    &mut xnorm,
                );
                infom = infom.max(info);
                if scaloc != ONE {
                    for j in 0..m {
                        dscal(n, scaloc, &mut c[j * ldc..], 1);
                    }
                    *scale *= scaloc;
                }
            } else {
                vec[0] = c[l1] - ddot(l1, &c[..], 1, &s[l1 * lds..], 1);
                if tbyt {
                    vec[2] = c[l1 + ldc] - ddot(l1, &c[ldc..], 1, &s[l1 * lds..], 1);
                }
                if dl != 1 {
                    vec[1] = c[l2] - ddot(l1, &c[..], 1, &s[l2 * lds..], 1);
                    if tbyt {
                        vec[3] = c[l2 + ldc]
                            - ddot(l1, &c[ldc..], 1, &s[l2 * lds..], 1);
                    }
                }
                let mut scaloc = ONE;
                let mut xnorm = ZERO;
                let info = dlasy2(
                    false,
                    false,
                    isgn,
                    m,
                    dl,
                    &at,
                    2,
                    &s[l1 + l1 * lds..],
                    lds,
                    &vec,
                    2,
                    &mut scaloc,
                    &mut x,
                    2,
                    &mut xnorm,
                );
                infom = infom.max(info);
                if scaloc != ONE {
                    for j in 0..m {
                        dscal(n, scaloc, &mut c[j * ldc..], 1);
                    }
                    *scale *= scaloc;
                }
            }
            c[l1] = x[0];
            if tbyt {
                c[l1 + ldc] = x[2];
            }
            if dl != 1 {
                c[l2] = x[1];
                if tbyt {
                    c[l2 + ldc] = x[3];
                }
            }
        }
    }
    infom
}

/// Compatibility wrapper: (n, a, x). Solves S*X + X*A = C with S = a (n×n), A = leading m×m of a, C = x; x overwritten with solution.
pub fn sb03or_compat(n: usize, a: &nalgebra::DMatrix<f64>, x: &mut nalgebra::DMatrix<f64>) -> i32 {
    if n == 0 || x.is_empty() {
        return 0;
    }
    let m = x.ncols();
    if m != 1 && m != 2 {
        return -4;
    }
    let mut scale = 1.0_f64;
    let s = a.as_slice();
    let lds = a.nrows();
    let ldc = x.nrows();
    let lda = a.nrows();
    let a_small: Vec<f64> = if m == 1 {
        vec![a[(0, 0)]]
    } else {
        vec![a[(0, 0)], a[(1, 0)], a[(0, 1)], a[(1, 1)]]
    };
    let info = sb03or(
        false,
        false,
        n,
        m,
        s,
        lds,
        &a_small,
        m,
        x.as_mut_slice(),
        ldc,
        &mut scale,
    );
    info
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03or_n0() {
        let s = vec![0.0_f64];
        let a = vec![1.0_f64];
        let mut c = vec![0.0_f64];
        assert_eq!(
            sb03or(false, false, 0, 1, &s, 1, &a, 1, &mut c, 1, &mut 1.0),
            0
        );
    }

    #[test]
    fn test_sb03or_1x1() {
        let s = vec![2.0_f64];
        let a = vec![3.0_f64];
        let mut c = vec![1.0_f64];
        let mut scale = 1.0_f64;
        let info = sb03or(false, false, 1, 1, &s, 1, &a, 1, &mut c, 1, &mut scale);
        assert_eq!(info, 0);
        assert!((scale - 1.0).abs() < 1e-10);
        assert!((c[0] - 0.2).abs() < 1e-10);
    }
}
