//! MB01UZ — T := alpha*op(T)*A or T := alpha*A*op(T) (complex) (SLICOT MB01UZ)
//
// op(T) = T, T', or conj(T'). Result overwrites the leading M×N part of T.
// Storage: t_re, t_im and a_re, a_im column-major; alpha_re, alpha_im; zwork length >= M (SIDE='L') or N (SIDE='R').

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UzSide {
    Left,
    Right,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UzUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01UzTrans {
    NoTrans,
    Trans,
    ConjTrans,
}

fn cmul(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    (ar * br - ai * bi, ar * bi + ai * br)
}
/// Overwrites the leading M×N part of T with the product. zwork length >= k (k = M if Left, N if Right).
pub fn mb01uz(
    side: Mb01UzSide,
    uplo: Mb01UzUplo,
    trans: Mb01UzTrans,
    m: usize,
    n: usize,
    alpha_re: f64,
    alpha_im: f64,
    t_re: &mut [f64],
    t_im: &mut [f64],
    ldt: usize,
    a_re: &[f64],
    a_im: &[f64],
    lda: usize,
    zwork_re: &mut [f64],
    zwork_im: &mut [f64],
) -> i32 {
    if !matches!(side, Mb01UzSide::Left | Mb01UzSide::Right) {
        return -1;
    }
    if !matches!(uplo, Mb01UzUplo::Upper | Mb01UzUplo::Lower) {
        return -2;
    }
    if !matches!(trans, Mb01UzTrans::NoTrans | Mb01UzTrans::Trans | Mb01UzTrans::ConjTrans) {
        return -3;
    }
    let k = if side == Mb01UzSide::Left { m } else { n };
    if ldt < k.max(m).max(n) || lda < m.max(1) {
        return -8;
    }
    let wrk_min = if (alpha_re != 0.0 || alpha_im != 0.0) && m > 0 && n > 0 {
        k
    } else {
        1
    };
    if zwork_re.len() < wrk_min || zwork_im.len() < wrk_min {
        return -14;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    if alpha_re == 0.0 && alpha_im == 0.0 {
        for j in 0..n {
            for i in 0..m {
                let idx = i + j * ldt;
                t_re[idx] = 0.0;
                t_im[idx] = 0.0;
            }
        }
        return 0;
    }

    let idx_t = |i: usize, j: usize| i + j * ldt;
    let idx_a = |i: usize, j: usize| i + j * lda;

    let get_t = |tre: &[f64], tim: &[f64], i: usize, j: usize| -> (f64, f64) {
        if (uplo == Mb01UzUplo::Upper && i <= j) || (uplo == Mb01UzUplo::Lower && i >= j) {
            (tre[idx_t(i, j)], tim[idx_t(i, j)])
        } else {
            (0.0, 0.0)
        }
    };

    if side == Mb01UzSide::Left {
        for i in 0..m {
            for p in 0..m {
                let (tr, ti) = match trans {
                    Mb01UzTrans::NoTrans => get_t(t_re, t_im, i, p),
                    Mb01UzTrans::Trans => get_t(t_re, t_im, p, i),
                    Mb01UzTrans::ConjTrans => {
                        let (r, im) = get_t(t_re, t_im, p, i);
                        (r, -im)
                    }
                };
                zwork_re[p] = tr;
                zwork_im[p] = ti;
            }
            for j in 0..n {
                let mut sum_re = 0.0;
                let mut sum_im = 0.0;
                for p in 0..m {
                    let (pr, pi) = cmul(zwork_re[p], zwork_im[p], a_re[idx_a(p, j)], a_im[idx_a(p, j)]);
                    sum_re += pr;
                    sum_im += pi;
                }
                let (sr, si) = cmul(alpha_re, alpha_im, sum_re, sum_im);
                t_re[idx_t(i, j)] = sr;
                t_im[idx_t(i, j)] = si;
            }
        }
    } else {
        for j in 0..n {
            for p in 0..n {
                let (tr, ti) = match trans {
                    Mb01UzTrans::NoTrans => get_t(t_re, t_im, p, j),
                    Mb01UzTrans::Trans => get_t(t_re, t_im, j, p),
                    Mb01UzTrans::ConjTrans => {
                        let (r, im) = get_t(t_re, t_im, j, p);
                        (r, -im)
                    }
                };
                zwork_re[p] = tr;
                zwork_im[p] = ti;
            }
            for i in 0..m {
                let mut sum_re = 0.0;
                let mut sum_im = 0.0;
                for p in 0..n {
                    let (pr, pi) = cmul(a_re[idx_a(i, p)], a_im[idx_a(i, p)], zwork_re[p], zwork_im[p]);
                    sum_re += pr;
                    sum_im += pi;
                }
                let (sr, si) = cmul(alpha_re, alpha_im, sum_re, sum_im);
                t_re[idx_t(i, j)] = sr;
                t_im[idx_t(i, j)] = si;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01uz_left_upper_notrans() {
        let m = 2;
        let n = 2;
        let mut t_re = [1.0, 0.0, 0.0, 1.0];
        let mut t_im = [0.0; 4];
        let a_re = [1.0, 0.0, 0.0, 1.0];
        let a_im = [0.0; 4];
        let mut z_re = vec![0.0; 2];
        let mut z_im = vec![0.0; 2];
        assert_eq!(
            mb01uz(
                Mb01UzSide::Left,
                Mb01UzUplo::Upper,
                Mb01UzTrans::NoTrans,
                m,
                n,
                1.0,
                0.0,
                &mut t_re,
                &mut t_im,
                2,
                &a_re,
                &a_im,
                2,
                &mut z_re,
                &mut z_im,
            ),
            0
        );
        assert!((t_re[0] - 1.0).abs() < 1e-14);
        assert!((t_re[3] - 1.0).abs() < 1e-14);
    }
}
