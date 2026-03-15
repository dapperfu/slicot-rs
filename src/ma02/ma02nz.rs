//! MA02NZ — Permute two rows and columns K and L of (skew-)symmetric/Hermitian complex matrix (SLICOT MA02NZ)
//
// Swap row K with row L and column K with column L. Only the referenced triangle is updated.
// Indices 0-based; require 0 <= k < l < n.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02NzUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02NzTrans {
    Trans,
    ConjTrans,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02NzSkew {
    Symmetric,
    Skew,
}

fn conj_re(ar: f64, _ai: f64) -> f64 {
    ar
}
fn conj_im(_ar: f64, ai: f64) -> f64 {
    -ai
}

/// Permutes rows and columns k and l. a_re, a_im column-major LDA×N. k and l 0-based, k < l.
pub fn ma02nz(
    uplo: Ma02NzUplo,
    trans: Ma02NzTrans,
    skew: Ma02NzSkew,
    n: usize,
    k: usize,
    l: usize,
    a_re: &mut [f64],
    a_im: &mut [f64],
    lda: usize,
) -> i32 {
    if n == 0 || k >= l || l >= n {
        return 0;
    }
    let idx = |i: usize, j: usize| i + j * lda;

    let (tk_re, tk_im) = (a_re[idx(k, k)], a_im[idx(k, k)]);
    a_re[idx(k, k)] = a_re[idx(l, l)];
    a_im[idx(k, k)] = a_im[idx(l, l)];
    a_re[idx(l, l)] = tk_re;
    a_im[idx(l, l)] = tk_im;

    match uplo {
        Ma02NzUplo::Lower => {
            for j in 0..k {
                a_re.swap(idx(k, j), idx(l, j));
                a_im.swap(idx(k, j), idx(l, j));
            }
            if l > k + 1 {
                match trans {
                    Ma02NzTrans::Trans => {
                        if skew == Ma02NzSkew::Symmetric {
                            for i in (k + 1)..l {
                                a_re.swap(idx(i, k), idx(l, i));
                                a_im.swap(idx(i, k), idx(l, i));
                            }
                        } else {
                            a_re[idx(l, k)] = -a_re[idx(l, k)];
                            a_im[idx(l, k)] = -a_im[idx(l, k)];
                            for i in (k + 1)..l {
                                let (tr, ti) = (a_re[idx(l, i)], a_im[idx(l, i)]);
                                a_re[idx(l, i)] = -a_re[idx(i, k)];
                                a_im[idx(l, i)] = -a_im[idx(i, k)];
                                a_re[idx(i, k)] = -tr;
                                a_im[idx(i, k)] = -ti;
                            }
                        }
                    }
                    Ma02NzTrans::ConjTrans => {
                        if skew == Ma02NzSkew::Symmetric {
                            let (lr, li) = (a_re[idx(l, k)], a_im[idx(l, k)]);
                            a_re[idx(l, k)] = conj_re(lr, li);
                            a_im[idx(l, k)] = conj_im(lr, li);
                            for i in (k + 1)..l {
                                let (lr, li) = (a_re[idx(l, i)], a_im[idx(l, i)]);
                                let (ikr, iki) = (a_re[idx(i, k)], a_im[idx(i, k)]);
                                a_re[idx(l, i)] = conj_re(ikr, iki);
                                a_im[idx(l, i)] = conj_im(ikr, iki);
                                a_re[idx(i, k)] = conj_re(lr, li);
                                a_im[idx(i, k)] = conj_im(lr, li);
                            }
                        } else {
                            a_re[idx(l, k)] = -conj_re(a_re[idx(l, k)], a_im[idx(l, k)]);
                            a_im[idx(l, k)] = conj_im(a_re[idx(l, k)], a_im[idx(l, k)]);
                            for i in (k + 1)..l {
                                let (lr, li) = (a_re[idx(l, i)], a_im[idx(l, i)]);
                                let (ikr, iki) = (a_re[idx(i, k)], a_im[idx(i, k)]);
                                a_re[idx(l, i)] = -conj_re(ikr, iki);
                                a_im[idx(l, i)] = conj_im(ikr, iki);
                                a_re[idx(i, k)] = -conj_re(lr, li);
                                a_im[idx(i, k)] = conj_im(lr, li);
                            }
                        }
                    }
                }
            }
            for i in (l + 1)..n {
                a_re.swap(idx(i, k), idx(i, l));
                a_im.swap(idx(i, k), idx(i, l));
            }
        }
        Ma02NzUplo::Upper => {
            for j in 0..k {
                a_re.swap(idx(j, k), idx(j, l));
                a_im.swap(idx(j, k), idx(j, l));
            }
            if l > k + 1 {
                match trans {
                    Ma02NzTrans::Trans => {
                        if skew == Ma02NzSkew::Symmetric {
                            for i in (k + 1)..l {
                                a_re.swap(idx(k, i), idx(i, l));
                                a_im.swap(idx(k, i), idx(i, l));
                            }
                        } else {
                            a_re[idx(k, l)] = -a_re[idx(k, l)];
                            a_im[idx(k, l)] = -a_im[idx(k, l)];
                            for i in (k + 1)..l {
                                let (tr, ti) = (a_re[idx(i, l)], a_im[idx(i, l)]);
                                a_re[idx(i, l)] = -a_re[idx(k, i)];
                                a_im[idx(i, l)] = -a_im[idx(k, i)];
                                a_re[idx(k, i)] = -tr;
                                a_im[idx(k, i)] = -ti;
                            }
                        }
                    }
                    Ma02NzTrans::ConjTrans => {
                        if skew == Ma02NzSkew::Symmetric {
                            let (kr, ki) = (a_re[idx(k, l)], a_im[idx(k, l)]);
                            a_re[idx(k, l)] = conj_re(kr, ki);
                            a_im[idx(k, l)] = conj_im(kr, ki);
                            for i in (k + 1)..l {
                                let (ilr, ili) = (a_re[idx(i, l)], a_im[idx(i, l)]);
                                let (kir, kii) = (a_re[idx(k, i)], a_im[idx(k, i)]);
                                a_re[idx(i, l)] = conj_re(kir, kii);
                                a_im[idx(i, l)] = conj_im(kir, kii);
                                a_re[idx(k, i)] = conj_re(ilr, ili);
                                a_im[idx(k, i)] = conj_im(ilr, ili);
                            }
                        } else {
                            a_re[idx(k, l)] = -conj_re(a_re[idx(k, l)], a_im[idx(k, l)]);
                            a_im[idx(k, l)] = conj_im(a_re[idx(k, l)], a_im[idx(k, l)]);
                            for i in (k + 1)..l {
                                let (ilr, ili) = (a_re[idx(i, l)], a_im[idx(i, l)]);
                                let (kir, kii) = (a_re[idx(k, i)], a_im[idx(k, i)]);
                                a_re[idx(i, l)] = -conj_re(kir, kii);
                                a_im[idx(i, l)] = conj_im(kir, kii);
                                a_re[idx(k, i)] = -conj_re(ilr, ili);
                                a_im[idx(k, i)] = conj_im(ilr, ili);
                            }
                        }
                    }
                }
            }
            for j in (l + 1)..n {
                a_re.swap(idx(k, j), idx(l, j));
                a_im.swap(idx(k, j), idx(l, j));
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02nz_upper_trans_sym() {
        // 3×3 upper triangle column-major: (0,0)=1,(0,1)=2,(0,2)=3,(1,1)=4,(1,2)=5,(2,2)=6. Permute row/col 0 and 1.
        let mut a_re = [1.0, 2.0, 3.0, 0.0, 4.0, 5.0, 0.0, 0.0, 6.0];
        let mut a_im = [0.0; 9];
        ma02nz(
            Ma02NzUplo::Upper,
            Ma02NzTrans::Trans,
            Ma02NzSkew::Symmetric,
            3,
            0,
            1,
            &mut a_re,
            &mut a_im,
            3,
        );
        assert_eq!(a_re[0], 4.0);
        assert_eq!(a_re[4], 1.0);
    }
}
