//! MA02EZ — Store by (skew-)symmetry/Hermitian the other triangle (complex) (SLICOT MA02EZ)
//
// Given upper or lower triangle, fill the other. TRANS: transpose (T) or conjugate transpose (C).
// SKEW: General (G), symmetric/Hermitian (N), skew-symmetric/Hermitian (S).

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02EzUplo {
    Upper,
    Lower,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02EzTrans {
    Trans,
    ConjTrans,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02EzSkew {
    General,
    Symmetric,
    Skew,
}

fn conj_re(ar: f64, _ai: f64) -> f64 {
    ar
}
fn conj_im(_ar: f64, ai: f64) -> f64 {
    -ai
}

/// Fills the other triangle. a_re, a_im column-major LDA×N.
pub fn ma02ez(
    uplo: Ma02EzUplo,
    trans: Ma02EzTrans,
    skew: Ma02EzSkew,
    n: usize,
    a_re: &mut [f64],
    a_im: &mut [f64],
    lda: usize,
) -> i32 {
    if n == 0 || lda < n {
        return 0;
    }
    let idx = |i: usize, j: usize| i + j * lda;

    match uplo {
        Ma02EzUplo::Lower => {
            for i in 0..n {
                for j in (i + 1)..n {
                    let (ri, rj) = (idx(i, j), idx(j, i));
                    match trans {
                        Ma02EzTrans::Trans => {
                            if skew == Ma02EzSkew::Skew {
                                a_re[ri] = -a_re[rj];
                                a_im[ri] = -a_im[rj];
                            } else {
                                a_re[ri] = a_re[rj];
                                a_im[ri] = a_im[rj];
                            }
                        }
                        Ma02EzTrans::ConjTrans => {
                            if skew == Ma02EzSkew::General {
                                a_re[ri] = conj_re(a_re[rj], a_im[rj]);
                                a_im[ri] = conj_im(a_re[rj], a_im[rj]);
                            } else if skew == Ma02EzSkew::Symmetric {
                                a_re[ri] = conj_re(a_re[rj], a_im[rj]);
                                a_im[ri] = conj_im(a_re[rj], a_im[rj]);
                            } else {
                                a_re[ri] = -conj_re(a_re[rj], a_im[rj]);
                                a_im[ri] = -conj_im(a_re[rj], a_im[rj]);
                            }
                        }
                    }
                }
            }
            if trans == Ma02EzTrans::ConjTrans {
                for i in 0..n {
                    let d = idx(i, i);
                    if skew == Ma02EzSkew::Symmetric {
                        a_im[d] = 0.0;
                    } else if skew == Ma02EzSkew::Skew {
                        a_re[d] = 0.0;
                    }
                }
            }
        }
        Ma02EzUplo::Upper => {
            for i in 0..n {
                for j in (i + 1)..n {
                    let (ri, rj) = (idx(j, i), idx(i, j));
                    match trans {
                        Ma02EzTrans::Trans => {
                            if skew == Ma02EzSkew::Skew {
                                a_re[ri] = -a_re[rj];
                                a_im[ri] = -a_im[rj];
                            } else {
                                a_re[ri] = a_re[rj];
                                a_im[ri] = a_im[rj];
                            }
                        }
                        Ma02EzTrans::ConjTrans => {
                            if skew == Ma02EzSkew::General {
                                a_re[ri] = conj_re(a_re[rj], a_im[rj]);
                                a_im[ri] = conj_im(a_re[rj], a_im[rj]);
                            } else if skew == Ma02EzSkew::Symmetric {
                                a_re[ri] = conj_re(a_re[rj], a_im[rj]);
                                a_im[ri] = conj_im(a_re[rj], a_im[rj]);
                            } else {
                                a_re[ri] = -conj_re(a_re[rj], a_im[rj]);
                                a_im[ri] = -conj_im(a_re[rj], a_im[rj]);
                            }
                        }
                    }
                }
            }
            if trans == Ma02EzTrans::ConjTrans {
                for i in 0..n {
                    let d = idx(i, i);
                    if skew == Ma02EzSkew::Symmetric {
                        a_im[d] = 0.0;
                    } else if skew == Ma02EzSkew::Skew {
                        a_re[d] = 0.0;
                    }
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02ez_upper_trans_sym() {
        // Column-major: (0,0)=1, (1,0)=0, (0,1)=2, (1,1)=3 => upper triangle 1,2,3
        let mut a_re = [1.0, 0.0, 2.0, 3.0];
        let mut a_im = [0.0; 4];
        ma02ez(
            Ma02EzUplo::Upper,
            Ma02EzTrans::Trans,
            Ma02EzSkew::Symmetric,
            2,
            &mut a_re,
            &mut a_im,
            2,
        );
        assert_eq!(a_re[1], 2.0);
        assert_eq!(a_re[2], 2.0);
    }
}
