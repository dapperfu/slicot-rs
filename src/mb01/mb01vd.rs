//! MB01VD — C = alpha*kron(op(A), op(B)) + beta*C (SLICOT MB01VD)

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01VdTrans {
    NoTrans,
    Trans,
}

/// C = alpha*kron(op(A), op(B)) + beta*C. op(A) is MA×NA, op(B) is MB×NB, C is (MA*MB)×(NA*NB).
pub fn mb01vd(
    trana: Mb01VdTrans,
    tranb: Mb01VdTrans,
    ma: usize,
    na: usize,
    mb: usize,
    nb: usize,
    alpha: f64,
    beta: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
) -> i32 {
    let (ar, ac) = if trana == Mb01VdTrans::NoTrans {
        (ma, na)
    } else {
        (na, ma)
    };
    let (br, bc) = if tranb == Mb01VdTrans::NoTrans {
        (mb, nb)
    } else {
        (nb, mb)
    };
    let mc = ar * br;
    let nc = ac * bc;
    if ldc < mc.max(1) {
        return -14;
    }
    if mc == 0 || nc == 0 {
        return 0;
    }
    for jc in 0..nc {
        let ja = jc / bc;
        let jb = jc % bc;
        for ic in 0..mc {
            let ia = ic / br;
            let ib = ic % br;
            let a_val = if trana == Mb01VdTrans::NoTrans {
                a[ia + ja * lda]
            } else {
                a[ja + ia * lda]
            };
            let b_val = if tranb == Mb01VdTrans::NoTrans {
                b[ib + jb * ldb]
            } else {
                b[jb + ib * ldb]
            };
            c[ic + jc * ldc] = alpha * a_val * b_val + beta * c[ic + jc * ldc];
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb01vd_1x1() {
        let mut c = vec![0.0, 0.0, 0.0, 0.0];
        assert_eq!(
            mb01vd(
                Mb01VdTrans::NoTrans,
                Mb01VdTrans::NoTrans,
                1,
                1,
                2,
                2,
                1.0,
                0.0,
                &[3.0],
                1,
                &[1.0, 0.0, 0.0, 1.0],
                2,
                &mut c,
                2,
            ),
            0
        );
        assert!((c[0] - 3.0).abs() < 1e-14);
        assert!((c[3] - 3.0).abs() < 1e-14);
    }
}
