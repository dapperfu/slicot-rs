//! MB04ND — RQ factorization of first block row and apply to second block row.
//!
//! [A R] * Q' = [0 R_bar], [C B] updated to [C_bar B_bar]. Uses DLARFG and MB04NY.

use crate::mb04::blas::dlarfg;
use crate::mb04::mb04ny::mb04ny;

/// UPLO: 'U' = A upper trapezoidal/triangular, 'F' = A full.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb04ndUplo {
    Upper,
    Full,
}

/// RQ factorization and apply. Column-major arrays.
/// R is N×N upper triangular (LDR), A is N×P (LDA), B is M×N (LDB), C is M×P (LDC).
/// TAU length N, DWORK length max(N-1, M).
#[allow(clippy::too_many_arguments)]
pub fn mb04nd(
    uplo: Mb04ndUplo,
    n: usize,
    m: usize,
    p: usize,
    r: &mut [f64],
    ldr: usize,
    a: &mut [f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
    c: &mut [f64],
    ldc: usize,
    tau: &mut [f64],
    dwork: &mut [f64],
) {
    if n == 0 || p == 0 {
        return;
    }
    let luplo = uplo == Mb04ndUplo::Upper;
    if luplo {
        for i in (0..n).rev() {
            let im = (n - i).min(p);
            let ip_1based = (p + i).saturating_sub(n).max(1);
            let ip = ip_1based.saturating_sub(1);
            let rii_idx = i + i * ldr;
            let v_start = i + ip * lda;
            dlarfg(im + 1, &mut r[rii_idx], &mut a[v_start..], lda, &mut tau[i]);
            if i > 0 {
                let a_col_len = i;
                let r_col_start = i * ldr;
                let v_copy: Vec<f64> = (0..im).map(|k| a[v_start + k * lda]).collect();
                mb04ny(
                    i,
                    im,
                    &v_copy,
                    1,
                    tau[i],
                    &mut r[r_col_start..r_col_start + a_col_len],
                    ldr,
                    &mut a[ip * lda..],
                    lda,
                    dwork,
                );
            }
            if m > 0 {
                let v_copy: Vec<f64> = (0..im).map(|k| a[v_start + k * lda]).collect();
                mb04ny(
                    m,
                    im,
                    &v_copy,
                    1,
                    tau[i],
                    &mut b[i * ldb..],
                    ldb,
                    &mut c[ip * ldc..],
                    ldc,
                    dwork,
                );
            }
        }
    } else {
        for i in (1..n).rev() {
            let rii_idx = i + i * ldr;
            let v_start = i * lda;
            dlarfg(
                p + 1,
                &mut r[rii_idx],
                &mut a[v_start..],
                lda,
                &mut tau[i],
            );
            let v_copy: Vec<f64> = (0..p).map(|k| a[v_start + k * lda]).collect();
            mb04ny(
                i,
                p,
                &v_copy,
                1,
                tau[i],
                &mut r[i * ldr..i * ldr + i],
                ldr,
                a,
                lda,
                dwork,
            );
        }
        dlarfg(p + 1, &mut r[0], &mut a[0..], lda, &mut tau[0]);
        if m > 0 {
            for i in (0..n).rev() {
                let v_start = i * lda;
                let v_copy: Vec<f64> = (0..p).map(|k| a[v_start + k * lda]).collect();
                mb04ny(
                    m,
                    p,
                    &v_copy,
                    1,
                    tau[i],
                    &mut b[i * ldb..],
                    ldb,
                    c,
                    ldc,
                    dwork,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb04nd_early_return() {
        let mut r = vec![0.0_f64];
        let mut a = vec![0.0_f64];
        let mut b = vec![0.0_f64];
        let mut c = vec![0.0_f64];
        let mut tau = vec![0.0_f64];
        let mut dwork = vec![0.0_f64; 2];
        mb04nd(
            Mb04ndUplo::Full,
            0,
            0,
            1,
            &mut r,
            1,
            &mut a,
            1,
            &mut b,
            1,
            &mut c,
            1,
            &mut tau,
            &mut dwork,
        );
    }
}
