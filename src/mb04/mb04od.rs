//! MB04OD — QR factorization of first block column and apply to second block column.
//!
//! Q'*[R B; A C] = [R_bar B_bar; 0 C_bar]. Uses DLARFG and MB04OY.

use crate::mb04::blas::dlarfg;
use crate::mb04::mb04oy::mb04oy;

/// UPLO: 'U' = A upper trapezoidal/triangular, 'F' = A full.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb04odUplo {
    Upper,
    Full,
}

/// QR factorization and apply. Column-major arrays.
/// R is N×N (LDR), A is P×N (LDA), B is N×M (LDB), C is P×M (LDC).
/// TAU length N, DWORK length max(N-1, M).
#[allow(clippy::too_many_arguments)]
pub fn mb04od(
    uplo: Mb04odUplo,
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
    let luplo = uplo == Mb04odUplo::Upper;
    if luplo {
        for i in 0..n {
            let im = i.min(p);
            let rii_idx = i + i * ldr;
            let v_start = i * lda;
            dlarfg(im + 1, &mut r[rii_idx], &mut a[v_start..], 1, &mut tau[i]);
            if n - i > 1 {
                let v_copy: Vec<f64> = (0..im).map(|k| a[v_start + k]).collect();
                mb04oy(
                    im,
                    n - i - 1,
                    &v_copy,
                    tau[i],
                    &mut r[i * ldr + i + 1..],
                    ldr,
                    &mut a[(i + 1) * lda..],
                    lda,
                    dwork,
                );
            }
            if m > 0 {
                let v_copy: Vec<f64> = (0..im).map(|k| a[v_start + k]).collect();
                mb04oy(
                    im,
                    m,
                    &v_copy,
                    tau[i],
                    &mut b[i * ldb..],
                    ldb,
                    c,
                    ldc,
                    dwork,
                );
            }
        }
    } else {
        for i in 0..(n - 1) {
            let rii_idx = i + i * ldr;
            let v_start = i * lda;
            dlarfg(p + 1, &mut r[rii_idx], &mut a[v_start..], 1, &mut tau[i]);
            let v_copy: Vec<f64> = (0..p).map(|k| a[v_start + k]).collect();
            mb04oy(
                p,
                n - i - 1,
                &v_copy,
                tau[i],
                &mut r[i * ldr + i + 1..],
                ldr,
                &mut a[(i + 1) * lda..],
                lda,
                dwork,
            );
        }
        let rnn_idx = n - 1 + (n - 1) * ldr;
        let v_start = (n - 1) * lda;
        dlarfg(p + 1, &mut r[rnn_idx], &mut a[v_start..], 1, &mut tau[n - 1]);
        if m > 0 {
            for i in 0..n {
                let v_start = i * lda;
                let v_copy: Vec<f64> = (0..p).map(|k| a[v_start + k]).collect();
                mb04oy(
                    p,
                    m,
                    &v_copy,
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
    fn test_mb04od_early_return() {
        let mut r = vec![0.0_f64];
        let mut a = vec![0.0_f64];
        let mut b = vec![0.0_f64];
        let mut c = vec![0.0_f64];
        let mut tau = vec![0.0_f64];
        let mut dwork = vec![0.0_f64; 2];
        mb04od(
            Mb04odUplo::Full,
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
