//! MB02CD — Cholesky factorization of positive definite block Toeplitz matrix (SLICOT).
//!
//! Builds full symmetric block Toeplitz T from first block row, then R'*R = T (TYPET='R').

use nalgebra::DMatrix;

/// Job: 'O' = Cholesky R only; 'R' = R and generator G; 'L' = L (inv); 'A' = both; 'G' = G only.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb02CdJob {
    OnlyR,
    RAndG,
    L,
    All,
    GOnly,
}

/// First block row ('R') or first block column ('C').
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum TypeT {
    Row,
    Col,
}

/// Forms symmetric block Toeplitz from first block row t (K×N*K), then Cholesky.
/// On success R contains upper Cholesky factor (N*K × N*K). Returns 0, or 1 if not positive definite.
pub fn mb02cd(
    _job: Mb02CdJob,
    typet: TypeT,
    k: usize,
    n: usize,
    t: &[f64],
    ldt: usize,
    r: &mut [f64],
    ldr: usize,
    _cs: &mut [f64],
    _dwork: &mut [f64],
) -> i32 {
    if k == 0 || n == 0 {
        return 0;
    }
    let nk = n * k;
    if typet == TypeT::Row && (ldt < k || t.len() < k * nk) {
        return -1;
    }
    let mut t_full = DMatrix::<f64>::zeros(nk, nk);
    match typet {
        TypeT::Row => {
            for i in 0..n {
                for j in 0..n {
                    let blk_idx = if j >= i { j - i } else { i - j };
                    for pi in 0..k {
                        for pj in 0..k {
                            let val = if j >= i {
                                t[pi + (blk_idx * k + pj) * ldt]
                            } else {
                                t[pj + (blk_idx * k + pi) * ldt]
                            };
                            t_full[(i * k + pi, j * k + pj)] = val;
                        }
                    }
                }
            }
        }
        TypeT::Col => {
            for i in 0..n {
                for j in 0..n {
                    let blk_idx = if j >= i { j - i } else { i - j };
                    for pi in 0..k {
                        for pj in 0..k {
                            let val = if j >= i {
                                t[(blk_idx * k + pi) + pj * ldt]
                            } else {
                                t[(blk_idx * k + pj) + pi * ldt]
                            };
                            t_full[(i * k + pi, j * k + pj)] = val;
                        }
                    }
                }
            }
        }
    }
    let ch = match t_full.cholesky() {
        Some(c) => c,
        None => return 1,
    };
    let lo = ch.l();
    for i in 0..nk {
        for j in 0..nk {
            r[i + j * ldr] = if i <= j { lo[(j, i)] } else { 0.0 };
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb02cd_trivial() {
        let t = vec![];
        let mut r = vec![0.0];
        let mut cs = vec![0.0];
        let mut dwork = vec![0.0];
        assert_eq!(
            mb02cd(
                Mb02CdJob::OnlyR,
                TypeT::Row,
                0,
                0,
                &t,
                0,
                &mut r,
                0,
                &mut cs,
                &mut dwork,
            ),
            0
        );
    }

    #[test]
    fn test_mb02cd_1x1_block() {
        let k = 1;
        let n = 2;
        let t = vec![2.0, 0.5, 0.5, 1.0];
        let mut r = vec![0.0; 4];
        let mut cs = vec![0.0; 1];
        let mut dwork = vec![0.0; 2];
        assert_eq!(
            mb02cd(
                Mb02CdJob::OnlyR,
                TypeT::Row,
                k,
                n,
                &t,
                k,
                &mut r,
                2,
                &mut cs,
                &mut dwork,
            ),
            0
        );
        assert!(r[0] > 0.0);
    }
}
