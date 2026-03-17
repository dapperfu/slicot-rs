//! MB01QD — Scale matrix A by scalar CTO/CFROM (SLICOT MB01QD)
//!
//! Multiplies A by CTO/CFROM with optional overflow protection. Supports TYPE 'G' (full),
//! 'L' (lower triangular), 'U' (upper triangular), and 'H' (upper Hessenberg).
//! Block structure (NBL > 0) supported for L, U, H.

use nalgebra::DMatrix;
use std::cmp::min;

/// Matrix storage type (1:1 with Fortran TYPE).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb01QdType {
    /// Full matrix ('G').
    General,
    /// Lower triangular ('L').
    Lower,
    /// Upper triangular ('U').
    Upper,
    /// Upper Hessenberg ('H').
    Hessenberg,
}

/// One step of overflow-safe scaling (1:1 with Fortran label 10). Returns (MUL, DONE, next CFROMC, next CTOC).
fn scale_mul_step(cfromc: f64, ctoc: f64) -> (f64, bool, f64, f64) {
    let smlnum = f64::MIN_POSITIVE;
    let bignum = 1.0 / smlnum;
    let cfrom1 = cfromc * smlnum;
    let cto1 = ctoc / bignum;
    if cfrom1.abs() > ctoc.abs() && ctoc != 0.0 {
        (smlnum, false, cfrom1, ctoc)
    } else if cto1.abs() > cfromc.abs() {
        (bignum, false, cfromc, cto1)
    } else {
        (ctoc / cfromc, true, cfromc, ctoc)
    }
}

/// Multiplies the M×N matrix A by the scalar CTO/CFROM (1:1 with SLICOT MB01QD).
/// CFROM must be nonzero. Supports General, Lower, Upper, Hessenberg; NBL=0 or block structure.
///
/// # Returns
/// 0 on success; < 0 if the i-th argument is invalid.
pub fn mb01qd(
    typ: Mb01QdType,
    m: usize,
    n: usize,
    _kl: i32,
    _ku: i32,
    cfrom: f64,
    cto: f64,
    nbl: i32,
    nrows: &[i32],
    a: &mut DMatrix<f64>,
) -> i32 {
    if cfrom == 0.0 {
        return -6;
    }
    if m == 0 || n == 0 {
        return 0;
    }
    let mut cfromc = cfrom;
    let mut ctoc = cto;
    loop {
        let (mul, done) = {
            let (mul, done, next_cf, next_ct) = scale_mul_step(cfromc, ctoc);
            cfromc = next_cf;
            ctoc = next_ct;
            (mul, done)
        };
        mb01qd_apply(typ, m, n, nbl, nrows, a, mul);
        if done {
            break;
        }
    }
    0
}

fn mb01qd_apply(
    typ: Mb01QdType,
    m: usize,
    n: usize,
    nbl: i32,
    nrows: &[i32],
    a: &mut DMatrix<f64>,
    mul: f64,
) -> i32 {
    let noblc = nbl == 0;
    match typ {
        Mb01QdType::General => {
            for j in 0..n {
                for i in 0..m {
                    a[(i, j)] *= mul;
                }
            }
        }
        Mb01QdType::Lower => {
            if noblc {
                for j in 0..n {
                    for i in j..m {
                        a[(i, j)] *= mul;
                    }
                }
            } else {
                let mut jfin = 0usize;
                for k in 0..(nbl as usize) {
                    if k >= nrows.len() {
                        break;
                    }
                    let jini = jfin;
                    jfin += nrows[k] as usize;
                    for j in jini..jfin.min(n) {
                        for i in jini..m {
                            a[(i, j)] *= mul;
                        }
                    }
                }
            }
        }
        Mb01QdType::Upper => {
            if noblc {
                for j in 0..n {
                    let imax = min(j + 1, m);
                    for i in 0..imax {
                        a[(i, j)] *= mul;
                    }
                }
            } else {
                let mut jfin = 0usize;
                for k in 0..(nbl as usize) {
                    if k >= nrows.len() {
                        break;
                    }
                    let jini = jfin;
                    jfin += nrows[k] as usize;
                    let jfin_n = if k == (nbl as usize).saturating_sub(1) {
                        n
                    } else {
                        jfin
                    };
                    for j in jini..jfin_n.min(n) {
                        let imax = min(jfin, m);
                        for i in 0..imax {
                            a[(i, j)] *= mul;
                        }
                    }
                }
            }
        }
        Mb01QdType::Hessenberg => {
            if noblc {
                for j in 0..n {
                    let imax = min(j + 2, m); // 1-based I=1..min(J+1,M) -> 0-based 0..min(j+2,m)
                    for i in 0..imax {
                        a[(i, j)] *= mul;
                    }
                }
            } else {
                let mut jfin = 0usize;
                for k in 0..(nbl as usize) {
                    if k >= nrows.len() {
                        break;
                    }
                    let jini = jfin;
                    jfin += nrows[k] as usize;
                    let ifin = if k + 1 < nrows.len() {
                        jfin + nrows[k + 1] as usize
                    } else {
                        n
                    };
                    for j in jini..jfin.min(n) {
                        let imax = min(ifin, m);
                        for i in 0..imax {
                            a[(i, j)] *= mul;
                        }
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
    fn test_mb01qd_general() {
        let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let nrows: [i32; 0] = [];
        assert_eq!(
            mb01qd(
                Mb01QdType::General,
                2,
                2,
                0,
                0,
                1.0,
                2.0,
                0,
                &nrows,
                &mut a
            ),
            0
        );
        assert_eq!(a[(0, 0)], 2.0);
        assert_eq!(a[(1, 1)], 8.0);
    }

    #[test]
    fn test_mb01qd_lower() {
        // Lower: scale only A(i,j) for i >= j.
        let mut a = DMatrix::from_row_slice(3, 3, &[1.0, 0.0, 0.0, 2.0, 3.0, 0.0, 4.0, 5.0, 6.0]);
        let nrows: [i32; 0] = [];
        assert_eq!(
            mb01qd(Mb01QdType::Lower, 3, 3, 0, 0, 1.0, 2.0, 0, &nrows, &mut a),
            0
        );
        assert_eq!(a[(0, 0)], 2.0);
        assert_eq!(a[(1, 0)], 4.0);
        assert_eq!(a[(1, 1)], 6.0);
        assert_eq!(a[(2, 2)], 12.0);
        assert_eq!(a[(0, 1)], 0.0); // unchanged (above diagonal)
    }

    #[test]
    fn test_mb01qd_upper() {
        let mut a = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 0.0, 4.0, 5.0, 0.0, 0.0, 6.0]);
        let nrows: [i32; 0] = [];
        assert_eq!(
            mb01qd(Mb01QdType::Upper, 3, 3, 0, 0, 1.0, 0.5, 0, &nrows, &mut a),
            0
        );
        assert_eq!(a[(0, 0)], 0.5);
        assert_eq!(a[(0, 1)], 1.0);
        assert_eq!(a[(1, 1)], 2.0);
        assert_eq!(a[(2, 2)], 3.0);
    }

    #[test]
    fn test_mb01qd_hessenberg() {
        // Upper Hessenberg: first subdiagonal + diagonal + above.
        let mut a = DMatrix::from_row_slice(3, 3, &[1.0, 2.0, 3.0, 7.0, 4.0, 5.0, 0.0, 8.0, 6.0]);
        let nrows: [i32; 0] = [];
        assert_eq!(
            mb01qd(Mb01QdType::Hessenberg, 3, 3, 0, 0, 1.0, 2.0, 0, &nrows, &mut a),
            0
        );
        assert_eq!(a[(0, 0)], 2.0);
        assert_eq!(a[(1, 0)], 14.0);
        assert_eq!(a[(0, 1)], 4.0);
        assert_eq!(a[(1, 1)], 8.0);
        assert_eq!(a[(2, 1)], 16.0);
        assert_eq!(a[(2, 2)], 12.0);
    }
}
