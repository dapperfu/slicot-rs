//! SB10UD — Normalization of a system for H2 controller design.
//! Reduces D12 and D21 to unit diagonal form and transforms B, C accordingly.

use nalgebra::DMatrix;

const EPS: f64 = 2.22e-16_f64;

/// Normalizes D12 to [0; I] and D21 to [0 I] via SVD, and updates B, C, D.
/// On entry: B (N×M), C (NP×N), D (NP×M). NCON = M2, NMEAS = NP2.
/// On exit: B and C are transformed; D12/D21 blocks overwritten; TU (M2×M2), TY (NP2×NP2) output.
///
/// # Returns
/// 0 success; 1 D12 not full column rank; 2 D21 not full row rank; 3 SVD did not converge; < 0 invalid argument.
pub fn sb10ud(
    n: usize,
    m: usize,
    np: usize,
    ncon: usize,
    nmeas: usize,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &mut DMatrix<f64>,
    tu: &mut DMatrix<f64>,
    ty: &mut DMatrix<f64>,
    rcond: &mut [f64],
    tol: f64,
) -> i32 {
    if b.nrows() != n || b.ncols() != m {
        return -7;
    }
    if c.nrows() != np || c.ncols() != n {
        return -9;
    }
    if d.nrows() != np || d.ncols() != m {
        return -11;
    }
    if ncon > m || nmeas > np || np - nmeas < ncon || m - ncon < nmeas {
        return -5;
    }
    let m1 = m - ncon;
    let np1 = np - nmeas;
    let m2 = ncon;
    let np2 = nmeas;
    if tu.nrows() != m2 || tu.ncols() != m2 {
        return -13;
    }
    if ty.nrows() != np2 || ty.ncols() != np2 {
        return -15;
    }
    if rcond.len() < 2 {
        return -17;
    }
    let tol_use = if tol > 0.0 { tol } else { EPS.sqrt() };
    rcond[0] = 0.0;
    rcond[1] = 0.0;

    if n == 0 || m2 == 0 || np2 == 0 {
        for i in 0..m2 {
            for j in 0..m2 {
                tu[(i, j)] = if i == j { 1.0 } else { 0.0 };
            }
        }
        for i in 0..np2 {
            for j in 0..np2 {
                ty[(i, j)] = if i == j { 1.0 } else { 0.0 };
            }
        }
        rcond[0] = 1.0;
        rcond[1] = 1.0;
        return 0;
    }

    // --- D12: NP1×M2, full column rank. D12 = U*S*V'. TU = V*inv(S). Then B2 <- B2*TU, and L*[C1;D12] so D12 <- [0;I]. L = [0;I]*U'.
    let mut d12 = DMatrix::zeros(np1, m2);
    for i in 0..np1 {
        for j in 0..m2 {
            d12[(i, j)] = d[(i, m1 + j)];
        }
    }
    let svd12 = d12.svd(true, true);
    let s12 = &svd12.singular_values;
    let min_s12 = s12.iter().cloned().fold(0.0_f64, f64::max);
    let rank12 = s12.iter().filter(|x| **x > tol_use).count();
    if rank12 < m2 {
        rcond[0] = if s12[0] > 0.0 { min_s12 / s12[0] } else { 0.0 };
        return 1;
    }
    let u12 = match &svd12.u {
        Some(u) => u,
        None => return 3,
    };
    let v_t12 = match &svd12.v_t {
        Some(vt) => vt,
        None => return 3,
    };
    let v12 = v_t12.transpose();
    rcond[0] = if s12[0] > 0.0 { s12[m2 - 1] / s12[0] } else { 0.0 };
    if rcond[0] < tol_use {
        return 1;
    }
    for i in 0..m2 {
        for j in 0..m2 {
            let s_inv = if s12[j] > tol_use { 1.0 / s12[j] } else { 0.0 };
            tu[(i, j)] = v12[(i, j)] * s_inv;
        }
    }
    let mut b2_new = DMatrix::zeros(n, m2);
    for i in 0..n {
        for j in 0..m2 {
            let mut sum = 0.0;
            for k in 0..m2 {
                sum += b[(i, m1 + k)] * tu[(k, j)];
            }
            b2_new[(i, j)] = sum;
        }
    }
    for i in 0..n {
        for j in 0..m2 {
            b[(i, m1 + j)] = b2_new[(i, j)];
        }
    }
    let mut l1 = DMatrix::zeros(np1, np1);
    for i in (np1 - m2)..np1 {
        for j in 0..np1 {
            l1[(i, j)] = u12[(j, i)];
        }
    }
    let mut c1 = DMatrix::zeros(np1, n);
    for i in 0..np1 {
        for j in 0..n {
            c1[(i, j)] = c[(i, j)];
        }
    }
    let c1_new = &l1 * &c1;
    for i in 0..np1 {
        for j in 0..n {
            c[(i, j)] = c1_new[(i, j)];
        }
    }
    for i in 0..np1 {
        for j in 0..m2 {
            d[(i, m1 + j)] = if i >= np1 - m2 && (i - (np1 - m2)) == j {
                1.0
            } else {
                0.0
            };
        }
    }

    // --- D21: NP2×M1, full row rank. TY*D21 = [0 I]. TY = [0 I]*pinv(D21). Then C2 <- TY*C2, D21 <- [0 I].
    if m1 > 0 {
        let mut d21 = DMatrix::zeros(np2, m1);
        for i in 0..np2 {
            for j in 0..m1 {
                d21[(i, j)] = d[(np1 + i, j)];
            }
        }
        let svd21 = d21.svd(true, true);
        let s21 = &svd21.singular_values;
        let rank21 = s21.iter().filter(|x| **x > tol_use).count();
        if rank21 < m1 {
            rcond[1] = if s21.len() > 0 && s21[0] > 0.0 {
                s21[s21.len() - 1] / s21[0]
            } else {
                0.0
            };
            return 2;
        }
        let u21 = match &svd21.u {
            Some(u) => u,
            None => return 3,
        };
        let v_t21 = match &svd21.v_t {
            Some(vt) => vt,
            None => return 3,
        };
        let v21 = v_t21.transpose();
        rcond[1] = if s21[0] > 0.0 {
            s21[s21.len().saturating_sub(1)] / s21[0]
        } else {
            0.0
        };
        if rcond[1] < tol_use {
            return 2;
        }
        for i in 0..np2 {
            for j in 0..np2 {
                ty[(i, j)] = 0.0;
            }
        }
        for i in 0..m1.min(np2) {
            for j in 0..np2 {
                let mut sum = 0.0;
                for k in 0..m1.min(np2) {
                    let s_inv = if k < s21.len() && s21[k] > tol_use {
                        1.0 / s21[k]
                    } else {
                        0.0
                    };
                    sum += v21[(k, i)] * s_inv * u21[(j, k)];
                }
                if np2 - m1 + i < np2 {
                    ty[(np2 - m1 + i, j)] = sum;
                }
            }
        }
        let mut c2 = DMatrix::zeros(np2, n);
        for i in 0..np2 {
            for j in 0..n {
                c2[(i, j)] = c[(np1 + i, j)];
            }
        }
        let c2_new = &*ty * &c2;
        for i in 0..np2 {
            for j in 0..n {
                c[(np1 + i, j)] = c2_new[(i, j)];
            }
        }
        for i in 0..np2 {
            for j in 0..m1 {
                d[(np1 + i, j)] = if i >= np2 - m1 && (i - (np2 - m1)) == j {
                    1.0
                } else {
                    0.0
                };
            }
        }
    } else {
        for i in 0..np2 {
            for j in 0..np2 {
                ty[(i, j)] = if i == j { 1.0 } else { 0.0 };
            }
        }
        rcond[1] = 1.0;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb10ud_n0() {
        let mut b = DMatrix::zeros(0, 2);
        let mut c = DMatrix::zeros(2, 0);
        let mut d = DMatrix::zeros(2, 2);
        let mut tu = DMatrix::zeros(1, 1);
        let mut ty = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 2];
        assert_eq!(
            sb10ud(0, 2, 2, 1, 1, &mut b, &mut c, &mut d, &mut tu, &mut ty, &mut rcond, 1e-10),
            0
        );
    }

    #[test]
    fn test_sb10ud_simple() {
        let n = 1usize;
        let m = 2usize;
        let np = 2usize;
        let ncon = 1usize;
        let nmeas = 1usize;
        let mut b = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let mut c = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let mut d = DMatrix::from_row_slice(2, 2, &[0.0, 0.0, 0.0, 1.0]);
        let mut tu = DMatrix::zeros(1, 1);
        let mut ty = DMatrix::zeros(1, 1);
        let mut rcond = [0.0; 2];
        let info = sb10ud(n, m, np, ncon, nmeas, &mut b, &mut c, &mut d, &mut tu, &mut ty, &mut rcond, 1e-10);
        assert!(info == 0, "info = {}", info);
        assert!((d[(1, 1)] - 1.0).abs() < 1e-10);
    }
}
