//! SG03BW — Generalized Sylvester equation (SLICOT SG03BW)
//!
//! Solves A'*X*C + E'*X*D = SCALE*Y or A*X*C' + E*X*D' = SCALE*Y
//! with A,E M×M (quasi-triangular/triangular), C,D N×N, N=1 or 2.

use nalgebra::DMatrix;

/// Transpose mode for the equation.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Trans {
    /// Solve A'*X*C + E'*X*D = SCALE*Y
    NoTrans,
    /// Solve A*X*C' + E*X*D' = SCALE*Y
    Trans,
}

/// LU with complete pivoting for small matrices (up to 4×4). Solves MAT*rhs = rhs in place.
fn dgetc2_gesc2(dim: usize, mat: &mut [f64], rhs: &mut [f64], _scale: &mut f64) -> i32 {
    if dim == 0 {
        return 0;
    }
    let n = dim;
    let mut piv2 = vec![0usize; n];
    let mut lu = mat.to_vec();

    for k in 0..n {
        let mut max_val = 0.0_f64;
        let mut ip = k;
        let mut jp = k;
        for i in k..n {
            for j in k..n {
                let v = lu[i * 4 + j].abs();
                if v > max_val {
                    max_val = v;
                    ip = i;
                    jp = j;
                }
            }
        }
        piv2[k] = jp;
        if max_val == 0.0 {
            return 1;
        }
        if ip != k {
            for j in 0..4 {
                lu.swap(ip * 4 + j, k * 4 + j);
            }
            rhs.swap(ip, k);
        }
        if jp != k {
            for i in 0..4 {
                lu.swap(i * 4 + jp, i * 4 + k);
            }
        }
        let pivot = lu[k * 4 + k];
        for i in k + 1..n {
            lu[i * 4 + k] /= pivot;
            for j in k + 1..n {
                lu[i * 4 + j] -= lu[i * 4 + k] * lu[k * 4 + j];
            }
        }
    }

    for k in 0..n {
        for i in k + 1..n {
            rhs[i] -= lu[i * 4 + k] * rhs[k];
        }
    }
    for k in (0..n).rev() {
        rhs[k] /= lu[k * 4 + k];
        for i in 0..k {
            rhs[i] -= lu[i * 4 + k] * rhs[k];
        }
    }

    for k in (0..n).rev() {
        let jp = piv2[k];
        if jp != k {
            rhs.swap(jp, k);
        }
    }
    0
}

/// Solves the generalized Sylvester equation (1) or (2).
///
/// * `trans` — NoTrans: A'*X*C + E'*X*D = SCALE*Y; Trans: A*X*C' + E*X*D' = SCALE*Y
/// * `a` — M×M upper quasitriangular
/// * `c` — N×N (N=1 or 2)
/// * `e` — M×M upper triangular
/// * `d` — N×N
/// * `x` — M×N, on entry Y, on exit solution X (scaled by SCALE)
/// * `scale` — output scale factor, 0 < scale <= 1
///
/// # Returns
/// 0 success; < 0 invalid argument (-i = i-th arg); 1 nearly singular.
pub fn sg03bw(
    trans: Trans,
    a: &DMatrix<f64>,
    c: &DMatrix<f64>,
    e: &DMatrix<f64>,
    d: &DMatrix<f64>,
    x: &mut DMatrix<f64>,
    scale: &mut f64,
) -> i32 {
    const MONE: f64 = -1.0;
    const ONE: f64 = 1.0;
    const ZERO: f64 = 0.0;

    let m = a.nrows();
    let n = c.nrows();
    if n != 1 && n != 2 {
        return -3;
    }
    if a.ncols() != m || e.nrows() != m || e.ncols() != m || c.ncols() != n || d.nrows() != n || d.ncols() != n || x.nrows() != m || x.ncols() != n {
        return -5;
    }

    *scale = ONE;
    if m == 0 {
        return 0;
    }

    let notrns = trans == Trans::NoTrans;
    let mut mat = [0.0_f64; 16];
    let mut rhs = [0.0_f64; 4];
    let mut tm = [[0.0_f64; 2]; 2];

    if notrns {
        let mut me = 0_usize;
        while me != m {
            let (ma, mb) = if me + 1 == m {
                (me, 1)
            } else if a[(me + 1, me)].abs() == ZERO {
                (me, 1)
            } else {
                (me, 2)
            };
            let ma = ma;
            let me_end = me + mb;

            for i in 0..mb {
                let mai = ma + i;
                for j in 0..mb {
                    let maj = ma + j;
                    let val = c[(0, 0)] * a[(maj, mai)];
                    let val = if maj <= mai { val + d[(0, 0)] * e[(maj, mai)] } else { val };
                    mat[i * 4 + j] = val;
                }
                rhs[i] = x[(mai, 0)];
            }
            if n == 2 {
                for i in 0..mb {
                    let mai = ma + i;
                    for j in 0..mb {
                        let maj = ma + j;
                        mat[i * 4 + 2 + j] = c[(0, 1)] * a[(maj, mai)];
                        mat[2 * 4 + i * 4 + j] = c[(1, 0)] * a[(maj, mai)];
                        mat[2 * 4 + i * 4 + 2 + j] = c[(1, 1)] * a[(maj, mai)];
                        if maj <= mai {
                            mat[i * 4 + 2 + j] += d[(0, 1)] * e[(maj, mai)];
                            mat[2 * 4 + i * 4 + j] += d[(1, 0)] * e[(maj, mai)];
                            mat[2 * 4 + i * 4 + 2 + j] += d[(1, 1)] * e[(maj, mai)];
                        }
                    }
                    rhs[2 + i] = x[(mai, 1)];
                }
            }

            let dimmat = mb * n;
            let info1 = dgetc2_gesc2(dimmat, &mut mat, &mut rhs, scale);
            if info1 != 0 {
                return 1;
            }

            for i in 0..mb {
                x[(ma + i, 0)] = rhs[i];
            }
            if n == 2 {
                for i in 0..mb {
                    x[(ma + i, 1)] = rhs[2 + i];
                }
            }

            if me_end < m {
                for i in 0..mb {
                    for j in 0..n {
                        tm[i][j] = (0..n).map(|k| x[(ma + i, k)] * c[(k, j)]).sum();
                    }
                }
                for i in me_end..m {
                    for j in 0..n {
                        x[(i, j)] += MONE * (0..mb).map(|ii| a[(ma + ii, i)] * tm[ii][j]).sum::<f64>();
                    }
                }
                for i in 0..mb {
                    for j in 0..n {
                        tm[i][j] = (0..n).map(|k| x[(ma + i, k)] * d[(k, j)]).sum();
                    }
                }
                for i in me_end..m {
                    for j in 0..n {
                        x[(i, j)] += MONE * (0..mb).map(|ii| e[(ma + ii, i)] * tm[ii][j]).sum::<f64>();
                    }
                }
            }
            me = me_end;
        }
    } else {
        let mut ma = m;
        while ma != 0 {
            let (ma_new, mb) = if ma == 1 {
                (0, 1)
            } else if a[(ma - 1, ma - 2)].abs() == ZERO {
                (ma - 1, 1)
            } else {
                (ma - 2, 2)
            };
            let me_end = ma - 1;

            for i in 0..mb {
                let mai = ma_new + i;
                for j in 0..mb {
                    let maj = ma_new + j;
                    let val = c[(0, 0)] * a[(mai, maj)];
                    let val = if maj >= mai { val + d[(0, 0)] * e[(mai, maj)] } else { val };
                    mat[i * 4 + j] = val;
                }
                rhs[i] = x[(mai, 0)];
            }
            if n == 2 {
                for i in 0..mb {
                    let mai = ma_new + i;
                    for j in 0..mb {
                        let maj = ma_new + j;
                        mat[i * 4 + 2 + j] = c[(1, 0)] * a[(mai, maj)];
                        mat[2 * 4 + i * 4 + j] = c[(0, 1)] * a[(mai, maj)];
                        mat[2 * 4 + i * 4 + 2 + j] = c[(1, 1)] * a[(mai, maj)];
                        if maj >= mai {
                            mat[i * 4 + 2 + j] += d[(1, 0)] * e[(mai, maj)];
                            mat[2 * 4 + i * 4 + j] += d[(0, 1)] * e[(mai, maj)];
                            mat[2 * 4 + i * 4 + 2 + j] += d[(1, 1)] * e[(mai, maj)];
                        }
                    }
                    rhs[2 + i] = x[(mai, 1)];
                }
            }

            let dimmat = mb * n;
            let info1 = dgetc2_gesc2(dimmat, &mut mat, &mut rhs, scale);
            if info1 != 0 {
                return 1;
            }

            for i in 0..mb {
                x[(ma_new + i, 0)] = rhs[i];
            }
            if n == 2 {
                for i in 0..mb {
                    x[(ma_new + i, 1)] = rhs[2 + i];
                }
            }

            if ma_new > 0 {
                for i in 0..mb {
                    for j in 0..n {
                        tm[i][j] = (0..n).map(|k| x[(ma_new + i, k)] * c[(j, k)]).sum();
                    }
                }
                for i in 0..ma_new {
                    for j in 0..n {
                        x[(i, j)] += MONE * (0..mb).map(|ii| a[(i, ma_new + ii)] * tm[ii][j]).sum::<f64>();
                    }
                }
                for i in 0..mb {
                    for j in 0..n {
                        tm[i][j] = (0..n).map(|k| x[(ma_new + i, k)] * d[(j, k)]).sum();
                    }
                }
                for i in 0..ma_new {
                    for j in 0..n {
                        x[(i, j)] += MONE * (0..mb).map(|ii| e[(i, ma_new + ii)] * tm[ii][j]).sum::<f64>();
                    }
                }
            }
            ma = ma_new;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sg03bw_m0() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let e = DMatrix::<f64>::zeros(0, 0);
        let d = DMatrix::from_row_slice(1, 1, &[1.0]);
        let mut x = DMatrix::<f64>::zeros(0, 1);
        let mut scale = 0.0;
        assert_eq!(sg03bw(Trans::NoTrans, &a, &c, &e, &d, &mut x, &mut scale), 0);
        assert_eq!(scale, 1.0);
    }

    #[test]
    fn sg03bw_n1_m1() {
        let a = DMatrix::from_row_slice(1, 1, &[-2.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let e = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[1.0]);
        let y0 = 1.0;
        let mut x = DMatrix::from_row_slice(1, 1, &[y0]);
        let mut scale = 0.0;
        assert_eq!(sg03bw(Trans::NoTrans, &a, &c, &e, &d, &mut x, &mut scale), 0);
        assert!(scale > 0.0 && scale <= 1.0);
        let lhs = a.transpose() * &x * &c + e.transpose() * &x * &d;
        assert!((lhs[(0, 0)] - scale * y0).abs() < 1e-10);
    }
}
