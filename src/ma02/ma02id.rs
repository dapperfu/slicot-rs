//! MA02ID — Norms of a real skew-Hamiltonian or Hamiltonian matrix (SLICOT MA02ID)
//
// X = [A G; Q A'] (skew-H, G=-G', Q=-Q') or [A G; Q -A] (H, G=G', Q=Q').
// QG: columns 0..n-1 = lower part of Q, columns 1..n = upper part of G.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02IdTyp {
    SkewHamiltonian,
    Hamiltonian,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02IdNorm {
    One,
    Frobenius,
    Infinity,
    MaxAbs,
}

/// Reconstructs Q from QG (lower triangle in cols 0..n) and returns full n×n Q.
fn unpack_q(qg: &DMatrix<f64>, n: usize, skew: bool) -> DMatrix<f64> {
    let mut q = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            q[(i, j)] = qg[(i, j)];
        }
        if skew {
            for j in (i + 1)..n {
                q[(i, j)] = -qg[(j, i)];
            }
        } else {
            for j in (i + 1)..n {
                q[(i, j)] = qg[(j, i)];
            }
        }
    }
    q
}

/// Reconstructs G from QG (upper triangle in cols 1..n) and returns full n×n G.
fn unpack_g(qg: &DMatrix<f64>, n: usize, skew: bool) -> DMatrix<f64> {
    let mut g = DMatrix::zeros(n, n);
    for j in 0..n {
        for i in 0..=j {
            g[(i, j)] = qg[(i, j + 1)];
        }
        if skew {
            for i in (j + 1)..n {
                g[(i, j)] = -qg[(j, i + 1)];
            }
        } else {
            for i in (j + 1)..n {
                g[(i, j)] = qg[(j, i + 1)];
            }
        }
    }
    g
}

/// Returns the specified norm of the Hamiltonian/skew-Hamiltonian matrix.
pub fn ma02id(
    typ: Ma02IdTyp,
    norm: Ma02IdNorm,
    a: &DMatrix<f64>,
    qg: &DMatrix<f64>,
    dwork: &mut [f64],
) -> f64 {
    let n = a.nrows();
    if n == 0 || a.ncols() != n || qg.nrows() != n || qg.ncols() < n + 1 {
        return 0.0;
    }
    let skew = typ == Ma02IdTyp::SkewHamiltonian;
    let q = unpack_q(qg, n, skew);
    let g = unpack_g(qg, n, skew);
    let at = a.transpose();
    let sign_a22 = if skew { 1.0 } else { -1.0 };

    match norm {
        Ma02IdNorm::MaxAbs => {
            let mut m = 0.0_f64;
            for i in 0..n {
                for j in 0..n {
                    m = m.max(a[(i, j)].abs()).max(g[(i, j)].abs()).max(q[(i, j)].abs());
                    m = m.max(if i == j { (sign_a22 * a[(i, j)]).abs() } else { at[(i, j)].abs() });
                }
            }
            m
        }
        Ma02IdNorm::Frobenius => {
            let mut sum = 0.0;
            for i in 0..n {
                for j in 0..n {
                    sum += a[(i, j)] * a[(i, j)] + g[(i, j)] * g[(i, j)] + q[(i, j)] * q[(i, j)];
                    sum += if i == j {
                        (sign_a22 * a[(i, j)]).powi(2)
                    } else {
                        at[(i, j)] * at[(i, j)]
                    };
                }
            }
            sum.sqrt()
        }
        Ma02IdNorm::One | Ma02IdNorm::Infinity => {
            if dwork.len() < 2 * n {
                return 0.0;
            }
            let (col_sums, row_sums) = dwork.split_at_mut(n);
            for j in 0..n {
                let mut s = 0.0;
                for i in 0..n {
                    s += a[(i, j)].abs() + q[(i, j)].abs();
                }
                col_sums[j] = s;
            }
            for j in 0..n {
                let mut s = 0.0;
                for i in 0..n {
                    s += g[(i, j)].abs();
                    s += if i == j { (sign_a22 * a[(i, j)]).abs() } else { at[(i, j)].abs() };
                }
                col_sums[j] = col_sums[j].max(s);
            }
            for i in 0..n {
                row_sums[i] = 0.0;
            }
            for j in 0..n {
                for i in 0..n {
                    row_sums[i] += a[(i, j)].abs() + g[(i, j)].abs();
                }
            }
            for j in 0..n {
                for i in 0..n {
                    row_sums[i] += q[(i, j)].abs();
                    row_sums[i] += if i == j {
                        (sign_a22 * a[(j, i)]).abs()
                    } else {
                        at[(j, i)].abs()
                    };
                }
            }
            let one = col_sums[..n].iter().fold(0.0_f64, |a, &b| a.max(b));
            let inf = row_sums.iter().fold(0.0_f64, |a, &b| a.max(b));
            if norm == Ma02IdNorm::One {
                one
            } else {
                inf
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02id_zero() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let qg = DMatrix::<f64>::zeros(0, 1);
        let mut dwork = vec![0.0; 1];
        assert_eq!(
            ma02id(Ma02IdTyp::SkewHamiltonian, Ma02IdNorm::Frobenius, &a, &qg, &mut dwork),
            0.0
        );
    }

    #[test]
    fn test_ma02id_frobenius_1x1() {
        let a = DMatrix::from_row_slice(1, 1, &[1.0]);
        let qg = DMatrix::from_row_slice(1, 2, &[0.0, 0.0]);
        let mut dwork = vec![0.0; 2];
        let f = ma02id(Ma02IdTyp::SkewHamiltonian, Ma02IdNorm::Frobenius, &a, &qg, &mut dwork);
        assert!((f - 2.0_f64.sqrt()).abs() < 1e-10);
    }
}
