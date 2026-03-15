//! MA02IZ — Norms of complex skew-Hamiltonian or Hamiltonian matrix (SLICOT MA02IZ)
//
// X = [A G; Q A'] (skew-H) or [A G; Q -A] (H). QG stores lower Q and upper G (complex).
// A_re, A_im, QG_re, QG_im column-major.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02IzTyp {
    SkewHamiltonian,
    Hamiltonian,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02IzNorm {
    One,
    Frobenius,
    Infinity,
    MaxAbs,
}

fn cabs(re: f64, im: f64) -> f64 {
    (re * re + im * im).sqrt()
}

/// Unpack complex Q from QG (lower triangle cols 0..n).
fn unpack_q_re_im(
    qg_re: &[f64],
    qg_im: &[f64],
    ldqg: usize,
    n: usize,
    skew: bool,
) -> (Vec<f64>, Vec<f64>) {
    let mut q_re = vec![0.0; n * n];
    let mut q_im = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            q_re[i + j * n] = qg_re[i + j * ldqg];
            q_im[i + j * n] = qg_im[i + j * ldqg];
        }
        for j in (i + 1)..n {
            if skew {
                q_re[i + j * n] = -qg_re[j + i * ldqg];
                q_im[i + j * n] = -qg_im[j + i * ldqg];
            } else {
                q_re[i + j * n] = qg_re[j + i * ldqg];
                q_im[i + j * n] = qg_im[j + i * ldqg];
            }
        }
    }
    (q_re, q_im)
}

/// Unpack complex G from QG (upper triangle cols 1..n).
fn unpack_g_re_im(
    qg_re: &[f64],
    qg_im: &[f64],
    ldqg: usize,
    n: usize,
    skew: bool,
) -> (Vec<f64>, Vec<f64>) {
    let mut g_re = vec![0.0; n * n];
    let mut g_im = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..=j {
            g_re[i + j * n] = qg_re[i + (j + 1) * ldqg];
            g_im[i + j * n] = qg_im[i + (j + 1) * ldqg];
        }
        for i in (j + 1)..n {
            if skew {
                g_re[i + j * n] = -qg_re[j + (i + 1) * ldqg];
                g_im[i + j * n] = -qg_im[j + (i + 1) * ldqg];
            } else {
                g_re[i + j * n] = qg_re[j + (i + 1) * ldqg];
                g_im[i + j * n] = qg_im[j + (i + 1) * ldqg];
            }
        }
    }
    (g_re, g_im)
}

/// Returns the specified norm of the complex Hamiltonian/skew-Hamiltonian matrix.
pub fn ma02iz(
    typ: Ma02IzTyp,
    norm: Ma02IzNorm,
    n: usize,
    a_re: &[f64],
    a_im: &[f64],
    lda: usize,
    qg_re: &[f64],
    qg_im: &[f64],
    ldqg: usize,
    dwork: &mut [f64],
) -> f64 {
    if n == 0 || lda < n || ldqg < n {
        return 0.0;
    }
    let skew = typ == Ma02IzTyp::SkewHamiltonian;
    let sign_a22 = if skew { 1.0 } else { -1.0 };
    let (q_re, q_im) = unpack_q_re_im(qg_re, qg_im, ldqg, n, skew);
    let (g_re, g_im) = unpack_g_re_im(qg_re, qg_im, ldqg, n, skew);

    match norm {
        Ma02IzNorm::MaxAbs => {
            let mut m = 0.0_f64;
            for i in 0..n {
                for j in 0..n {
                    m = m.max(cabs(a_re[i + j * lda], a_im[i + j * lda]));
                    m = m.max(cabs(q_re[i + j * n], q_im[i + j * n]));
                    m = m.max(cabs(g_re[i + j * n], g_im[i + j * n]));
                    let a22_re = if i == j { sign_a22 * a_re[i + j * lda] } else { a_re[j + i * lda] };
                    let a22_im = if i == j { sign_a22 * a_im[i + j * lda] } else { -a_im[j + i * lda] };
                    m = m.max(cabs(a22_re, a22_im));
                }
            }
            m
        }
        Ma02IzNorm::Frobenius => {
            let mut sum = 0.0_f64;
            for i in 0..n {
                for j in 0..n {
                    sum += cabs(a_re[i + j * lda], a_im[i + j * lda]).powi(2);
                    sum += cabs(q_re[i + j * n], q_im[i + j * n]).powi(2);
                    sum += cabs(g_re[i + j * n], g_im[i + j * n]).powi(2);
                    let a22_re = if i == j { sign_a22 * a_re[i + j * lda] } else { a_re[j + i * lda] };
                    let a22_im = if i == j { sign_a22 * a_im[i + j * lda] } else { -a_im[j + i * lda] };
                    sum += cabs(a22_re, a22_im).powi(2);
                }
            }
            sum.sqrt()
        }
        Ma02IzNorm::One | Ma02IzNorm::Infinity => {
            if dwork.len() < 2 * n {
                return 0.0;
            }
            let (col_sums, row_sums) = dwork.split_at_mut(n);
            for j in 0..n {
                let mut s = 0.0;
                for i in 0..n {
                    s += cabs(a_re[i + j * lda], a_im[i + j * lda]) + cabs(q_re[i + j * n], q_im[i + j * n]);
                }
                col_sums[j] = s;
            }
            for j in 0..n {
                let mut s = 0.0;
                for i in 0..n {
                    s += cabs(g_re[i + j * n], g_im[i + j * n]);
                    let a22_re = if i == j { sign_a22 * a_re[i + j * lda] } else { a_re[j + i * lda] };
                    let a22_im = if i == j { sign_a22 * a_im[i + j * lda] } else { -a_im[j + i * lda] };
                    s += cabs(a22_re, a22_im);
                }
                col_sums[j] = col_sums[j].max(s);
            }
            for i in 0..n {
                row_sums[i] = 0.0;
            }
            for j in 0..n {
                for i in 0..n {
                    row_sums[i] += cabs(a_re[i + j * lda], a_im[i + j * lda]) + cabs(g_re[i + j * n], g_im[i + j * n]);
                }
            }
            for j in 0..n {
                for i in 0..n {
                    row_sums[i] += cabs(q_re[i + j * n], q_im[i + j * n]);
                    let a22_re = if i == j { sign_a22 * a_re[j + i * lda] } else { a_re[i + j * lda] };
                    let a22_im = if i == j { sign_a22 * a_im[j + i * lda] } else { -a_im[i + j * lda] };
                    row_sums[i] += cabs(a22_re, a22_im);
                }
            }
            let one = col_sums[..n].iter().fold(0.0_f64, |a, &b| a.max(b));
            let inf = row_sums.iter().fold(0.0_f64, |a, &b| a.max(b));
            if norm == Ma02IzNorm::One {
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
    fn test_ma02iz_zero() {
        let mut dwork = vec![0.0; 2];
        assert_eq!(
            ma02iz(
                Ma02IzTyp::SkewHamiltonian,
                Ma02IzNorm::Frobenius,
                0,
                &[],
                &[],
                0,
                &[],
                &[],
                0,
                &mut dwork,
            ),
            0.0
        );
    }

    #[test]
    fn test_ma02iz_frobenius_1x1() {
        let a_re = vec![1.0];
        let a_im = vec![0.0];
        // QG is n×(n+1) = 1×2, column-major ldqg=1
        let qg_re = vec![0.0, 0.0];
        let qg_im = vec![0.0, 0.0];
        let mut dwork = vec![0.0; 2];
        let f = ma02iz(
            Ma02IzTyp::SkewHamiltonian,
            Ma02IzNorm::Frobenius,
            1,
            &a_re,
            &a_im,
            1,
            &qg_re,
            &qg_im,
            1, // ldqg >= n; need (n+1) columns => 2 elements
            &mut dwork,
        );
        assert!((f - 2.0_f64.sqrt()).abs() < 1e-10);
    }
}
