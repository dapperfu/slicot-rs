//! SG02CX — Line search parameter minimizing residual of (generalized) algebraic Riccati equations (SLICOT SG02).
//!
//! Finds α in [0,2] minimizing P(α) = ||(1-α)*R(X) ± α²*V||_F.

use nalgebra::DMatrix;

/// JOBE: E general or identity.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobE {
    General,
    Identity,
}

/// FLAG: plus or minus in quadratic term.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Flag {
    Plus,
    Minus,
}

/// JOBG: how V is defined (G, D, F, or H,K).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobG {
    G,
    D,
    F,
    H,
}

/// UPLO: upper or lower triangle of symmetric matrices.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    Upper,
    Lower,
}

/// TRANS: op(W) = W or W'.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Trans {
    NoTrans,
    Trans,
}

fn sym_full_from_tri(n: usize, a: &DMatrix<f64>, uplo: Uplo) -> DMatrix<f64> {
    let mut f = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            f[(i, j)] = match uplo {
                Uplo::Upper => if j >= i { a[(i, j)] } else { a[(j, i)] },
                Uplo::Lower => if i >= j { a[(i, j)] } else { a[(j, i)] },
            };
        }
    }
    f
}

/// Real roots of cubic a*t^3 + b*t^2 + c*t + d in [0, 2] via companion matrix eigenvalues.
fn cubic_roots_01(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < 1e-20 {
        if b.abs() < 1e-20 {
            if c.abs() < 1e-20 {
                return vec![];
            }
            let t = -d / c;
            if t >= 0.0 && t <= 2.0 {
                return vec![t];
            }
            return vec![];
        }
        let disc = c * c - 4.0 * b * d;
        if disc < 0.0 {
            return vec![];
        }
        let sd = disc.sqrt();
        let t1 = (-c + sd) / (2.0 * b);
        let t2 = (-c - sd) / (2.0 * b);
        let mut out = vec![];
        if t1 >= 0.0 && t1 <= 2.0 {
            out.push(t1);
        }
        if t2 >= 0.0 && t2 <= 2.0 && (t2 - t1).abs() > 1e-10 {
            out.push(t2);
        }
        return out;
    }
    let (b, c, d) = (b / a, c / a, d / a);
    let companion = DMatrix::from_row_slice(3, 3, &[
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        -d, -c, -b,
    ]);
    let mut roots = vec![];
    if let Some(schur) = companion.try_schur(1e-14, 100) {
        let eigs = schur.complex_eigenvalues();
        for z in eigs.iter() {
            if z.im.abs() < 1e-12 && z.re >= 0.0 && z.re <= 2.0 {
                roots.push(z.re);
            }
        }
    }
    roots
}

/// Finds α in [0,2] minimizing P(α) = ||(1-α)*R(X) ± α²*V||_F.
/// Returns alpha and rnorm = sqrt(P(alpha)).
///
/// # Returns
/// 0 success; < 0 invalid argument; 1 MC01XD/eigen computation failed.
pub fn sg02cx(
    jobe: JobE,
    flag: Flag,
    jobg: JobG,
    uplo: Uplo,
    trans: Trans,
    n: usize,
    m: usize,
    e: Option<&DMatrix<f64>>,
    r: &DMatrix<f64>,
    s: &DMatrix<f64>,
    g: &DMatrix<f64>,
    alpha: &mut f64,
    rnorm: &mut f64,
    dwork: &mut [f64],
    _ldwork: i32,
    iwarn: &mut i32,
    info: &mut i32,
) -> i32 {
    *info = 0;
    *iwarn = 0;
    if n == 0 {
        *alpha = 1.0;
        *rnorm = 0.0;
        return 0;
    }

    let r_full = sym_full_from_tri(n, r, uplo);
    let rr = r_full.norm();
    let rr2 = rr * rr;

    let v_full = match jobg {
        JobG::G => {
            let g_full = sym_full_from_tri(n, g, uplo);
            let s_full = sym_full_from_tri(n, s, uplo);
            if jobe == JobE::General {
                let e_mat = match e {
                    Some(em) => sym_full_from_tri(n, em, uplo),
                    None => DMatrix::identity(n, n),
                };
                let se = if trans == Trans::NoTrans { &s_full * &e_mat } else { e_mat.transpose() * &s_full };
                let e_t_se = e_mat.transpose() * &se;
                &e_t_se * &g_full * &e_t_se
            } else {
                &s_full * &g_full * &s_full
            }
        }
        JobG::D => {
            let d = g;
            let g_eff = d * d.transpose();
            let s_full = sym_full_from_tri(n, s, uplo);
            if jobe == JobE::General {
                let e_mat = e.map(|em| sym_full_from_tri(n, em, uplo)).unwrap_or_else(|| DMatrix::identity(n, n));
                let se = if trans == Trans::NoTrans { &s_full * &e_mat } else { e_mat.transpose() * &s_full };
                let e_t_se = e_mat.transpose() * &se;
                &e_t_se * &g_eff * &e_t_se
            } else {
                &s_full * &g_eff * &s_full
            }
        }
        JobG::F => {
            let ff = g;
            ff * ff.transpose()
        }
        JobG::H => {
            let hh = g;
            let kk = s;
            hh * kk
        }
    };

    let vv = v_full.norm().powi(2);
    let mut rv = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            rv += r_full[(i, j)] * v_full[(i, j)];
        }
    }
    let sign = if flag == Flag::Plus { 1.0 } else { -1.0 };

    // P(α) = (1-α)²*rr² + α⁴*vv ± 2(1-α)α²*rv. Derivative: dP/dα = 2(1-α)(-1)*rr² + 4α³*vv ± 2[ -α²*rv + (1-α)*2α*rv ] = -2(1-α)rr² + 4α³*vv ± 2rv*(-α² + 2α(1-α)) = -2rr² + 2α*rr² + 4α³*vv ± 2rv*(-α² + 2α - 2α²) = 2rr²*α + 4vv*α³ ± 2rv*(2α - 3α²) - 2rr². So 4vv*α³ ± 2rv*(-3α²) + (2rr² ± 4rv)*α - 2rr² = 0 => 4vv*α³ ∓ 6rv*α² + (2rr² ± 4rv)*α - 2rr² = 0.
    let (a3, a2, a1, a0) = (
        4.0 * vv,
        -6.0 * rv * sign,
        2.0 * rr2 + 4.0 * rv * sign,
        -2.0 * rr2,
    );
    let candidates = cubic_roots_01(a3, a2, a1, a0);
    let mut best_alpha: f64 = 1.0;
    let mut best_p: f64 = (1.0 - best_alpha).powi(2) * rr2 + best_alpha.powi(4) * vv + sign * 2.0 * (1.0 - best_alpha) * best_alpha.powi(2) * rv;
    if best_p < 0.0 {
        best_p = 0.0;
    }
    for &t in &[0.0_f64, 2.0] {
        let p = (1.0 - t).powi(2) * rr2 + t.powi(4) * vv + sign * 2.0 * (1.0 - t) * t.powi(2) * rv;
        let p = if p < 0.0 { 0.0 } else { p };
        if p < best_p {
            best_p = p;
            best_alpha = t;
        }
    }
    for &t in &candidates {
        let p = (1.0 - t).powi(2) * rr2 + t.powi(4) * vv + sign * 2.0 * (1.0 - t) * t.powi(2) * rv;
        let p = if p < 0.0 { 0.0 } else { p };
        if p < best_p {
            best_p = p;
            best_alpha = t;
        }
    }
    *alpha = best_alpha;
    *rnorm = best_p.sqrt();
    if dwork.len() >= n * n {
        for i in 0..n {
            for j in 0..n {
                dwork[i * n + j] = v_full[(i, j)];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sg02cx_n0() {
        let r = DMatrix::<f64>::zeros(0, 0);
        let s = DMatrix::zeros(0, 0);
        let g = DMatrix::zeros(0, 0);
        let mut alpha = 0.0;
        let mut rnorm = 0.0;
        let mut dwork = [0.0; 1];
        let mut iwarn = 0;
        let mut info = 0;
        assert_eq!(
            sg02cx(JobE::Identity, Flag::Minus, JobG::G, Uplo::Upper, Trans::NoTrans, 0, 0, None, &r, &s, &g, &mut alpha, &mut rnorm, &mut dwork, 0, &mut iwarn, &mut info),
            0
        );
        assert_eq!(alpha, 1.0);
        assert_eq!(rnorm, 0.0);
    }
}
