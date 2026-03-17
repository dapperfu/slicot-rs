//! AB13HD — L-infinity norm of descriptor system (SLICOT).
//!
//! Full SLICOT-equivalent API. Continuous/discrete, standard or descriptor.
//! Simplified path implemented: JOBE=Identity, EQUIL=NoScale, CKPROP=NoCheck;
//! frequency sweep for standard system (continuous or discrete).

use nalgebra::DMatrix;
use num_complex::Complex64;
use std::f64::consts::PI;

/// System type: continuous or discrete.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// E matrix form.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobE {
    Identity,
    General,
    Compressed,
}

/// Equilibration.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    Scale,
    NoScale,
}

/// D matrix presence.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobD {
    Present,
    Zero,
    FullRank,
}

/// Check properness (descriptor).
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum CkProp {
    Check,
    NoCheck,
}

/// Reduce order before norm.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Reduce {
    Reduce,
    NoReduce,
}

/// Poles/midpoints selection.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Poles {
    All,
    Part,
}

/// Full AB13HD: L-infinity norm of (A-lambda*E,B,C,D).
///
/// # Returns
/// INFO: 0 = success (quick return); 1 = not implemented (unsupported options); <0 = invalid argument.
/// Outputs: NR (reduced order), GPEAK (norm, frequency), FPEAK (frequency), IWARN.
///
/// Simplified path: only JOBE=Identity, EQUIL=NoScale, CKPROP=NoCheck are implemented.
pub fn ab13hd(
    dico: Dico,
    jobe: JobE,
    equil: Equil,
    jobd: JobD,
    ckprop: CkProp,
    _reduce: Reduce,
    _poles: Poles,
    n: usize,
    m: usize,
    p: usize,
    _ranke: usize,
    fpeak: &mut [f64],
    a: &DMatrix<f64>,
    e: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    nr: &mut usize,
    gpeak: &mut [f64],
    _tol: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
    _ldwork: i32,
    _zwork: &mut [num_complex::Complex64],
    _lzwork: i32,
    _bwork: &mut [bool],
    iwarn: &mut i32,
) -> i32 {
    if n == 0 && m == 0 && p == 0 {
        *nr = 0;
        if gpeak.len() >= 2 {
            gpeak[0] = 0.0;
            gpeak[1] = 0.0;
        }
        if fpeak.len() >= 2 {
            fpeak[0] = 0.0;
            fpeak[1] = 0.0;
        }
        *iwarn = 0;
        return 0;
    }
    if n > 0 && (a.nrows() != n || a.ncols() != n) {
        return -12;
    }
    if n > 0 && (e.nrows() != n || e.ncols() != n) {
        return -14;
    }
    if n > 0 && m > 0 && (b.nrows() != n || b.ncols() != m) {
        return -16;
    }
    if p > 0 && n > 0 && (c.nrows() != p || c.ncols() != n) {
        return -18;
    }
    if p > 0 && m > 0 && (d.nrows() != p || d.ncols() != m) {
        return -20;
    }

    // Simplified path: only standard system, no equilibration, no properness check.
    if jobe != JobE::Identity || equil != Equil::NoScale || ckprop != CkProp::NoCheck {
        return 1;
    }

    *nr = n;
    *iwarn = 0;

    // N = 0: static system; norm = sigma_max(D) or 0
    if n == 0 {
        if matches!(jobd, JobD::Zero) || (p == 0 || m == 0) {
            if gpeak.len() >= 1 {
                gpeak[0] = 0.0;
            }
            if gpeak.len() >= 2 {
                gpeak[1] = 0.0;
            }
            if fpeak.len() >= 1 {
                fpeak[0] = 0.0;
            }
            if fpeak.len() >= 2 {
                fpeak[1] = 1.0;
            }
            return 0;
        }
        let svd = d.clone().svd(false, false);
        let sigma_max = if svd.singular_values.len() > 0 {
            svd.singular_values[0]
        } else {
            0.0
        };
        if gpeak.len() >= 1 {
            gpeak[0] = sigma_max;
        }
        if gpeak.len() >= 2 {
            gpeak[1] = 1.0;
        }
        if fpeak.len() >= 1 {
            fpeak[0] = 0.0;
        }
        if fpeak.len() >= 2 {
            fpeak[1] = 0.0; // infinite frequency for static
        }
        return 0;
    }

    // Frequency sweep: G(freq) = C*inv(freq*I - A)*B + D, norm = max_omega sigma_max(G)
    let anrm = a.norm();
    let (omega_max, n_grid) = match dico {
        Dico::Continuous => ((anrm * 2.0 + 1.0).max(1.0), 128),
        Dico::Discrete => (2.0 * PI, 128),
    };
    let mut best_norm = 0.0_f64;
    let mut best_omega = 0.0_f64;

    let mut a_work = a.clone();
    let mut g_complex = DMatrix::<Complex64>::zeros(p, m);

    for k in 0..=n_grid {
        let omega = (k as f64 / n_grid as f64) * omega_max;
        let freq = match dico {
            Dico::Continuous => Complex64::new(0.0, omega),
            Dico::Discrete => Complex64::new(omega.cos(), omega.sin()),
        };
        let info = crate::tb05::tb05ad::tb05ad(
            crate::tb05::tb05ad::Baleig::N,
            crate::tb05::tb05ad::Inita::G,
            &mut a_work,
            b,
            c,
            freq,
            &mut g_complex,
            None,
            None,
            None,
            None,
        );
        if info != 0 {
            continue; // skip singular frequency
        }
        if matches!(jobd, JobD::Present | JobD::FullRank) && p > 0 && m > 0 {
            for i in 0..p {
                for j in 0..m {
                    g_complex[(i, j)] += Complex64::new(d[(i, j)], 0.0);
                }
            }
        }
        let svd = g_complex.clone().svd(false, false);
        let sigma = if svd.singular_values.len() > 0 {
            svd.singular_values[0]
        } else {
            0.0
        };
        if sigma > best_norm {
            best_norm = sigma;
            best_omega = omega;
        }
    }

    if gpeak.len() >= 1 {
        gpeak[0] = best_norm;
    }
    if gpeak.len() >= 2 {
        gpeak[1] = 1.0;
    }
    if fpeak.len() >= 1 {
        fpeak[0] = best_omega;
    }
    if fpeak.len() >= 2 {
        fpeak[1] = 1.0;
    }
    0
}

/// Compatibility wrapper: (n, m) -> INFO. P = m for square D.
#[inline]
pub fn ab13hd_nm(n: usize, m: usize) -> i32 {
    let p = m;
    if n == 0 && m == 0 {
        return 0;
    }
    let a = DMatrix::zeros(n.max(1), n.max(1));
    let e = DMatrix::zeros(n.max(1), n.max(1));
    let b = DMatrix::zeros(n.max(1), m.max(1));
    let c = DMatrix::zeros(p.max(1), n.max(1));
    let d = DMatrix::zeros(p.max(1), m.max(1));
    let mut nr = 0_usize;
    let mut gpeak = [0.0_f64; 2];
    let mut fpeak = [0.0_f64; 2];
    let mut iwarn = 0i32;
    ab13hd(
        Dico::Continuous,
        JobE::Identity,
        Equil::NoScale,
        JobD::Zero,
        CkProp::NoCheck,
        Reduce::NoReduce,
        Poles::All,
        n,
        m,
        p,
        0,
        &mut fpeak,
        &a,
        &e,
        &b,
        &c,
        &d,
        &mut nr,
        &mut gpeak,
        0.0,
        &mut [],
        &mut [],
        0,
        &mut [],
        0,
        &mut [],
        &mut iwarn,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13hd_trivial() {
        let a = DMatrix::zeros(0, 0);
        let e = DMatrix::zeros(0, 0);
        let b = DMatrix::zeros(0, 0);
        let c = DMatrix::zeros(0, 0);
        let d = DMatrix::zeros(0, 0);
        let mut nr = 0_usize;
        let mut gpeak = [0.0; 2];
        let mut fpeak = [0.0; 2];
        let mut iwarn = 0i32;
        assert_eq!(
            ab13hd(
                Dico::Continuous,
                JobE::Identity,
                Equil::NoScale,
                JobD::Zero,
                CkProp::NoCheck,
                Reduce::NoReduce,
                Poles::All,
                0,
                0,
                0,
                0,
                &mut fpeak,
                &a,
                &e,
                &b,
                &c,
                &d,
                &mut nr,
                &mut gpeak,
                0.0,
                &mut [],
                &mut [],
                0,
                &mut [],
                0,
                &mut [],
                &mut iwarn,
            ),
            0
        );
        assert_eq!(nr, 0);
    }

    #[test]
    fn test_ab13hd_continuous_n1() {
        // G(s) = 1/(s+1) => ||G||_inf = 1 at omega=0
        let a = DMatrix::from_row_slice(1, 1, &[-1.0]);
        let e = DMatrix::from_row_slice(1, 1, &[1.0]);
        let b = DMatrix::from_row_slice(1, 1, &[1.0]);
        let c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut nr = 0_usize;
        let mut gpeak = [0.0; 2];
        let mut fpeak = [0.0; 2];
        let mut iwarn = 0i32;
        let info = ab13hd(
            Dico::Continuous,
            JobE::Identity,
            Equil::NoScale,
            JobD::Zero,
            CkProp::NoCheck,
            Reduce::NoReduce,
            Poles::All,
            1,
            1,
            1,
            0,
            &mut fpeak,
            &a,
            &e,
            &b,
            &c,
            &d,
            &mut nr,
            &mut gpeak,
            0.0,
            &mut [],
            &mut [],
            0,
            &mut [],
            0,
            &mut [],
            &mut iwarn,
        );
        assert_eq!(info, 0);
        assert_eq!(nr, 1);
        assert!(gpeak[0] > 0.9 && gpeak[0] <= 1.1);
    }
}
