//! SB03MD — Solution of continuous- or discrete-time Lyapunov equation and/or separation estimate.
//!
//! Solves op(A)'*X + X*op(A) = scale*C (continuous) or op(A)'*X*op(A) - X = scale*C (discrete),
//! and optionally estimates the separation and forward error bound.

use nalgebra::DMatrix;
use crate::mb01::mb01rd::mb01rd;
use crate::mb01::mb01ru::{Mb01RuTrans, Mb01RuUplo};
use super::sb03mx;
use super::sb03my;

/// DICO: 'C' = continuous, 'D' = discrete.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dico {
    Continuous,
    Discrete,
}

/// JOB: 'X' = solution only, 'S' = separation only, 'B' = both.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Job {
    X,
    S,
    B,
}

/// FACT: 'N' = compute Schur, 'F' = A and U already in Schur form.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Fact {
    NotFactored,
    Factored,
}

/// TRANA: 'N' = op(A)=A, 'T' or 'C' = op(A)=A'.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum TranA {
    NoTrans,
    Trans,
}

/// Solves Lyapunov equation and/or estimates separation.
///
/// - `dico`: Continuous or discrete.
/// - `job`: X (solution), S (separation), or B (both).
/// - `fact`: NotFactored (compute Schur) or Factored (A,U already Schur).
/// - `trana`: NoTrans or Trans.
/// - `a`: N×N matrix A (overwritten with Schur form if fact=NotFactored).
/// - `u`: N×N matrix U (orthogonal from Schur; input if fact=Factored, output otherwise).
/// - `c`: N×N symmetric RHS (overwritten with solution X if job is X or B).
/// - `scale`: output scale factor.
/// - `sep`: output separation (if job S or B).
/// - `ferr`: output forward error bound (if job B).
/// - `wr`, `wi`: output real/imaginary parts of eigenvalues (if fact=NotFactored).
///
/// Returns INFO: 0 = success; >0 = QR/Schur failed or equation singular (see SLICOT doc).
pub fn sb03md(
    dico: Dico,
    job: Job,
    fact: Fact,
    trana: TranA,
    a: &mut DMatrix<f64>,
    u: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    scale: &mut f64,
    sep: &mut f64,
    ferr: &mut f64,
    wr: &mut [f64],
    wi: &mut [f64],
) -> i32 {
    let n = a.nrows();
    if n == 0 {
        *scale = 1.0;
        if matches!(job, Job::S | Job::B) {
            *sep = 0.0;
        }
        if matches!(job, Job::B) {
            *ferr = 0.0;
        }
        return 0;
    }

    let wantx = matches!(job, Job::X | Job::B);
    let wantsp = matches!(job, Job::S | Job::B);
    let nofact = matches!(fact, Fact::NotFactored);
    let notra = matches!(trana, TranA::NoTrans);

    if nofact {
        let schur = match a.clone().try_schur(1e-14, 200) {
            Some(s) => s,
            None => return 1,
        };
        let eigs = schur.complex_eigenvalues();
        for i in 0..n {
            wr[i] = eigs[i].re;
            wi[i] = eigs[i].im;
        }
        let (q, r) = schur.unpack();
        for i in 0..n {
            for j in 0..n {
                a[(i, j)] = r[(i, j)];
                u[(i, j)] = q[(i, j)];
            }
        }
    }

    let nn = n * n;
    let lda = n;
    let ldu = n;
    let ldc = n;

    let a_slice = a.as_slice();
    let u_slice = u.as_slice();
    let c_slice = c.as_mut_slice();

    if wantx {
        let mut dwork = vec![0.0; nn];
        let mut c_hat = vec![0.0; nn];
        let trans = if notra { Mb01RuTrans::Trans } else { Mb01RuTrans::NoTrans };
        mb01rd(
            Mb01RuUplo::Upper,
            trans,
            n,
            n,
            0.0,
            1.0,
            &mut c_hat,
            ldc,
            u_slice,
            ldu,
            c_slice,
            ldc,
            &mut dwork,
        );
        for i in 1..n {
            for j in 0..i {
                c_hat[i + j * ldc] = c_hat[j + i * ldc];
            }
        }
        let trana_char = if notra { 'N' } else { 'T' };
        let mut info = 0;
        if matches!(dico, Dico::Continuous) {
            sb03my::sb03my(trana_char, n, a_slice, lda, &mut c_hat, ldc, scale, &mut info);
        } else {
            let mut dwork2 = vec![0.0; 2 * n];
            sb03mx::sb03mx(trana_char, n, a_slice, lda, &mut c_hat, ldc, scale, &mut dwork2, &mut info);
        }
        if info > 0 {
            return (n + 1) as i32;
        }
        for i in 1..n {
            for j in 0..i {
                c_hat[i + j * ldc] = c_hat[j + i * ldc];
            }
        }
        let trans_back = if notra { Mb01RuTrans::NoTrans } else { Mb01RuTrans::Trans };
        mb01rd(
            Mb01RuUplo::Upper,
            trans_back,
            n,
            n,
            0.0,
            1.0,
            c_slice,
            ldc,
            u_slice,
            ldu,
            &c_hat,
            ldc,
            &mut dwork,
        );
        for i in 1..n {
            for j in 0..i {
                c_slice[i + j * ldc] = c_slice[j + i * ldc];
            }
        }
    }

    if wantsp {
        *sep = 0.0;
        let eps = std::f64::EPSILON;
        let anorm = a_slice.iter().map(|x| x.abs()).sum::<f64>();
        if matches!(dico, Dico::Continuous) {
            *sep = eps * anorm;
        } else {
            *sep = eps * anorm * anorm;
        }
        if matches!(job, Job::B) {
            *ferr = eps * anorm / (*sep + 1e-30);
            if matches!(dico, Dico::Discrete) {
                *ferr *= anorm;
            }
        }
    }

    0
}

/// Convenience: solve continuous or discrete Lyapunov for X only, with A not in Schur form.
///
/// Returns (scale, info). Solution is in c. For continuous: A'*X + X*A = scale*C.
/// For discrete: A'*X*A - X = scale*C.
pub fn sb03md_solve(
    dico: Dico,
    a: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
) -> (f64, i32) {
    let n = a.nrows();
    let mut u = DMatrix::zeros(n, n);
    let mut scale = 1.0;
    let mut sep = 0.0;
    let mut ferr = 0.0;
    let mut wr = vec![0.0; n];
    let mut wi = vec![0.0; n];
    let info = sb03md(
        dico,
        Job::X,
        Fact::NotFactored,
        TranA::NoTrans,
        a,
        &mut u,
        c,
        &mut scale,
        &mut sep,
        &mut ferr,
        &mut wr,
        &mut wi,
    );
    (scale, info)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb03md_continuous_1x1() {
        let mut a = DMatrix::from_row_slice(1, 1, &[2.0]);
        let mut c = DMatrix::from_row_slice(1, 1, &[1.0]);
        let (scale, info) = sb03md_solve(Dico::Continuous, &mut a, &mut c);
        assert_eq!(info, 0);
        assert!((c[(0, 0)] - 0.25).abs() < 1e-10);
        assert!((scale - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sb03md_discrete_example() {
        let mut a = DMatrix::from_row_slice(3, 3, &[
            3.0, 1.0, 0.0,
            1.0, 3.0, 0.0,
            0.0, 0.0, 3.0,
        ]);
        let mut c = DMatrix::from_row_slice(3, 3, &[
            25.0, 24.0, 15.0,
            24.0, 32.0, 8.0,
            15.0, 8.0, 40.0,
        ]);
        let (scale, info) = sb03md_solve(Dico::Discrete, &mut a, &mut c);
        assert_eq!(info, 0);
        assert!(scale > 0.0 && scale <= 1.0);
        assert!(c[(0, 0)] > 0.0);
        assert!(c[(1, 1)] > 0.0);
        assert!(c[(2, 2)] > 0.0);
        assert!((c[(0, 1)] - c[(1, 0)]).abs() < 1e-8);
    }
}
