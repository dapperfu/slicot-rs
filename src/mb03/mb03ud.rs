//! MB03UD — SVD of real upper triangular matrix.
//!
//! A = Q*S*P' with S diagonal (singular values). Uses nalgebra SVD.
//! JOBQ='V' => compute Q; JOBP='V' => compute P' (returned in A).

use nalgebra::DMatrix;
use std::cmp::Ordering;

/// Job for left singular vectors.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb03udJobQ {
    /// Compute Q.
    Compute,
    /// Do not compute Q.
    No,
}

/// Job for right singular vectors.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mb03udJobP {
    /// Compute P' (output in A).
    Compute,
    /// Do not compute P.
    No,
}

/// SVD of upper triangular A. A is N×N (column-major LDA).
/// On exit: SV(1:N) = singular values descending; if JOBP=Compute, A contains P'; if JOBQ=Compute, Q is in Q (LDQ×N).
/// DWORK length at least 5*N; on exit DWORK(1) = optimal LDWORK.
/// Returns INFO: 0 = success, <0 = invalid argument, >0 = DBDSQR did not converge (not used with nalgebra).
pub fn mb03ud(
    jobq: Mb03udJobQ,
    jobp: Mb03udJobP,
    n: usize,
    a: &mut [f64],
    lda: usize,
    q: &mut [f64],
    ldq: usize,
    sv: &mut [f64],
    dwork: &mut [f64],
) -> i32 {
    if n == 0 {
        if !dwork.is_empty() {
            dwork[0] = 1.0;
        }
        return 0;
    }
    if n > a.len() / lda.max(1) || n > sv.len() {
        return -3;
    }
    let wantq = jobq == Mb03udJobQ::Compute;
    let wantp = jobp == Mb03udJobP::Compute;
    if wantq && (ldq < n || q.len() < n * ldq) {
        return -7;
    }
    if dwork.len() < 5 * n {
        return -10;
    }

    let mut mat = DMatrix::from_fn(n, n, |i, j| {
        if i <= j {
            a[i + j * lda]
        } else {
            0.0
        }
    });

    let svd = mat.svd(true, wantp);
    let mut s = svd.singular_values.iter().copied().collect::<Vec<_>>();
    s.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));
    for (i, &v) in s.iter().enumerate() {
        sv[i] = v;
    }

    if wantq {
        let u = svd.u.unwrap();
        for j in 0..n {
            for i in 0..n {
                q[i + j * ldq] = u[(i, j)];
            }
        }
    }

    if wantp {
        let v_t = svd.v_t.unwrap();
        for i in 0..n {
            for j in 0..n {
                a[i + j * lda] = v_t[(i, j)];
            }
        }
    }

    dwork[0] = (5 * n) as f64;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb03ud_1x1() {
        let mut a = vec![2.0];
        let mut q = vec![0.0];
        let mut sv = vec![0.0];
        let mut dwork = vec![0.0; 10];
        let info = mb03ud(
            Mb03udJobQ::No,
            Mb03udJobP::No,
            1,
            &mut a,
            1,
            &mut q,
            1,
            &mut sv,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert!((sv[0] - 2.0).abs() < 1e-10);
    }
}
