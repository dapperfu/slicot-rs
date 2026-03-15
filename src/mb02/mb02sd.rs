//! MB02SD — LU factorization of an upper Hessenberg matrix (SLICOT).
//!
//! H = P*L*U with partial pivoting; L unit lower bidiagonal, U upper triangular.
//! H is overwritten: L (without unit diagonal) in lower part, U in upper.

use nalgebra::DMatrix;

/// Computes LU factorization of upper Hessenberg H (column-major, LDH×N).
/// H is overwritten with L and U; ipiv[i] = 1-based index of row swapped with row i.
/// Returns 0, or >0 if U(i,i)=0 (singular).
pub fn mb02sd(n: usize, h: &mut [f64], ldh: usize, ipiv: &mut [i32]) -> i32 {
    if n == 0 {
        return 0;
    }
    if ldh < n || h.len() < n * ldh || ipiv.len() < n {
        return -1;
    }
    for i in 0..n {
        ipiv[i] = (i + 1) as i32;
    }
    for k in 0..n {
        let mut imax = k;
        let mut amax = h[k + k * ldh].abs();
        if k + 1 < n {
            let a1 = h[k + 1 + k * ldh].abs();
            if a1 > amax {
                imax = k + 1;
                amax = a1;
            }
        }
        if amax == 0.0 {
            return (k + 1) as i32;
        }
        if imax != k {
            ipiv[k] = (imax + 1) as i32;
            for j in 0..n {
                h.swap(k + j * ldh, imax + j * ldh);
            }
        }
        if k + 1 < n {
            let mult = h[k + 1 + k * ldh] / h[k + k * ldh];
            h[k + 1 + k * ldh] = mult;
            for j in (k + 1)..n {
                h[k + 1 + j * ldh] -= mult * h[k + j * ldh];
            }
        }
    }
    0
}

/// DMatrix wrapper: overwrites H with L/U and fills ipiv (1-based). Returns 0 or >0 if singular.
pub fn mb02sd_matrix(h: &mut DMatrix<f64>, ipiv: &mut [i32]) -> i32 {
    let n = h.nrows();
    if n != h.ncols() || ipiv.len() < n {
        return -1;
    }
    let ldh = h.nrows();
    let data = h.as_mut_slice();
    mb02sd(n, data, ldh, ipiv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mb02sd_trivial() {
        let mut h = vec![0.0];
        let mut ipiv = vec![0];
        assert_eq!(mb02sd(0, &mut h, 0, &mut ipiv), 0);
    }

    #[test]
    fn test_mb02sd_2x2() {
        let mut h = vec![1.0, 2.0, 1.0, 3.0];
        let mut ipiv = vec![0, 0];
        assert_eq!(mb02sd(2, &mut h, 2, &mut ipiv), 0);
        assert!(h[0] != 0.0 && h[3] != 0.0);
    }
}
