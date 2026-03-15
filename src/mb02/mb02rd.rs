//! MB02RD — Solve H*X = B with upper Hessenberg H factored by MB02SD (SLICOT).
//!
//! Uses factored H = P*L*U and IPIV from MB02SD; overwrites B with X.

use nalgebra::DMatrix;

/// Solves H*X = B. H is the N×N factored matrix from MB02SD (L in lower, U in upper).
/// ipiv is the 1-based pivot vector from MB02SD. B is N×NRHS, overwritten with X.
/// Returns 0 on success, or <0 if invalid input.
pub fn mb02rd(
    n: usize,
    h: &DMatrix<f64>,
    ipiv: &[i32],
    b: &mut DMatrix<f64>,
) -> i32 {
    if n == 0 {
        return 0;
    }
    if h.nrows() != n || h.ncols() != n || ipiv.len() < n || b.nrows() != n {
        return -1;
    }
    let nrhs = b.ncols();

    // Apply P to B: row i of P*B = row (ipiv(i)-1) of B
    let mut work: Vec<f64> = (0..n * nrhs).map(|_| 0.0).collect();
    for j in 0..nrhs {
        for i in 0..n {
            let src = (ipiv[i] as usize) - 1;
            work[i + j * n] = b[(src, j)];
        }
    }
    for j in 0..nrhs {
        for i in 0..n {
            b[(i, j)] = work[i + j * n];
        }
    }

    // Solve L*y = P*B (forward). L is unit lower bidiagonal: L(i,i)=1, L(i+1,i)=h(i+1,i)
    for j in 0..nrhs {
        for i in 1..n {
            let mult = h[(i, i - 1)];
            b[(i, j)] -= mult * b[(i - 1, j)];
        }
    }

    // Solve U*X = y (back substitution)
    for j in 0..nrhs {
        for i in (0..n).rev() {
            let mut t = b[(i, j)];
            for k in (i + 1)..n {
                t -= h[(i, k)] * b[(k, j)];
            }
            let uii = h[(i, i)];
            if uii == 0.0 {
                return 1;
            }
            b[(i, j)] = t / uii;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mb02::mb02sd::mb02sd_matrix;

    #[test]
    fn test_mb02rd_trivial() {
        let a = DMatrix::<f64>::zeros(0, 0);
        let mut b = DMatrix::<f64>::zeros(0, 0);
        let ipiv: Vec<i32> = vec![];
        assert_eq!(mb02rd(0, &a, &ipiv, &mut b), 0);
    }

    #[test]
    fn test_mb02rd_solve() {
        let orig = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 1.0, 3.0]);
        let mut h = orig.clone();
        let mut ipiv = vec![0i32; 2];
        assert_eq!(mb02sd_matrix(&mut h, &mut ipiv), 0);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 2.0]);
        assert_eq!(mb02rd(2, &h, &ipiv, &mut b), 0);
        let x = b[(0, 0)];
        let y = b[(1, 0)];
        let r0 = orig[(0, 0)] * x + orig[(0, 1)] * y - 1.0;
        let r1 = orig[(1, 0)] * x + orig[(1, 1)] * y - 2.0;
        assert!(r0.abs() < 1e-10 && r1.abs() < 1e-10, "residuals {} {}", r0, r1);
    }
}
