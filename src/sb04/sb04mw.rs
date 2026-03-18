//! SB04MW — Solve linear system with upper Hessenberg matrix (compact storage) (SLICOT SB04MW).
//!
//! Used by SB04MY for the Sylvester equation back substitution. Gaussian elimination
//! with partial pivoting; matrix stored compactly row-wise.

const ZERO: f64 = 0.0;

/// Solves M×M system with upper Hessenberg matrix stored in D compactly row-wise,
/// and RHS in the next M elements. On exit, solution (permuted) overwrites the RHS part of D.
/// IPR(0..M) will contain solution component indices (which D index holds the i-th solution).
///
/// D layout: first M*(M+1)/2 + M elements = Hessenberg matrix compact, then M elements = RHS.
/// Total length >= M*(M+1)/2 + 2*M.
///
/// # Returns
/// 0 on success; 1 if singular.
pub fn sb04mw(m: usize, d: &mut [f64], ipr: &mut [i32]) -> i32 {
    if m == 0 {
        return 0;
    }
    let mind = (m * (m + 1)) / 2 + 2 * m;
    if d.len() < mind || ipr.len() < 2 * m {
        return -1;
    }

    // Fortran: first RHS at 1-based index (M*(M+3))/2; 0-based = that minus 1
    let rhs_start = (m * (m + 3)) / 2 - 1;
    let mut m1 = m;
    let mut i1 = 0usize; // start of row i in compact storage (0-based)

    for i in 0..m {
        let mpi = m + i;
        ipr[mpi] = i1 as i32;
        ipr[i] = (rhs_start + i) as i32;
        i1 += m1;
        if i > 0 {
            m1 -= 1;
        }
    }

    let m1_count = m - 1;

    for i in 0..m1_count {
        let i1 = i + 1;
        let mpi = m + i;
        let mpi1 = mpi + 1;
        let mut iprm_idx = ipr[mpi] as usize;
        let mut iprm1_idx = ipr[mpi1] as usize;
        let mut d1 = d[iprm_idx];
        let d2 = d[iprm1_idx];

        if d1.abs() <= d2.abs() {
            ipr[mpi] = iprm1_idx as i32;
            ipr[mpi1] = iprm_idx as i32;
            ipr.swap(i, i1);
            iprm_idx = ipr[mpi] as usize;
            iprm1_idx = ipr[mpi1] as usize;
            d1 = d2;
        }

        if d1 == ZERO {
            return 1;
        }

        let mult = -d[iprm1_idx] / d1;
        let ipr_i1 = ipr[i1] as usize;
        let ipr_i = ipr[i] as usize;
        d[ipr_i1] += mult * d[ipr_i];

        let len = m - i - 1;
        for k in 0..len {
            d[iprm1_idx + 1 + k] += mult * d[iprm_idx + 1 + k];
        }
        ipr[mpi1] = (iprm1_idx + 1) as i32;
    }

    let m2 = 2 * m - 1;
    let ipr_m2 = ipr[m2] as usize;
    if d[ipr_m2] == ZERO {
        return 1;
    }

    let ipr_m = ipr[m - 1] as usize;
    d[ipr_m] /= d[ipr_m2];

    let mut mpi = m2;
    for ii in 0..m1_count {
        let i = m1_count - 1 - ii;
        mpi -= 1;
        let iprm = ipr[mpi] as usize;
        let mut mult = ZERO;
        let mut iprm1 = iprm;
        for k in (i + 1)..m {
            iprm1 += 1;
            mult += d[ipr[k] as usize] * d[iprm1];
        }
        let ipr_i = ipr[i] as usize;
        d[ipr_i] = (d[ipr_i] - mult) / d[iprm];
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sb04mw_m1() {
        let m = 1;
        // D length >= M*(M+1)/2 + 2*M = 3. Matrix [0], RHS at index 1.
        let mut d = vec![1.0, 6.0, 0.0];
        let mut ipr = vec![0i32; 2 * m];
        assert_eq!(sb04mw(m, &mut d, &mut ipr), 0);
        assert!((d[ipr[0] as usize] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_sb04mw_m2() {
        let m = 2;
        // Compact: row0 at 0,1 (2 elts), row1 at 2,3 (2 elts); RHS at 4,5. len >= 7.
        let mut d = vec![
            1.0, 2.0,   // row 0: (0,0), (0,1)
            0.0, 3.0,   // row 1: (1,0), (1,1) subdiag then diag
            1.0, 2.0,   // RHS
            0.0,
        ];
        let mut ipr = vec![0i32; 2 * m];
        assert_eq!(sb04mw(m, &mut d, &mut ipr), 0);
        // 1*x0 + 2*x1 = 1, 0*x0 + 3*x1 = 2 => x1 = 2/3, x0 = -1/3
        let x0 = d[ipr[0] as usize];
        let x1 = d[ipr[1] as usize];
        assert!((x1 - 2.0 / 3.0).abs() < 1e-10);
        assert!((x0 - (-1.0 / 3.0)).abs() < 1e-10);
    }
}
