//! AB13MD — Upper bound on structured singular value (SLICOT).
//!
//! Full SLICOT-equivalent API. Complex matrix Z, block structure (NBLOCK, ITYPE).
//! Returns an upper bound on mu (structured singular value). Special cases (Z=0,
//! single full block) exact; general case uses spectral norm (valid upper bound).

use nalgebra::DMatrix;
use num_complex::Complex64;

/// Whether X from previous call is supplied.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Fact {
    /// X contains information from previous call.
    FromPrevious,
    /// X does not contain information from previous call.
    New,
}

/// Full AB13MD: upper bound on structured singular value (complex matrix).
///
/// # Arguments
/// * `fact` - use previous X or not
/// * `n` - order of Z
/// * `z` - N×N complex matrix (column-major, ldz >= n)
/// * `ldz` - leading dimension of Z
/// * `m` - number of diagonal blocks (>= 1 for non-trivial)
/// * `nblock` - block sizes, length M
/// * `itype` - block types (1=real, 2=complex), length M
/// * `x` - workspace input/output, length M + MR - 1 (MR = number of real blocks)
/// * `bound` - output upper bound
/// * `d` - output diagonal D, length N
/// * `g` - output diagonal G, length N
/// * `iwork` - integer workspace
/// * `dwork` - real workspace
/// * `ldwork` - length of DWORK (or -1 for query)
/// * `zwork` - complex workspace
/// * `lzwork` - length of ZWORK (or -1 for query)
///
/// # Returns
/// INFO: 0 = success (quick return); 1 = not implemented; <0 = invalid argument.
pub fn ab13md(
    _fact: Fact,
    n: usize,
    z: &[Complex64],
    ldz: usize,
    m: usize,
    nblock: &[i32],
    itype: &[i32],
    _x: &mut [f64],
    bound: &mut f64,
    d: &mut [f64],
    g: &mut [f64],
    _iwork: &mut [i32],
    _dwork: &mut [f64],
    ldwork: i32,
    _zwork: &mut [Complex64],
    lzwork: i32,
) -> i32 {
    if n == 0 {
        *bound = 0.0;
        return 0;
    }
    if ldz < n {
        return -4;
    }
    if z.len() < ldz * n {
        return -3;
    }
    if m < 1 {
        return -6;
    }
    if nblock.len() < m || itype.len() < m {
        return -7;
    }
    let n_i = n as i32;
    let m_i = m as i32;
    let min_ldwork = 2 * n_i * n_i * m_i - n_i * n_i + 9 * m_i * m_i + n_i * m_i + 11 * n_i + 33 * m_i - 11;
    let min_lzwork = 6 * n_i * n_i * m_i + 12 * n_i * n_i + 6 * m_i + 6 * n_i - 3;
    if ldwork >= 0 && ldwork < min_ldwork {
        return -16;
    }
    if lzwork >= 0 && lzwork < min_lzwork {
        return -18;
    }

    // Validate block structure (SLICOT INFO 1,2,3,4)
    let mut nsum: i32 = 0;
    for i in 0..m {
        if nblock[i] < 1 {
            return 1;
        }
        if itype[i] == 1 && nblock[i] != 1 {
            return 3;
        }
        if itype[i] != 1 && itype[i] != 2 {
            return 4;
        }
        nsum += nblock[i];
    }
    if nsum != n_i {
        return 2;
    }

    // D = 1, G = 0 (output)
    for i in 0..n {
        d[i] = 1.0;
        g[i] = 0.0;
    }

    // Build Z as N×N complex matrix (column-major from z)
    let mut z_mat = DMatrix::<Complex64>::zeros(n, n);
    for j in 0..n {
        for i in 0..n {
            z_mat[(i, j)] = z[i + j * ldz];
        }
    }

    let znorm = z_mat.norm();
    if znorm == 0.0 {
        *bound = 0.0;
        return 0;
    }

    // Special case: NBLOCK(1) = N (single block)
    if nblock[0] == n_i {
        if itype[0] == 1 {
            // 1-by-1 real block
            *bound = z_mat[(0, 0)].re.abs();
            return 0;
        }
        // N×N complex block: bound = max singular value
        let svd = z_mat.clone().svd(false, false);
        *bound = if svd.singular_values.len() > 0 {
            svd.singular_values[0]
        } else {
            0.0
        };
        return 0;
    }

    // General case: spectral norm is a valid upper bound on mu
    let svd = z_mat.svd(false, false);
    *bound = if svd.singular_values.len() > 0 {
        svd.singular_values[0]
    } else {
        0.0
    };
    0
}

/// Compatibility wrapper: (n, m) -> INFO. Uses minimal workspace sizes.
#[inline]
pub fn ab13md_nm(n: usize, m: usize) -> i32 {
    if n == 0 {
        return 0;
    }
    if m < 1 {
        return 0; // trivial
    }
    let ldz = n;
    let z: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); ldz * n];
    let nblock: Vec<i32> = vec![1; m];
    let itype: Vec<i32> = vec![2; m];
    let mut x = vec![0.0_f64; m + 1];
    let mut bound = 0.0;
    let mut d = vec![0.0; n];
    let mut g = vec![0.0; n];
    let ldw = (2 * n * n * m + 9 * m * m + n * m + 11 * n + 33 * m).max(1);
    let lzw = (6 * n * n * m + 12 * n * n + 6 * m + 6 * n).max(1);
    let mut iwork = vec![0i32; (4 * m + n).max(1)];
    let mut dwork = vec![0.0; ldw];
    let mut zwork = vec![Complex64::new(0.0, 0.0); lzw];
    ab13md(
        Fact::New,
        n,
        &z,
        ldz,
        m,
        &nblock,
        &itype,
        &mut x,
        &mut bound,
        &mut d,
        &mut g,
        &mut iwork,
        &mut dwork,
        ldw as i32,
        &mut zwork,
        lzw as i32,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ab13md_trivial() {
        let mut bound = 0.0;
        let mut d = vec![0.0; 0];
        let mut g = vec![0.0; 0];
        assert_eq!(
            ab13md(
                Fact::New,
                0,
                &[],
                1,
                1,
                &[1],
                &[1],
                &mut [],
                &mut bound,
                &mut d,
                &mut g,
                &mut [],
                &mut [],
                1,
                &mut [],
                1,
            ),
            0
        );
        assert_eq!(bound, 0.0);
    }

    #[test]
    fn test_ab13md_n1_single_block() {
        let z = [Complex64::new(3.0, 4.0)]; // 1x1 complex, spectral norm = 5
        let mut bound = 0.0;
        let mut d = vec![0.0; 1];
        let mut g = vec![0.0; 1];
        let mut x = vec![0.0; 1];
        let mut iwork = vec![0i32; 4];
        let mut dwork = vec![0.0; 100];
        let mut zwork = vec![Complex64::new(0.0, 0.0); 50];
        assert_eq!(
            ab13md(
                Fact::New,
                1,
                &z,
                1,
                1,
                &[1],
                &[2],
                &mut x,
                &mut bound,
                &mut d,
                &mut g,
                &mut iwork,
                &mut dwork,
                100,
                &mut zwork,
                50,
            ),
            0
        );
        assert!((bound - 5.0).abs() < 1e-10);
        assert_eq!(d[0], 1.0);
        assert_eq!(g[0], 0.0);
    }

    #[test]
    fn test_ab13md_nm_returns_zero() {
        assert_eq!(ab13md_nm(1, 1), 0);
    }
}
