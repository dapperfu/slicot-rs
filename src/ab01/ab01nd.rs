//! AB01ND — Find controllable realization for multi-input system (SLICOT AB01ND)
//!
//! Reduces (A,B) to orthogonal canonical form (block Hessenberg). Not yet implemented.

use nalgebra::DMatrix;

/// Whether to accumulate the orthogonal transformation matrix Z.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobZ {
    No,
    Factored,
    Init,
}

/// Finds a controllable realization for the multi-input system (A,B). Stub: validates inputs and returns INFO=1.
///
/// # Returns
/// 0 if N=M=0; 1 = not implemented; < 0 = invalid argument index.
pub fn ab01nd(
    _jobz: JobZ,
    n: usize,
    m: usize,
    _a: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    ncont: &mut usize,
    _indcon: &mut usize,
    _nblk: &mut [i32],
    _z: Option<&mut DMatrix<f64>>,
    _tol: f64,
) -> i32 {
    if n == 0 && m == 0 {
        *ncont = 0;
        return 0;
    }
    if _a.nrows() != n || _a.ncols() != n {
        return -4;
    }
    if _b.nrows() != n || _b.ncols() != m {
        return -6;
    }
    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_ab01nd_n0_m0() {
        let mut a = DMatrix::zeros(0, 0);
        let mut b = DMatrix::zeros(0, 0);
        let mut ncont = 1;
        let mut indcon = 0;
        let mut nblk = [0i32; 1];
        assert_eq!(
            ab01nd(JobZ::No, 0, 0, &mut a, &mut b, &mut ncont, &mut indcon, &mut nblk, None, 0.0),
            0
        );
        assert_eq!(ncont, 0);
    }
}
