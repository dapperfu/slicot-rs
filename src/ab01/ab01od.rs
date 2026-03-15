//! AB01OD — Staircase form for (A,B) with optional forward/backward stages (SLICOT AB01OD)
//!
//! Reduces (A,B) to upper staircase form. Not yet implemented.

use nalgebra::DMatrix;

/// Reduction stages: Forward only, Backward only, or All.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Stages {
    Forward,
    Backward,
    All,
}

/// Whether to accumulate transformation matrix.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobUV {
    No,
    Init,
}

/// Reduces (A,B) to staircase form. Stub: validates inputs and returns INFO=1.
///
/// # Returns
/// 0 if N=M=0; 1 = not implemented; < 0 = invalid argument index.
pub fn ab01od(
    _stages: Stages,
    _jobu: JobUV,
    _jobv: JobUV,
    n: usize,
    m: usize,
    _a: &mut DMatrix<f64>,
    _b: &mut DMatrix<f64>,
    _u: Option<&mut DMatrix<f64>>,
    _v: Option<&mut DMatrix<f64>>,
    ncont: &mut usize,
    _indcon: &mut usize,
    _kstair: &mut [i32],
    _tol: f64,
) -> i32 {
    if n == 0 && m == 0 {
        *ncont = 0;
        return 0;
    }
    if _a.nrows() != n || _a.ncols() != n {
        return -6;
    }
    if _b.nrows() != n || _b.ncols() != m {
        return -8;
    }
    1
}
