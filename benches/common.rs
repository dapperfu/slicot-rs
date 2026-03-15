//! Shared size ladders and matrix/vector builders for SLICOT benchmarks.
//! No tiny sizes (e.g. 2×2); dimensions scale so speed differences are visible.

use nalgebra::{DMatrix, DVector};

/// Primary dimension ladder for state size / matrix rows (no tiny sizes).
pub const SIZE_LADDER_N: &[usize] = &[32, 64, 128, 256, 512, 1024];

/// Power-of-two ladder for FFT/signal routines (DE01OD, DG01*, etc.).
pub const SIZE_LADDER_POW2: &[usize] = &[64, 128, 256, 512, 1024, 2048];

/// Build an n×n matrix with deterministic values for benchmarking.
#[inline]
pub fn matrix_nn(n: usize) -> DMatrix<f64> {
    DMatrix::from_fn(n, n, |i, j| (i + j) as f64 * 0.1)
}

/// Build an n×m matrix with deterministic values.
#[inline]
pub fn matrix_nm(n: usize, m: usize) -> DMatrix<f64> {
    DMatrix::from_fn(n, m, |i, j| (i * 2 + j) as f64 * 0.1)
}

/// Build a p×n matrix (e.g. C in state-space).
#[inline]
pub fn matrix_pn(p: usize, n: usize) -> DMatrix<f64> {
    DMatrix::from_fn(p, n, |i, j| (i + j * 2) as f64 * 0.1)
}

/// Build a p×m matrix (e.g. D in state-space).
#[inline]
pub fn matrix_pm(p: usize, m: usize) -> DMatrix<f64> {
    DMatrix::from_fn(p, m, |i, j| (i + j) as f64 * 0.1)
}

/// Build a vector of length n with deterministic values.
#[inline]
pub fn vector_n(n: usize) -> DVector<f64> {
    DVector::from_fn(n, |i, _| (i + 1) as f64 * 0.1)
}

/// Build a slice-backed vector of length n (for routines that take &mut [f64]).
#[inline]
pub fn vec_n(n: usize) -> Vec<f64> {
    (0..n).map(|i| (i + 1) as f64 * 0.1).collect()
}

/// State-space (A,B,C,D) for dimensions n, m, p. A n×n, B n×m, C p×n, D p×m.
pub fn state_space_matrices(
    n: usize,
    m: usize,
    p: usize,
) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    let a = matrix_nn(n);
    let b = matrix_nm(n, m);
    let c = matrix_pn(p, n);
    let d = matrix_pm(p, m);
    (a, b, c, d)
}

/// Derive m from n (e.g. m = n/2) for (n, m) routines.
#[inline]
pub fn m_from_n(n: usize) -> usize {
    (n / 2).max(1)
}

/// Derive p from n for (n, m, p) routines.
#[inline]
pub fn p_from_n(n: usize) -> usize {
    (n / 2).max(1)
}
