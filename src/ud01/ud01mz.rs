//! UD01MZ — Print complex matrix (SLICOT). Output to 7 significant figures, L columns per line (1..3).

use nalgebra::DMatrix;
use num_complex::Complex64;
use std::fmt::Write;

/// Prints an M×N complex matrix. L is the number of columns per line (1..=3).
/// Writes to the given string (appends). Returns INFO: 0 = success, < 0 = invalid argument.
pub fn ud01mz<W: Write>(
    m: usize,
    n: usize,
    l: usize,
    a: &DMatrix<Complex64>,
    text: &str,
    out: &mut W,
) -> i32 {
    if m < 1 {
        return -1;
    }
    if n < 1 {
        return -2;
    }
    if l < 1 || l > 3 {
        return -3;
    }
    if a.nrows() < m || a.ncols() < n {
        return -5;
    }
    let _ = writeln!(out, " {} ( {}X {})", text.trim(), m, n);
    let _ = write!(out, "\n            ");
    for c in 1..=n.min(l) {
        let _ = write!(out, "{:>14} ", c);
    }
    let _ = writeln!(out);
    let mut col_start = 0_usize;
    while col_start < n {
        let cols_this = (n - col_start).min(l);
        for i in 0..m {
            let _ = write!(out, "{:3}   ", i + 1);
            for j in 0..cols_this {
                let z = a[(i, col_start + j)];
                let _ = write!(out, "({:14.7e},{:14.7e}) ", z.re, z.im);
            }
            let _ = writeln!(out);
        }
        col_start += cols_this;
        if col_start < n {
            let _ = write!(out, "\n            ");
            for c in 1..=((n - col_start).min(l)) {
                let _ = write!(out, "{:>14} ", col_start + c);
            }
            let _ = writeln!(out);
        }
    }
    0
}

/// Convenience: format complex matrix to a String.
pub fn ud01mz_string(
    m: usize,
    n: usize,
    l: usize,
    a: &DMatrix<Complex64>,
    text: &str,
) -> Result<String, i32> {
    let mut s = String::new();
    let info = ud01mz(m, n, l, a, text, &mut s);
    if info != 0 {
        return Err(info);
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ud01mz_invalid_l() {
        let a = DMatrix::from_fn(2, 2, |_, _| Complex64::new(0.0, 0.0));
        let mut out = String::new();
        assert_eq!(ud01mz(2, 2, 0, &a, "A", &mut out), -3);
        assert_eq!(ud01mz(2, 2, 4, &a, "A", &mut out), -3);
    }

    #[test]
    fn test_ud01mz_small() {
        let a = DMatrix::from_fn(2, 2, |i, j| Complex64::new((i + j) as f64, (i as f64) * 0.5));
        let s = ud01mz_string(2, 2, 2, &a, "Matrix A").unwrap();
        assert!(s.contains("Matrix A"));
        assert!(s.contains("( 2X 2)"));
    }
}
