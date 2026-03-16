//! UD01MD — Print real matrix (SLICOT). Output to 7 significant figures, L columns per line (1..5).

use nalgebra::DMatrix;
use std::fmt::Write;

/// Prints an M×N real matrix. L is the number of columns per line (1..=5).
/// Writes to the given string (appends). Returns INFO: 0 = success, < 0 = invalid argument.
pub fn ud01md<W: Write>(
    m: usize,
    n: usize,
    l: usize,
    a: &DMatrix<f64>,
    text: &str,
    out: &mut W,
) -> i32 {
    if m < 1 {
        return -1;
    }
    if n < 1 {
        return -2;
    }
    if l < 1 || l > 5 {
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
                let x = a[(i, col_start + j)];
                let _ = write!(out, "{:14.7e} ", x);
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

/// Convenience: format matrix to a String.
pub fn ud01md_string(m: usize, n: usize, l: usize, a: &DMatrix<f64>, text: &str) -> Result<String, i32> {
    let mut s = String::new();
    let info = ud01md(m, n, l, a, text, &mut s);
    if info != 0 {
        return Err(info);
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ud01md_invalid_l() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let mut out = String::new();
        assert_eq!(ud01md(2, 2, 0, &a, "A", &mut out), -3);
        assert_eq!(ud01md(2, 2, 6, &a, "A", &mut out), -3);
    }

    #[test]
    fn test_ud01md_small() {
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let s = ud01md_string(2, 2, 2, &a, "Matrix A").unwrap();
        assert!(s.contains("Matrix A"));
        assert!(s.contains("1.0") && s.contains("4.0"));
    }

    #[test]
    fn test_ud01md_4x4() {
        let a = DMatrix::from_row_slice(4, 4, &[
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        ]);
        let s = ud01md_string(4, 4, 4, &a, "Matrix A").unwrap();
        assert!(s.contains("Matrix A"));
        assert!(s.contains("( 4X 4)"));
        assert!(s.len() > 100);
    }
}
