//! UD01ND — Print matrix polynomial (SLICOT). Prints P(0), P(1), ..., P(DP) to 7 significant figures.

use nalgebra::DMatrix;
use std::fmt::Write;

/// Prints the coefficient matrices of the matrix polynomial. L = columns per line (1..=5).
/// For each degree k = 0..=DP prints TEXT followed by (k) and the MP×NP matrix P(k).
pub fn ud01nd<W: Write>(
    mp: usize,
    np: usize,
    dp: usize,
    l: usize,
    p: &[DMatrix<f64>],
    text: &str,
    out: &mut W,
) -> i32 {
    if mp < 1 {
        return -1;
    }
    if np < 1 {
        return -2;
    }
    if dp == usize::MAX {
        return -3;
    }
    if l < 1 || l > 5 {
        return -4;
    }
    if p.len() < dp + 1 {
        return -7;
    }
    for k in 0..=dp {
        if p[k].nrows() < mp || p[k].ncols() < np {
            return -7;
        }
        let title = if text.trim().is_empty() {
            format!("  P( {} ) ( {}X {})", k, mp, np)
        } else {
            format!(" {} ( {} ) ( {}X {})", text.trim(), k, mp, np)
        };
        let _ = writeln!(out, "\n {}", title);
        let _ = write!(out, "\n            ");
        let n_show = np.min(l);
        for c in 1..=n_show {
            let _ = write!(out, "{:>14} ", c);
        }
        let _ = writeln!(out);
        let mut col_start = 0_usize;
        while col_start < np {
            let cols_this = (np - col_start).min(l);
            for i in 0..mp {
                let _ = write!(out, "{:3}   ", i + 1);
                for j in 0..cols_this {
                    let x = p[k][(i, col_start + j)];
                    let _ = write!(out, "{:14.7e} ", x);
                }
                let _ = writeln!(out);
            }
            col_start += cols_this;
            if col_start < np {
                let _ = write!(out, "\n            ");
                for c in 1..=((np - col_start).min(l)) {
                    let _ = write!(out, "{:>14} ", col_start + c);
                }
                let _ = writeln!(out);
            }
        }
    }
    0
}

/// Convenience: format matrix polynomial to a String.
pub fn ud01nd_string(
    mp: usize,
    np: usize,
    dp: usize,
    l: usize,
    p: &[DMatrix<f64>],
    text: &str,
) -> Result<String, i32> {
    let mut s = String::new();
    let info = ud01nd(mp, np, dp, l, p, text, &mut s);
    if info != 0 {
        return Err(info);
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ud01nd_invalid_l() {
        let p = vec![DMatrix::zeros(2, 2); 1];
        let mut out = String::new();
        assert_eq!(ud01nd(2, 2, 0, 0, &p, "P", &mut out), -4);
    }

    #[test]
    fn test_ud01nd_dp0() {
        let p = vec![DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0])];
        let s = ud01nd_string(2, 2, 0, 2, &p, " P").unwrap();
        assert!(s.contains("P") && s.contains("( 0 )"));
        assert!(s.contains("2X 2"));
    }
}
