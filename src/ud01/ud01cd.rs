//! UD01CD — Read sparse matrix polynomial from input (SLICOT).
//!
//! Zero-fills P then reads nonzero (i,j,d) entries: one line (i, j, d), next line d+1 coefficients.
//! INFO = 0 success, 1 warning (index out of range), < 0 invalid argument.

use nalgebra::DMatrix;

/// Reads a sparse matrix polynomial. All elements are first set to zero; then each nonzero
/// (i,j)-element is given by one line "i j d" (1-based) and the next line with d+1 coefficients
/// for s^0, s^1, ..., s^d.
///
/// # Returns
/// 0 = success; 1 = at least one (i,j,d) was out of range (warning); < 0 = invalid argument.
pub fn ud01cd<I, S>(
    mp: usize,
    np: usize,
    dp: usize,
    mut lines: I,
    p: &mut [DMatrix<f64>],
) -> i32
where
    I: Iterator<Item = S>,
    S: AsRef<str>,
{
    if mp < 1 {
        return -1;
    }
    if np < 1 {
        return -2;
    }
    if dp == usize::MAX {
        return -3;
    }
    if p.len() < dp + 1 {
        return -4;
    }
    for k in 0..=dp {
        p[k] = DMatrix::zeros(mp, np);
    }
    let mut info = 0_i32;
    while let Some(line1) = lines.next() {
        let line1 = line1.as_ref().trim();
        if line1.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line1.split_whitespace().collect();
        if parts.len() < 3 {
            break;
        }
        let i: usize = match parts[0].parse() {
            Ok(x) => x,
            Err(_) => break,
        };
        let j: usize = match parts[1].parse() {
            Ok(x) => x,
            Err(_) => break,
        };
        let d: usize = match parts[2].parse() {
            Ok(x) => x,
            Err(_) => break,
        };
        if i < 1 || i > mp || j < 1 || j > np || d > dp {
            info = 1;
            if let Some(coeff_line) = lines.next() {
                let _ = coeff_line;
            }
            continue;
        }
        let coeff_line = match lines.next() {
            Some(l) => l,
            None => break,
        };
        let coeff_line = coeff_line.as_ref().trim();
        let coeffs: Vec<f64> = coeff_line
            .split_whitespace()
            .filter_map(|w| w.parse().ok())
            .collect();
        if coeffs.len() < d + 1 {
            break;
        }
        let row = i - 1;
        let col = j - 1;
        for k in 0..=d {
            p[k][(row, col)] = coeffs[k];
        }
    }
    info
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ud01cd_invalid_mp() {
        let mut p = vec![DMatrix::zeros(1, 1); 1];
        assert_eq!(ud01cd(0, 1, 0, std::iter::empty::<&str>(), &mut p), -1);
    }

    #[test]
    fn test_ud01cd_zero_fill() {
        let mut p = vec![DMatrix::from_element(2, 2, 1.0); 2];
        let lines: Vec<&str> = vec![];
        assert_eq!(ud01cd(2, 2, 1, lines.into_iter(), &mut p), 0);
        assert_eq!(p[0][(0, 0)], 0.0);
        assert_eq!(p[1][(1, 1)], 0.0);
    }

    #[test]
    fn test_ud01cd_one_entry() {
        let mut p = vec![DMatrix::zeros(4, 3); 3];
        let data = ["1  1  1", "1.0  1.0"];
        assert_eq!(ud01cd(4, 3, 2, data.iter().map(|s| *s), &mut p), 0);
        assert_eq!(p[0][(0, 0)], 1.0);
        assert_eq!(p[1][(0, 0)], 1.0);
    }

    #[test]
    fn test_ud01cd_example_like() {
        let mut p = vec![DMatrix::zeros(4, 3); 3];
        let data = [
            "1  1  1",
            "1.0  1.0",
            "2  2  2",
            "2.0  0.0  1.0",
            "3  3  2",
            "0.0  3.0  1.0",
            "4  1  0",
            "4.0",
        ];
        assert_eq!(ud01cd(4, 3, 2, data.iter().map(|s| *s), &mut p), 0);
        assert_eq!(p[0][(0, 0)], 1.0);
        assert_eq!(p[1][(0, 0)], 1.0);
        assert_eq!(p[0][(1, 1)], 2.0);
        assert_eq!(p[2][(1, 1)], 1.0);
        assert_eq!(p[1][(2, 2)], 3.0);
        assert_eq!(p[2][(2, 2)], 1.0);
        assert_eq!(p[0][(3, 0)], 4.0);
    }

    #[test]
    fn test_ud01cd_warning_out_of_range() {
        let mut p = vec![DMatrix::zeros(2, 2); 1];
        let data = ["10  10  0", "1.0"];
        assert_eq!(ud01cd(2, 2, 0, data.iter().map(|s| *s), &mut p), 1);
    }
}
