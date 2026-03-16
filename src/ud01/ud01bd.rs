//! UD01BD — Read matrix polynomial coefficients from input (SLICOT).
//!
//! Reads P(0), P(1), ..., P(DP) from lines: each coefficient matrix is preceded by a
//! text line, then MP rows of NP numbers (row-major). Returns INFO: 0 = success, < 0 = invalid argument.

use nalgebra::DMatrix;

/// Reads the coefficients of the matrix polynomial P(s) from an iterator of lines.
///
/// Each of the (DP+1) coefficient matrices must be preceded by one text line (e.g. "P0", "P1"),
/// then MP lines each containing NP whitespace-separated numbers.
///
/// # Arguments
/// - `mp`, `np`: rows and columns of each coefficient matrix; >= 1.
/// - `dp`: degree of the polynomial; >= 0.
/// - `lines`: iterator of input lines (without trailing newline).
/// - `p`: output; must hold (DP+1) matrices of size MP×NP. Stored as slice of matrices: p[k] = P(k).
///
/// # Returns
/// 0 on success; -1 invalid MP, -2 invalid NP, -3 invalid DP, -4 p length < DP+1, -5 insufficient lines.
pub fn ud01bd<I, S>(
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
        let _title = match lines.next() {
            Some(l) => l,
            None => return -5,
        };
        let mut row_data = vec![0.0_f64; mp * np];
        for i in 0..mp {
            let line = match lines.next() {
                Some(l) => l,
                None => return -5,
            };
            let line = line.as_ref().trim();
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.len() < np {
                return -5;
            }
            for j in 0..np {
                row_data[i * np + j] = match tokens[j].parse() {
                    Ok(x) => x,
                    Err(_) => return -5,
                };
            }
        }
        p[k] = DMatrix::from_row_slice(mp, np, &row_data);
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ud01bd_invalid_mp() {
        let mut p = vec![DMatrix::zeros(1, 1); 1];
        assert_eq!(ud01bd(0, 1, 0, std::iter::empty::<&str>(), &mut p), -1);
    }

    #[test]
    fn test_ud01bd_invalid_np() {
        let mut p = vec![DMatrix::zeros(1, 1); 1];
        assert_eq!(ud01bd(1, 0, 0, std::iter::empty::<&str>(), &mut p), -2);
    }

    #[test]
    fn test_ud01bd_insufficient_output() {
        let mut p = vec![DMatrix::zeros(1, 1)]; // only 1 matrix, need 2 for dp=1
        let lines = ["P0", "1.0", "P1", "2.0"].iter().map(|s| *s);
        assert_eq!(ud01bd(1, 1, 1, lines, &mut p), -4);
    }

    #[test]
    fn test_ud01bd_dp0_one_coeff() {
        let mut p = vec![DMatrix::zeros(2, 2); 1];
        let lines = ["P0", "1.0 0.0", "0.0 1.0"].iter().map(|s| *s);
        assert_eq!(ud01bd(2, 2, 0, lines, &mut p), 0);
        assert_eq!(p[0][(0, 0)], 1.0);
        assert_eq!(p[0][(1, 1)], 1.0);
    }

    #[test]
    fn test_ud01bd_example_like() {
        let mut p = vec![DMatrix::zeros(4, 3); 3];
        let data = [
            "P0",
            "1.0  0.0  0.0",
            "0.0  2.0  4.0",
            "0.0  4.0  8.0",
            "0.0  6.0  12.0",
            "P1",
            "0.0  1.0  2.0",
            "1.0  0.0  0.0",
            "2.0  0.0  0.0",
            "3.0  0.0  0.0",
            "P2",
            "1.0  0.0  0.0",
            "0.0  0.0  0.0",
            "0.0  0.0  0.0",
            "0.0  0.0  0.0",
        ];
        assert_eq!(ud01bd(4, 3, 2, data.iter().map(|s| *s), &mut p), 0);
        assert_eq!(p[0][(0, 0)], 1.0);
        assert_eq!(p[0][(3, 2)], 12.0);
        assert_eq!(p[1][(1, 0)], 1.0);
        assert_eq!(p[2][(0, 0)], 1.0);
    }
}
