//! UD01DD — Read sparse real matrix from input (SLICOT).
//!
//! Zero-fills A then reads lines "i j A(i,j)" (1-based). INFO = 0 success, 1 warning (index out of range).

use nalgebra::DMatrix;

/// Reads a sparse matrix. All elements are set to zero; then each line must contain "i j value"
/// (1-based indices). Stops at first invalid or empty line.
///
/// # Returns
/// 0 = success; 1 = at least one (i,j) was out of range (warning); < 0 = invalid argument.
pub fn ud01dd<I, S>(m: usize, n: usize, mut lines: I, a: &mut DMatrix<f64>) -> i32
where
    I: Iterator<Item = S>,
    S: AsRef<str>,
{
    if a.nrows() < m || a.ncols() < n {
        return -4;
    }
    a.fill(0.0);
    let mut info = 0_i32;
    for line in lines {
        let line = line.as_ref().trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
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
        let val: f64 = match parts[2].parse() {
            Ok(x) => x,
            Err(_) => break,
        };
        if i < 1 || i > m || j < 1 || j > n {
            info = 1;
            continue;
        }
        a[(i - 1, j - 1)] = val;
    }
    info
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ud01dd_zero_fill() {
        let mut a = DMatrix::from_element(3, 3, 1.0);
        let lines: Vec<&str> = vec![];
        assert_eq!(ud01dd(3, 3, lines.into_iter(), &mut a), 0);
        assert_eq!(a[(0, 0)], 0.0);
        assert_eq!(a[(2, 2)], 0.0);
    }

    #[test]
    fn test_ud01dd_example_like() {
        let mut a = DMatrix::zeros(6, 5);
        let data = [
            "1   1   -1.1",
            "6   1    1.5",
            "2   2   -2.2",
            "6   2    2.5",
            "3   3   -3.3",
            "6   3    3.5",
            "4   4   -4.4",
            "6   4    4.5",
            "5   5   -5.5",
            "6   5    5.5",
        ];
        assert_eq!(ud01dd(6, 5, data.iter().map(|s| *s), &mut a), 0);
        assert!((a[(0, 0)] + 1.1).abs() < 1e-10);
        assert!((a[(5, 0)] - 1.5).abs() < 1e-10);
        assert!((a[(5, 4)] - 5.5).abs() < 1e-10);
    }

    #[test]
    fn test_ud01dd_warning() {
        let mut a = DMatrix::zeros(2, 2);
        let data = ["1 1 1.0", "10 10 2.0"];
        assert_eq!(ud01dd(2, 2, data.iter().map(|s| *s), &mut a), 1);
        assert_eq!(a[(0, 0)], 1.0);
    }
}
