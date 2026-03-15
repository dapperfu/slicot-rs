//! MA02PD — Count zero rows and zero columns (SLICOT MA02PD)

use nalgebra::DMatrix;

/// Returns the number of zero rows and zero columns of A. Tolerance: exact zero.
pub fn ma02pd(a: &DMatrix<f64>, nzr: &mut usize, nzc: &mut usize) -> i32 {
    let m = a.nrows();
    let n = a.ncols();
    *nzr = 0;
    *nzc = 0;
    for i in 0..m {
        let mut row_zero = true;
        for j in 0..n {
            if a[(i, j)] != 0.0 {
                row_zero = false;
                break;
            }
        }
        if row_zero {
            *nzr += 1;
        }
    }
    for j in 0..n {
        let mut col_zero = true;
        for i in 0..m {
            if a[(i, j)] != 0.0 {
                col_zero = false;
                break;
            }
        }
        if col_zero {
            *nzc += 1;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02pd() {
        let a = DMatrix::from_row_slice(2, 2, &[0.0, 0.0, 1.0, 0.0]);
        let mut nzr = 0;
        let mut nzc = 0;
        assert_eq!(ma02pd(&a, &mut nzr, &mut nzc), 0);
        assert_eq!(nzr, 1);
        assert_eq!(nzc, 1);
    }
}
