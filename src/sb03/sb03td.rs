//! SB03TD — SLICOT Lyapunov. Thin wrapper: discrete solve A'*X*A - X = C.
use nalgebra::DMatrix;

use super::sb03md::{sb03md_solve, Dico};

pub fn sb03td(n: usize, a: &DMatrix<f64>, x: &mut DMatrix<f64>) -> i32 {
    if n == 0 {
        return 0;
    }
    let mut a_mat = a.clone();
    let mut c_mat = x.clone();
    let (_scale, info) = sb03md_solve(Dico::Discrete, &mut a_mat, &mut c_mat);
    if info == 0 {
        for i in 0..n {
            for j in 0..n {
                x[(i, j)] = c_mat[(i, j)];
            }
        }
    }
    info
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sb03td() {
        let a = DMatrix::zeros(1, 1);
        let mut x = DMatrix::zeros(1, 1);
        assert_eq!(sb03td(1, &a, &mut x), 0);
    }
}
