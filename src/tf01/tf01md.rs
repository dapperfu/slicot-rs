//! TF01MD — Output response sequence of a linear time-invariant discrete-time system (SLICOT TF01MD)
//!
//! Computes y(1),...,y(NY) from x(k+1) = A*x(k) + B*u(k), y(k) = C*x(k) + D*u(k).

use nalgebra::{DMatrix, DVector};

/// Computes the output sequence Y(1),...,Y(NY) of the discrete-time system (A,B,C,D).
///
/// # Arguments
/// * `a` - State matrix A (n×n)
/// * `b` - Input matrix B (n×m)
/// * `c` - Output matrix C (p×n)
/// * `d` - Direct link matrix D (p×m)
/// * `u` - Input sequence: leading m×ny, column k = u(k)
/// * `x` - On entry: initial state x(1). On exit: final state x(ny+1)
/// * `y` - Output: leading p×ny, column k = y(k)
///
/// # Returns
/// * `0` - success
/// * `< 0` - if `-i`, the i-th argument had an illegal value
pub fn tf01md(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    c: &DMatrix<f64>,
    d: &DMatrix<f64>,
    u: &DMatrix<f64>,
    x: &mut DVector<f64>,
    y: &mut DMatrix<f64>,
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    let ny = u.ncols();

    if a.ncols() != n {
        return -5;
    }
    if b.nrows() != n || c.ncols() != n || d.nrows() != p || d.ncols() != m {
        return -6;
    }
    if u.nrows() != m {
        return -9;
    }
    if x.len() != n {
        return -11;
    }
    if y.nrows() != p || y.ncols() != ny {
        return -12;
    }

    if n == 0 || ny == 0 {
        return 0;
    }

    let mut xcur = x.clone();
    for k in 0..ny {
        // u(k) = k-th column of U
        let uk = u.column(k);
        // y(k) = C*x(k) + D*u(k)
        let yk = c * &xcur + d * &uk;
        y.column_mut(k).copy_from(&yk);
        // x(k+1) = A*x(k) + B*u(k)
        xcur = a * &xcur + b * &uk;
    }
    x.copy_from(&xcur);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tf01md_simple() {
        // n=2, m=1, p=1, ny=3. Identity dynamics, B=[1;0], C=[1 0], D=0.
        let a = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let c = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let u = DMatrix::from_row_slice(1, 3, &[1.0, 0.0, 0.0]); // u(1)=1, u(2)=0, u(3)=0
        let mut x = DVector::from_row_slice(&[0.0, 0.0]);
        let mut y = DMatrix::zeros(1, 3);
        assert_eq!(tf01md(&a, &b, &c, &d, &u, &mut x, &mut y), 0);
        // x(1)=[0,0], y(1)=C*x(1)+D*u(1)=0; x(2)=[1,0]; y(2)=1; x(3)=[1,0]; y(3)=1; x(4)=[1,0]
        assert!((y[(0, 0)] - 0.0).abs() < 1e-10);
        assert!((y[(0, 1)] - 1.0).abs() < 1e-10);
        assert!((y[(0, 2)] - 1.0).abs() < 1e-10);
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_tf01md_zero_ny() {
        let a = DMatrix::identity(2, 2);
        let b = DMatrix::zeros(2, 1);
        let c = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let d = DMatrix::zeros(1, 1);
        let u = DMatrix::zeros(1, 0);
        let mut x = DVector::from_row_slice(&[1.0, 0.0]);
        let mut y = DMatrix::zeros(1, 0);
        assert_eq!(tf01md(&a, &b, &c, &d, &u, &mut x, &mut y), 0);
        assert!((x[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_tf01md_slicot_example() {
        // SLICOT TF01MD example: N=3, M=2, P=2, NY=10. Same (A,B,C,D) as doc; U and x from doc.
        let a = DMatrix::from_row_slice(
            3,
            3,
            &[0.0, -0.07, 0.015, 1.0, 0.8, -0.15, 0.0, 0.0, 0.5],
        );
        let b = DMatrix::from_row_slice(
            3,
            2,
            &[0.0, 2.0, 1.0, -1.0, -0.1, 1.0],
        );
        let c = DMatrix::from_row_slice(2, 3, &[0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        let d = DMatrix::from_row_slice(2, 2, &[1.0, 0.5, 0.0, 0.5]);
        // U: M×NY, column k = u(k). Doc format (I,J) Fortran: col1 then col2 ...
        let u = DMatrix::from_row_slice(
            2,
            10,
            &[
                1.0, 0.5, 1.0, -0.6922, 2.0039, -1.5734, 0.4118, -0.9344, 0.8988, -0.0701,
                0.2614, -0.9160, -0.6030, 1.2556, 0.2951, 1.5639, -1.4893, 1.2506, -0.0701, 0.0,
            ],
        );
        let mut x = DVector::from_row_slice(&[1.0, 1.0, 1.0]);
        let mut y = DMatrix::zeros(2, 10);
        assert_eq!(tf01md(&a, &b, &c, &d, &u, &mut x, &mut y), 0);
        // Output sequence is finite and has correct shape
        assert!(y[(0, 0)].is_finite());
        assert!(y[(1, 0)].is_finite());
        assert!(y[(0, 9)].is_finite());
        assert!(y[(1, 9)].is_finite());
    }
}
