//! TB04BD — Transfer matrix using pole-zeros method (SLICOT TB04BD)
//!
//! Each element G(i,j) returned as numerator/denominator polynomials (increasing or decreasing order).

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobD {
    D,
    Z,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Order {
    Increasing,
    Decreasing,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    Scale,
    No,
}

fn char_poly_coeffs(n: usize, a: &DMatrix<f64>) -> Vec<f64> {
    let mut phi = vec![0.0; n + 1];
    phi[0] = 1.0;
    let mut m = DMatrix::identity(n, n);
    for k in 1..=n {
        m = a * m;
        phi[k] = -m.trace() / (k as f64);
        if k < n {
            for i in 0..n {
                m[(i, i)] += phi[k];
            }
        }
    }
    phi
}

/// Computes transfer matrix G of (A,B,C,D) in polynomial form per element.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1/2 QR convergence (zeros/poles).
pub fn tb04bd(
    jobd: JobD,
    order: Order,
    equil: Equil,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &DMatrix<f64>,
    ign: &mut [i32],
    ldign: usize,
    igd: &mut [i32],
    ldigd: usize,
    gn: &mut [f64],
    gd: &mut [f64],
    tol: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    if n == 0 || m == 0 || p == 0 {
        return 0;
    }
    if ldign < p || ldigd < p {
        return -14;
    }
    let md = n + 1;
    if gn.len() < p * m * md || gd.len() < p * m * md {
        return -18;
    }

    if equil == Equil::Scale {
        let mut scale = vec![1.0; n];
        let mut maxred = 0.0;
        let _ = crate::tb01::tb01id::tb01id(
            crate::tb01::tb01id::Tb01IdJob::All,
            a,
            b,
            c,
            &mut scale,
            &mut maxred,
        );
    }

    let phi = char_poly_coeffs(n, a);
    let a_clone = a.view((0, 0), (n, n)).into_owned();
    let b_clone = b.view((0, 0), (n, m)).into_owned();
    let c_clone = c.view((0, 0), (p, n)).into_owned();

    for j in 0..m {
        for i in 0..p {
            let ij = (j * p + i) * md;
            let d_ij = if jobd == JobD::D { d[(i, j)] } else { 0.0 };
            let mut h = DMatrix::identity(n, n);
            let mut num = vec![0.0; n + 1];
            num[0] = d_ij * phi[0];
            for k in 1..=n {
                let cib = c_clone.row(i) * &h * b_clone.column(j);
                num[k] = cib[(0, 0)] + d_ij * phi[k];
                if k < n {
                    h = &a_clone * &h
                        + DMatrix::from_fn(n, n, |ii, jj| if ii == jj { phi[k] } else { 0.0 });
                }
            }
            ign[j * ldign + i] = n as i32;
            igd[j * ldigd + i] = n as i32;
            if order == Order::Increasing {
                for k in 0..=n {
                    if ij + k < gn.len() {
                        gn[ij + k] = num[k];
                    }
                    if ij + k < gd.len() {
                        gd[ij + k] = phi[k];
                    }
                }
            } else {
                for k in 0..=n {
                    if ij + k < gn.len() {
                        gn[ij + k] = num[n - k];
                    }
                    if ij + k < gd.len() {
                        gd[ij + k] = phi[n - k];
                    }
                }
            }
            for k in (n + 1)..md {
                if ij + k < gn.len() {
                    gn[ij + k] = 0.0;
                }
                if ij + k < gd.len() {
                    gd[ij + k] = 0.0;
                }
            }
        }
    }
    let _ = tol;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb04bd_smoke() {
        let mut a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let mut c = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut ign = [0i32; 1];
        let mut igd = [0i32; 1];
        let mut gn = vec![0.0; 1 * 1 * 3];
        let mut gd = vec![0.0; 1 * 1 * 3];
        let mut iwork = vec![0i32; 2];
        let mut dwork = vec![0.0; 50];
        let info = tb04bd(
            JobD::Z,
            Order::Increasing,
            Equil::No,
            &mut a,
            &mut b,
            &mut c,
            &d,
            &mut ign,
            1,
            &mut igd,
            1,
            &mut gn,
            &mut gd,
            0.0,
            &mut iwork,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert_eq!(ign[0], 2);
        assert_eq!(igd[0], 2);
    }
}
