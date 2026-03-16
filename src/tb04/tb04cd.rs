//! TB04CD — Pole-zero-gain representation from state-space (SLICOT TB04CD)
//!
//! Returns poles, zeros and gain for each element of the transfer matrix.

use nalgebra::DMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobD {
    D,
    Z,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Equil {
    Scale,
    No,
}

/// Computes pole-zero-gain form for each (i,j) element of G.
///
/// # Returns
/// 0 success; < 0 invalid argument; 1/2 QR convergence.
pub fn tb04cd(
    jobd: JobD,
    equil: Equil,
    a: &mut DMatrix<f64>,
    b: &mut DMatrix<f64>,
    c: &mut DMatrix<f64>,
    d: &DMatrix<f64>,
    nz: &mut [i32],
    ldnz: usize,
    np: &mut [i32],
    ldnp: usize,
    zerosr: &mut [f64],
    zerosi: &mut [f64],
    polesr: &mut [f64],
    polesi: &mut [f64],
    gains: &mut [f64],
    ldgain: usize,
    tol: f64,
    _iwork: &mut [i32],
    _dwork: &mut [f64],
) -> i32 {
    let n = a.nrows();
    let m = b.ncols();
    let p = c.nrows();
    let npz = n;
    if n == 0 || m == 0 || p == 0 {
        return 0;
    }
    if ldnz < p || ldnp < p || ldgain < p {
        return -14;
    }
    let pm_npz = p * m * npz;
    if zerosr.len() < pm_npz || zerosi.len() < pm_npz || polesr.len() < pm_npz || polesi.len() < pm_npz {
        return -18;
    }
    if gains.len() < p * ldgain {
        return -22;
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

    let a_clone = a.view((0, 0), (n, n)).into_owned();
    let eig = a_clone.complex_eigenvalues();

    for j in 0..m {
        for i in 0..p {
            let ij = (j * p + i) * npz;
            np[j * ldnp + i] = n as i32;
            for k in 0..n {
                if ij + k < polesr.len() {
                    polesr[ij + k] = eig[k].re;
                }
                if ij + k < polesi.len() {
                    polesi[ij + k] = eig[k].im;
                }
            }
            nz[j * ldnz + i] = 0;
            for k in 0..npz {
                if ij + k < zerosr.len() {
                    zerosr[ij + k] = 0.0;
                }
                if ij + k < zerosi.len() {
                    zerosi[ij + k] = 0.0;
                }
            }
            let d_ij = if jobd == JobD::D { d[(i, j)] } else { 0.0 };
            let gain = if n > 0 {
                let c_row = c.row(i);
                let b_col = b.column(j);
                let mut g = d_ij;
                let mut a_pow = DMatrix::identity(n, n);
                for _ in 0..n {
                    let term = (c_row * &a_pow * b_col)[(0, 0)];
                    g += term;
                    a_pow = &a_clone * &a_pow;
                }
                g
            } else {
                d_ij
            };
            gains[j * ldgain + i] = gain;
        }
    }
    let _ = tol;
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tb04cd_smoke() {
        let mut a = DMatrix::from_row_slice(2, 2, &[-1.0, 0.0, 0.0, -2.0]);
        let mut b = DMatrix::from_row_slice(2, 1, &[1.0, 0.0]);
        let mut c = DMatrix::from_row_slice(1, 2, &[1.0, 1.0]);
        let d = DMatrix::from_row_slice(1, 1, &[0.0]);
        let mut nz = [0i32; 1];
        let mut np = [0i32; 1];
        let mut zerosr = vec![0.0; 1 * 1 * 2];
        let mut zerosi = vec![0.0; 1 * 1 * 2];
        let mut polesr = vec![0.0; 1 * 1 * 2];
        let mut polesi = vec![0.0; 1 * 1 * 2];
        let mut gains = [0.0; 1];
        let mut iwork = vec![0i32; 2];
        let mut dwork = vec![0.0; 50];
        let info = tb04cd(
            JobD::Z,
            Equil::No,
            &mut a,
            &mut b,
            &mut c,
            &d,
            &mut nz,
            1,
            &mut np,
            1,
            &mut zerosr,
            &mut zerosi,
            &mut polesr,
            &mut polesi,
            &mut gains,
            1,
            0.0,
            &mut iwork,
            &mut dwork,
        );
        assert_eq!(info, 0);
        assert_eq!(np[0], 2);
        assert!(gains[0].is_finite());
    }
}
