//! DG01MD — Discrete Fourier transform or inverse of a complex signal (SLICOT DG01MD)
//
// INDI: 'D' = direct FFT, 'I' = inverse FFT. N must be a power of 2, N >= 2.
// XR, XI overwritten with result (column-major complex).

use std::f64::consts::PI;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dg01MdIndi {
    Direct,
    Inverse,
}

fn is_power_of_two(n: usize) -> bool {
    n >= 2 && (n & (n - 1)) == 0
}

/// In-place FFT or inverse FFT. xr and xi must have length >= N.
pub fn dg01md(indi: Dg01MdIndi, n: usize, xr: &mut [f64], xi: &mut [f64]) -> i32 {
    if !matches!(indi, Dg01MdIndi::Direct | Dg01MdIndi::Inverse) {
        return -1;
    }
    if !is_power_of_two(n) {
        return -2;
    }
    if xr.len() < n || xi.len() < n {
        return -3;
    }

    let pi2 = if indi == Dg01MdIndi::Direct { -2.0 * PI } else { 2.0 * PI };

    // Bit-reversal permutation
    let mut j = 1_usize;
    for i in 1..=n {
        if i < j {
            xr.swap(i - 1, j - 1);
            xi.swap(i - 1, j - 1);
        }
        let mut k = n / 2;
        while k >= 1 && j > k {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    // Radix-2 decimation-in-time
    let mut i = 1_usize;
    while i < n {
        let l = 2 * i;
        let whelp = pi2 / (l as f64);
        let wstpi = whelp.sin();
        let wstpr = -2.0 * (whelp / 2.0).sin().powi(2);
        let mut wr = 1.0_f64;
        let mut wi = 0.0_f64;

        for j in 1..=i {
            let mut k = j;
            while k <= n {
                let m = k + i;
                if m <= n {
                    let tr = wr * xr[m - 1] - wi * xi[m - 1];
                    let ti = wr * xi[m - 1] + wi * xr[m - 1];
                    xr[m - 1] = xr[k - 1] - tr;
                    xi[m - 1] = xi[k - 1] - ti;
                    xr[k - 1] += tr;
                    xi[k - 1] += ti;
                }
                k += l;
            }
            let help = wr;
            wr = wr + wr * wstpr - wi * wstpi;
            wi = wi + help * wstpi + wi * wstpr;
        }
        i = l;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dg01md_power_of_two() {
        let mut xr = [1.0, 1.0, 0.0, 0.0];
        let mut xi = [0.0; 4];
        assert_eq!(dg01md(Dg01MdIndi::Direct, 4, &mut xr, &mut xi), 0);
        assert_eq!(dg01md(Dg01MdIndi::Inverse, 4, &mut xr, &mut xi), 0);
        // Round-trip scales by N (inverse does not normalize): original [1,1,0,0] -> [4,4,0,0]
        assert!((xr[0] - 4.0).abs() < 1e-6);
        assert!((xr[1] - 4.0).abs() < 1e-6);
    }
}
