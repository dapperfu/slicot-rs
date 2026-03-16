//! DE01OD — Convolution or deconvolution of two real signals using FFT (SLICOT DE01OD).
//!
//! Uses DG01MD (FFT). N must be a power of 2, N >= 2.

use crate::dg01::dg01md::dg01md;
use crate::dg01::dg01md::Dg01MdIndi;

fn is_power_of_two(n: usize) -> bool {
    n >= 2 && (n & (n - 1)) == 0
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum De01OdConv {
    Convolution,
    Deconvolution,
}

/// Convolution or deconvolution of real signals A and B. Result overwrites A; B is overwritten.
/// Returns 0 on success; <0 invalid argument.
pub fn de01od(conv: De01OdConv, n: usize, a: &mut [f64], b: &mut [f64]) -> i32 {
    if !matches!(conv, De01OdConv::Convolution | De01OdConv::Deconvolution) {
        return -1;
    }
    if !is_power_of_two(n) {
        return -2;
    }
    if a.len() < n || b.len() < n {
        return -3;
    }
    let lconv = conv == De01OdConv::Convolution;
    let mut info = dg01md(Dg01MdIndi::Direct, n, a, b);
    if info != 0 {
        return info;
    }
    let nd2p1 = n / 2 + 1;
    if lconv {
        a[0] *= b[0];
    } else {
        a[0] = if b[0] == 0.0 { 0.0 } else { a[0] / b[0] };
    }
    for i in 0..=(n - nd2p1) {
        let j = nd2p1 - i;
        let kj = nd2p1 + i;
        if j == 0 || kj > n {
            break;
        }
        let j0 = j - 1;
        let kj0 = kj - 1;
        let ac = 0.5 * (a[j0] + a[kj0]);
        let as_ = 0.5 * (b[j0] - b[kj0]);
        let bc = 0.5 * (b[kj0] + b[j0]);
        let bs = 0.5 * (a[kj0] - a[j0]);
        let (cr, ci) = if lconv {
            (ac * bc - as_ * bs, as_ * bc + ac * bs)
        } else {
            if bc.abs().max(bs.abs()) == 0.0 {
                (0.0, 0.0)
            } else {
                let denom = bc * bc + bs * bs;
                ((ac * bc + as_ * bs) / denom, (as_ * bc - ac * bs) / denom)
            }
        };
        a[j0] = cr;
        b[j0] = ci;
        a[kj0] = cr;
        b[kj0] = -ci;
    }
    b[0] = 0.0;
    info = dg01md(Dg01MdIndi::Inverse, n, a, b);
    if info != 0 {
        return info;
    }
    let scale = 1.0 / (n as f64);
    for i in 0..n {
        a[i] *= scale;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_de01od_conv_n4() {
        let mut a = [1.0, 0.0, 0.0, 0.0];
        let mut b = [1.0, 0.0, 0.0, 0.0];
        let info = de01od(De01OdConv::Convolution, 4, &mut a, &mut b);
        assert_eq!(info, 0);
        assert!((a[0] - 1.0).abs() < 1e-10);
    }
}
