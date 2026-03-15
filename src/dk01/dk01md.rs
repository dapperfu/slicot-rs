//! DK01MD — Apply anti-aliasing window to a real signal (SLICOT DK01MD)
//
// TYPE: 'M' Hamming, 'N' Hann, 'Q' Quadratic. A is overwritten with windowed signal.

use std::f64::consts::PI;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dk01MdType {
    /// Hamming: 0.54 + 0.46*cos(pi*(i-1)/(N-1))
    Hamming,
    /// Hann: 0.5*(1 + cos(pi*(i-1)/(N-1)))
    Hann,
    /// Quadratic (piecewise formula)
    Quadratic,
}

/// Applies the selected window to A in place. Returns 0 on success; < 0 invalid argument.
pub fn dk01md(typ: Dk01MdType, n: usize, a: &mut [f64]) -> i32 {
    if n == 0 {
        return -2;
    }
    if a.len() < n {
        return -3;
    }

    let fn_ = (n - 1) as f64;
    if fn_ <= 0.0 {
        return 0;
    }

    match typ {
        Dk01MdType::Hamming => {
            let temp = PI / fn_;
            for i in 0..n {
                a[i] *= 0.54 + 0.46 * (temp * (i as f64)).cos();
            }
        }
        Dk01MdType::Hann => {
            let temp = PI / fn_;
            for i in 0..n {
                a[i] *= 0.5 * (1.0 + (temp * (i as f64)).cos());
            }
        }
        Dk01MdType::Quadratic => {
            let n1 = (n - 1) / 2 + 1;
            for i in 0..n {
                let buf = (i as f64) / fn_;
                let temp = buf * buf;
                if i < n1 {
                    a[i] *= (1.0 - 2.0 * temp) * (1.0 - buf);
                } else {
                    a[i] *= 2.0 * (1.0 - buf * temp);
                }
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dk01md_hamming() {
        let mut a = [1.0, 1.0, 1.0];
        assert_eq!(dk01md(Dk01MdType::Hamming, 3, &mut a), 0);
        assert!((a[0] - 1.0).abs() < 1e-14);
        assert!((a[1] - (0.54 + 0.46 * (PI / 2.0_f64).cos())).abs() < 1e-14);
        assert!((a[2] - (0.54 + 0.46 * PI.cos())).abs() < 1e-14);
    }

    #[test]
    fn test_dk01md_hann() {
        let mut a = [1.0, 1.0];
        assert_eq!(dk01md(Dk01MdType::Hann, 2, &mut a), 0);
        assert!((a[0] - 1.0).abs() < 1e-14);
        assert!((a[1] - 0.0).abs() < 1e-14);
    }

    #[test]
    fn test_dk01md_quadratic() {
        let mut a = [1.0, 1.0, 1.0];
        assert_eq!(dk01md(Dk01MdType::Quadratic, 3, &mut a), 0);
        assert!(a[0] >= 0.0 && a[2] >= 0.0);
    }
}
