//! Example: DE01OD — convolution of two real sequences using FFT.
//!
//! N must be a power of 2 (e.g. 8). Result overwrites the first sequence.

use slicot_rs::de01::de01od::{de01od, De01OdConv};

fn main() {
    let n = 8;
    let mut a = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // impulse
    let mut b = vec![1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // step-like
    let info = de01od(De01OdConv::Convolution, n, &mut a, &mut b);
    assert_eq!(info, 0);
    println!("Convolution (first 8 samples): {:?}", &a[..n]);
}
