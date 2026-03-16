//! # DE01OD — Convolution of two real sequences via FFT
//!
//! ## What the SLICOT routine does
//!
//! **DE01OD** computes the **convolution** (or correlation) of two real sequences of length \(N\)
//! using the FFT. The result is returned in the first sequence; the second is used as workspace.
//! \(N\) must be a power of two. This is the standard FFT-based convolution: linear convolution
//! of two length-\(N\) sequences has length \(2N-1\); the routine uses a length-\(2N\) FFT and
//! returns the first \(N\) samples (or the full linear convolution depending on the routine’s
//! convention).
//!
//! ## Why it exists
//!
//! In signal processing and control, convolution appears in filtering, system response
//! computation, and correlation-based identification. Doing convolution in the frequency domain
//! via FFT is \(O(N \log N)\) instead of \(O(N^2)\) for the direct sum, so DE01OD is the
//! preferred building block in SLICOT for signal-processing and identification routines.
//!
//! ## When to use it
//!
//! - You need the convolution of two real sequences and have power-of-two length.
//! - You are building filters or computing impulse/step responses in the discrete-time domain.
//!
//! This example convolves a simple impulse with a step-like sequence and prints the first
//! eight samples.

use slicot_rs::de01::de01od::{de01od, De01OdConv};

fn main() {
    let n = 8;
    let mut a = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // impulse
    let mut b = vec![1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // step-like
    let info = de01od(De01OdConv::Convolution, n, &mut a, &mut b);
    assert_eq!(info, 0);
    println!("Convolution (first {} samples): {:?}", n, &a[..n]);
}
