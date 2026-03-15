//! Integration tests: validate Rust implementations against Fortran routine specifications.
//!
//! Each test checks mathematical invariants or 1:1 behavior documented in the SLICOT Fortran
//! source (SLICOT-Reference/src/*.f). Running `cargo test` runs both unit tests and these
//! spec-consistency tests.

use nalgebra::DMatrix;
use slicot_rs::ma01::ma01ad::ma01ad;
use slicot_rs::ma01::ma01dd::ma01dd;
use slicot_rs::ma02::ma02ad::{self, Ma02AdJob};
use slicot_rs::mb01::mb01ru::{self, Mb01RuTrans, Mb01RuUplo};

/// MA01AD: (YR + i*YI)^2 must equal XR + i*XI (Fortran PURPOSE).
#[test]
fn ma01ad_spec_complex_square_root() {
    let cases = [(4.0, 0.0), (0.0, 1.0), (1.0, 1.0), (-1.0, 0.0), (2.0, -1.0)];
    for (xr, xi) in cases {
        let mut yr = 0.0;
        let mut yi = 0.0;
        assert_eq!(ma01ad(xr, xi, &mut yr, &mut yi), 0);
        // (yr + i*yi)^2 = (yr^2 - yi^2) + i*2*yr*yi
        let re = yr * yr - yi * yi;
        let im = 2.0 * yr * yi;
        assert!(
            (re - xr).abs() < 1e-13 && (im - xi).abs() < 1e-13,
            "MA01AD (xr={}, xi={}) -> (yr={}, yi={}), (yr+i*yi)^2 = ({}, {})",
            xr, xi, yr, yi, re, im
        );
        assert!(yr >= -1e-15, "MA01AD: YR >= 0");
        if xi != 0.0 {
            assert_eq!((yi > 0.0), (xi > 0.0), "MA01AD: SIGN(YI) = SIGN(XI)");
        }
    }
}

/// MA01DD: D = min(|A1-A2|, |1/A1-1/A2|) (Fortran METHOD).
#[test]
fn ma01dd_spec_chordal_metric() {
    let eps = 1e-15_f64;
    let safemn = 1e-308_f64;
    // Same point -> D = 0
    let mut d = -1.0;
    ma01dd(1.0, 0.0, 1.0, 0.0, eps, safemn, &mut d);
    assert!(d >= 0.0 && d < 1e-13);
    // Different: D should be min of direct distance and chordal
    ma01dd(1.0, 0.0, 2.0, 0.0, eps, safemn, &mut d);
    let direct = 1.0_f64; // |1-2|
    assert!(d > 0.0 && d <= direct + 1e-10);
}

/// MA02AD: transpose then transpose again recovers original (Fortran: B = A or A').
#[test]
fn ma02ad_spec_transpose_twice() {
    let a = DMatrix::from_row_slice(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let mut b = DMatrix::zeros(3, 2);
    let mut c = DMatrix::zeros(2, 3);
    ma02ad::ma02ad(Ma02AdJob::All, &a, &mut b);
    ma02ad::ma02ad(Ma02AdJob::All, &b, &mut c);
    for i in 0..2 {
        for j in 0..3 {
            assert!((c[(i, j)] - a[(i, j)]).abs() < 1e-15, "MA02AD transpose twice");
        }
    }
}

/// MB01RU: R = alpha*R + beta*A*X*A' (Fortran formula, TRANS='N').
#[test]
fn mb01ru_spec_symmetric_update() {
    let m = 2_usize;
    let n = 2_usize;
    let mut r = [1.0, 0.0, 0.0, 1.0];
    let a = [1.0, 0.0, 0.0, 1.0];
    let x = [1.0, 0.0, 0.0, 1.0];
    let mut dwork = vec![0.0; m * n];
    let info = mb01ru::mb01ru(
        Mb01RuUplo::Upper,
        Mb01RuTrans::NoTrans,
        m,
        n,
        1.0,
        1.0,
        &mut r,
        2,
        &a,
        2,
        &x,
        2,
        &mut dwork,
    );
    assert_eq!(info, 0);
    // R := R + A*X*A' = I + I*I*I = 2*I (upper stored: r00=2, r01=0, r11=2)
    assert!((r[0] - 2.0).abs() < 1e-14);
    assert!((r[2] - 0.0).abs() < 1e-14);
    assert!((r[3] - 2.0).abs() < 1e-14);
}
