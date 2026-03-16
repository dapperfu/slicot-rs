//! Integration test: Rust SLICOT routines must match Fortran output 1:1 when the Fortran
//! reference is available.
//!
//! Set `SLICOT_EXAMPLES_DIR` to the directory containing the built Fortran example executables
//! (e.g. TAB01ND). If unset or the executable is missing, the test is skipped so CI without
//! Fortran still passes.
//!
//! Run with Fortran comparison:
//!   SLICOT_EXAMPLES_DIR=/path/to/examples cargo test --test fortran_1to1_compare
//!
//! Run without (Rust-only, test skips comparison):
//!   cargo test --test fortran_1to1_compare

use std::io::Cursor;
use std::path::Path;
use std::process::Command;

const REL_TOL: f64 = 1e-10;

/// AB01ND: run Rust, optionally run Fortran TAB01ND, assert 1:1 match (INFO, NCONT, INDCON, A_cont, B_cont).
#[test]
fn ab01nd_1to1_fortran() {
    // Same input as fuzz corpus seed (SLICOT-style .dat)
    let dat_content = r#" AB01ND EXAMPLE PROGRAM DATA
   3     2     0.0     I
  -1.0   0.0   0.0
  -2.0  -2.0  -2.0
  -1.0   0.0  -3.0
   1.0   0.0   0.0
   0.0   2.0   1.0
"#;

    let dat = match slicot_rs::slicot_io::parse_ab01nd_dat(Cursor::new(dat_content)) {
        Ok(d) => d,
        Err(e) => panic!("parse_ab01nd_dat failed: {}", e),
    };

    let mut a = dat.a.clone();
    let mut b = dat.b.clone();
    let mut ncont = 0usize;
    let mut indcon = 0usize;
    let mut nblk = vec![0i32; dat.n];
    let jobz = match dat.jobz {
        'N' | 'n' => slicot_rs::ab01::ab01nd::JobZ::No,
        'F' | 'f' => slicot_rs::ab01::ab01nd::JobZ::Factored,
        'I' | 'i' => slicot_rs::ab01::ab01nd::JobZ::Init,
        _ => slicot_rs::ab01::ab01nd::JobZ::No,
    };

    let info = slicot_rs::ab01::ab01nd::ab01nd(
        jobz,
        dat.n,
        dat.m,
        &mut a,
        &mut b,
        &mut ncont,
        &mut indcon,
        &mut nblk,
        None,
        dat.tol,
    );

    let rust_res = slicot_rs::slicot_io::Ab01ndRes {
        info,
        ncont,
        indcon,
        nblk: nblk.clone(),
        a_cont: a.clone(),
        b_cont: b.clone(),
        z: None,
    };

    let fortran_res = match run_fortran_ab01nd(&dat_content) {
        Some(r) => r,
        None => {
            eprintln!("Fortran TAB01ND not available (set SLICOT_EXAMPLES_DIR and build Fortran); skipping 1:1 comparison");
            return;
        }
    };

    assert_eq!(fortran_res.info, rust_res.info, "INFO must match Fortran 1:1");
    assert_eq!(fortran_res.ncont, rust_res.ncont, "NCONT must match Fortran 1:1");
    assert_eq!(fortran_res.indcon, rust_res.indcon, "INDCON must match Fortran 1:1");

    if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fortran_res.a_cont, &rust_res.a_cont, REL_TOL) {
        panic!("A_cont differs from Fortran at ({}, {})", i, j);
    }
    if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fortran_res.b_cont, &rust_res.b_cont, REL_TOL) {
        panic!("B_cont differs from Fortran at ({}, {})", i, j);
    }
}

fn run_fortran_ab01nd(dat_content: &str) -> Option<slicot_rs::slicot_io::Ab01ndRes> {
    let exe_dir = std::env::var("SLICOT_EXAMPLES_DIR").ok()?;
    let tab01nd = Path::new(&exe_dir).join("TAB01ND");
    if !tab01nd.is_file() {
        return None;
    }

    let mut dat_path = std::env::temp_dir();
    dat_path.push("slicot_test_ab01nd.dat");
    let mut res_path = std::env::temp_dir();
    res_path.push("slicot_test_ab01nd.res");

    std::fs::write(&dat_path, dat_content).ok()?;
    let out = Command::new(&tab01nd)
        .arg(&dat_path)
        .arg(&res_path)
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let res_bytes = std::fs::read(&res_path).ok()?;
    slicot_rs::slicot_io::parse_ab01nd_res(Cursor::new(res_bytes)).ok()
}
