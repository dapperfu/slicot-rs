//! Fuzz target: parse bytes as AB01ND .dat, run Rust ab01nd, optionally run Fortran and compare.
//! Corpus: seed with SLICOT-style .dat content in fuzz/corpus/ab01nd/.
//! Set SLICOT_EXAMPLES_DIR to the directory containing TAB01ND for Fortran comparison.

#![no_main]

use libfuzzer_sys::fuzz_target;
use slicot_rs::ab01::ab01nd::ab01nd;
use slicot_rs::ab01::ab01nd::JobZ;
use slicot_rs::slicot_io::{parse_ab01nd_dat, rel_tol_eq_matrix, Ab01ndRes};
use std::io::Cursor;
use std::process::Command;

const MAX_N: usize = 20;
const MAX_M: usize = 20;
const REL_TOL: f64 = 1e-10;

fuzz_target!(|data: &[u8]| {
    fuzz_ab01nd(data);
});

fn fuzz_ab01nd(data: &[u8]) {
    let Ok(dat) = parse_ab01nd_dat(Cursor::new(data)) else {
        return;
    };
    if dat.n > MAX_N || dat.m > MAX_M {
        return;
    }
    let mut a = dat.a.clone();
    let mut b = dat.b.clone();
    let mut ncont = 0usize;
    let mut indcon = 0usize;
    let mut nblk = vec![0i32; dat.n];
    let jobz = match dat.jobz {
        'N' | 'n' => JobZ::No,
        'F' | 'f' => JobZ::Factored,
        'I' | 'i' => JobZ::Init,
        _ => JobZ::No,
    };
    let info = ab01nd(
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
    if info != 0 {
        return;
    }
    if let Some((fortran_res, rust_res)) = run_fortran_and_build_rust_res(&dat, info, ncont, indcon, &nblk, &a, &b) {
        if fortran_res.info != rust_res.info {
            panic!("INFO mismatch: fortran {} rust {}", fortran_res.info, rust_res.info);
        }
        if fortran_res.ncont != rust_res.ncont {
            panic!("NCONT mismatch: fortran {} rust {}", fortran_res.ncont, rust_res.ncont);
        }
        if let Some((i, j)) = rel_tol_eq_matrix(&fortran_res.a_cont, &rust_res.a_cont, REL_TOL) {
            panic!("A_cont diff at ({}, {})", i, j);
        }
        if let Some((i, j)) = rel_tol_eq_matrix(&fortran_res.b_cont, &rust_res.b_cont, REL_TOL) {
            panic!("B_cont diff at ({}, {})", i, j);
        }
    }
}

fn run_fortran_and_build_rust_res(
    dat: &slicot_rs::slicot_io::Ab01ndDat,
    info: i32,
    ncont: usize,
    indcon: usize,
    nblk: &[i32],
    a_cont: &nalgebra::DMatrix<f64>,
    b_cont: &nalgebra::DMatrix<f64>,
) -> Option<(Ab01ndRes, Ab01ndRes)> {
    let exe_dir = std::env::var("SLICOT_EXAMPLES_DIR").ok()?;
    let tab01nd = std::path::Path::new(&exe_dir).join("TAB01ND");
    if !tab01nd.is_file() {
        return None;
    }
    let mut dat_file = std::env::temp_dir();
    dat_file.push("slicot_fuzz_ab01nd.dat");
    let mut out_file = std::env::temp_dir();
    out_file.push("slicot_fuzz_ab01nd.res");
    let mut f = std::fs::File::create(&dat_file).ok()?;
    slicot_rs::slicot_io::write_ab01nd_dat(&mut f, dat).ok()?;
    drop(f);
    let out = Command::new(&tab01nd)
        .arg(&dat_file)
        .arg(&out_file)
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let fortran_res = slicot_rs::slicot_io::parse_ab01nd_res(Cursor::new(std::fs::read(&out_file).ok()?)).ok()?;
    let rust_res = Ab01ndRes {
        info,
        ncont,
        indcon,
        nblk: nblk.to_vec(),
        a_cont: a_cont.clone(),
        b_cont: b_cont.clone(),
        z: None,
    };
    Some((fortran_res, rust_res))
}
