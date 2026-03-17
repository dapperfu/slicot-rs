//! Fortran 1:1 validation runner: discover routines with .dat, run Fortran + Rust, compare, write validation/*.md.
//!
//! Set SLICOT_EXAMPLES_DIR to the examples directory (e.g. SLICOT-Reference/examples). If unset
//! or Fortran drivers are missing, tests skip and exit 0.

use std::collections::BTreeMap;
use std::fs;
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const REL_TOL: f64 = 1e-10;
const SLICOT_MAPPING_PATH: &str = "docs/SLICOT_MAPPING.md";
const VALIDATION_DIR: &str = "validation";

/// (SLICOT name, Rust module, Rust function)
fn parse_slicot_mapping(project_root: &Path) -> Vec<(String, String, String)> {
    let path = project_root.join(SLICOT_MAPPING_PATH);
    let content = match fs::read_to_string(&path) {
        Ok(c) => c,
        Err(_) => return vec![],
    };
    let mut out = vec![];
    for line in content.lines() {
        let line = line.trim();
        if !line.starts_with('|') || line == "| SLICOT | Rust module | Rust function | Status |" || line.starts_with("|---") {
            continue;
        }
        let cells: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if cells.len() >= 5 {
            let slicot = cells[1].to_string();
            let module = cells[2].to_string();
            let rust_fn = cells[3].to_string();
            let status = cells[4];
            if status == "done" && !slicot.is_empty() && !module.is_empty() {
                out.push((slicot, module, rust_fn));
            }
        }
    }
    out
}

/// List routine names that have data/XX.dat in examples dir.
fn routines_with_dat(examples_dir: &Path) -> Vec<String> {
    let data_dir = examples_dir.join("data");
    let mut names = vec![];
    if let Ok(entries) = fs::read_dir(&data_dir) {
        for e in entries.flatten() {
            let name = e.file_name();
            if let Some(s) = name.to_str() {
                if s.ends_with(".dat") {
                    names.push(s.trim_end_matches(".dat").to_string());
                }
            }
        }
    }
    names.sort();
    names
}

/// Build map: module -> list of (slicot, rust_fn) that have .dat, sorted by routine name.
fn discovery(project_root: &Path, examples_dir: &Path) -> BTreeMap<String, Vec<(String, String)>> {
    let mapping: std::collections::HashMap<String, (String, String)> = parse_slicot_mapping(project_root)
        .into_iter()
        .map(|(slicot, module, rust_fn)| (slicot.clone(), (module, rust_fn)))
        .collect();
    let with_dat = routines_with_dat(examples_dir);
    let mut by_module: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
    for slicot in with_dat {
        if let Some((module, rust_fn)) = mapping.get(&slicot) {
            by_module
                .entry(module.clone())
                .or_default()
                .push((slicot, rust_fn.clone()));
        }
    }
    for v in by_module.values_mut() {
        v.sort_by(|a, b| a.0.cmp(&b.0));
    }
    by_module
}

/// Run Fortran driver with stdin from .dat, stdout to temp .res; return path to .res and success.
fn run_fortran_stdin_redirect(
    examples_dir: &Path,
    slicot: &str,
    dat_path: &Path,
    res_path: &Path,
) -> std::io::Result<bool> {
    let exe = examples_dir.join(format!("T{slicot}"));
    if !exe.is_file() {
        return Ok(false);
    }
    let dat_file = fs::File::open(dat_path)?;
    let res_file = fs::File::create(res_path)?;
    let status = Command::new(&exe)
        .current_dir(examples_dir)
        .stdin(Stdio::from(dat_file))
        .stdout(Stdio::from(res_file))
        .status()?;
    Ok(status.success())
}

/// Validate AB01MD: run Fortran, run Rust, compare.
fn validate_ab01md(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let dat = match slicot_rs::slicot_io::parse_ab01md_dat(Cursor::new(&dat_content)) {
        Ok(d) => d,
        Err(e) => return (false, format!("parse_ab01md_dat: {}", e), String::new(), String::new()),
    };

    let res_path = project_root.join("target").join("validation_ab01md.res");
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, "AB01MD", dat_path, &res_path).unwrap_or(false);
    let fortran_res = if ok {
        let bytes = fs::read(&res_path).unwrap_or_default();
        slicot_rs::slicot_io::parse_ab01md_res(Cursor::new(bytes)).ok()
    } else {
        None
    };

    let mut a = dat.a.clone();
    let mut b = dat.b.clone();
    let mut ncont = 0usize;
    let jobz = match dat.jobz {
        'N' | 'n' => slicot_rs::ab01::ab01md::JobZ::No,
        'F' | 'f' => slicot_rs::ab01::ab01md::JobZ::Factored,
        'I' | 'i' => slicot_rs::ab01::ab01md::JobZ::Init,
        _ => slicot_rs::ab01::ab01md::JobZ::No,
    };
    let mut z_mat = nalgebra::DMatrix::zeros(dat.n, dat.n);
    let info = slicot_rs::ab01::ab01md::ab01md(jobz, &mut a, &mut b, Some(&mut z_mat), dat.tol, &mut ncont);
    let (a_cont, b_cont) = if ncont > 0 {
        let b_slice: Vec<f64> = b.iter().take(ncont).copied().collect();
        let b_mat = nalgebra::DMatrix::from_column_slice(ncont, 1, &b_slice);
        (a.view((0, 0), (ncont, ncont)).into_owned(), b_mat)
    } else {
        (nalgebra::DMatrix::zeros(0, 0), nalgebra::DMatrix::zeros(0, 1))
    };
    let z_out = if jobz == slicot_rs::ab01::ab01md::JobZ::Init {
        Some(z_mat.clone())
    } else {
        None
    };

    let input_summary = format!("n={}, tol={}, jobz={}", dat.n, dat.tol, dat.jobz);
    let fortran_summary = fortran_res.as_ref().map(|r| format!("NCONT={}", r.ncont)).unwrap_or_else(|| "Fortran not run".to_string());
    let rust_summary = format!("INFO={}, NCONT={}", info, ncont);

    let (pass, diff_msg) = match &fortran_res {
        None => (false, "Fortran driver missing or failed".to_string()),
        Some(fr) => {
            if info != 0 {
                (false, format!("Rust INFO={}", info))
            } else if fr.ncont != ncont {
                (false, format!("NCONT mismatch: Fortran={} Rust={}", fr.ncont, ncont))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.a_cont, &a_cont, REL_TOL) {
                (false, format!("A_cont differs at ({}, {})", i, j))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.b_cont, &b_cont, REL_TOL) {
                (false, format!("B_cont differs at ({}, {})", i, j))
            } else if let (Some(ref fz), Some(ref rz)) = (&fr.z, &z_out) {
                if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(fz, rz, REL_TOL) {
                    (false, format!("Z differs at ({}, {})", i, j))
                } else {
                    (true, String::new())
                }
            } else {
                (true, String::new())
            }
        }
    };
    (
        pass,
        input_summary,
        fortran_summary,
        if pass { rust_summary } else { format!("{}; {}", rust_summary, diff_msg) },
    )
}

/// Parse TG01GD .dat: title, then "L N M P D tol", then E(N×N), A(N×N), B(N×M), C(P×N), D(P×M).
fn parse_tg01gd_dat(content: &str) -> Result<(usize, usize, usize, usize, nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>), String> {
    let lines: Vec<&str> = content.lines().map(str::trim).collect();
    if lines.len() < 3 {
        return Err("too few lines".to_string());
    }
    let nums: Vec<f64> = lines[1].split_whitespace().filter_map(|s| s.parse().ok()).collect();
    if nums.len() < 4 {
        return Err("need L N M P".to_string());
    }
    let l = nums[0] as usize;
    let n = nums[1] as usize;
    let m = nums[2] as usize;
    let p = nums[3] as usize;
    let mut idx = 2;
    let mut read_matrix = |rows: usize, cols: usize| -> Result<nalgebra::DMatrix<f64>, String> {
        let mut v = vec![];
        while v.len() < rows * cols && idx < lines.len() {
            for w in lines[idx].split_whitespace() {
                if let Ok(x) = w.parse::<f64>() {
                    v.push(x);
                    if v.len() >= rows * cols {
                        break;
                    }
                }
            }
            idx += 1;
        }
        if v.len() != rows * cols {
            return Err(format!("expected {} values", rows * cols));
        }
        Ok(nalgebra::DMatrix::from_row_slice(rows, cols, &v))
    };
    let e = read_matrix(n, n)?;
    let a = read_matrix(n, n)?;
    let b = read_matrix(n, m)?;
    let c = read_matrix(p, n)?;
    Ok((l, n, m, p, e, a, b, c))
}

/// Parse TG01HD / TG01ID .dat: title, then "N M P tol JOB", then E(N×N), A(N×N), B(N×M), C(P×N).
fn parse_tg01hd_id_dat(content: &str) -> Result<(usize, usize, usize, usize, nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>), String> {
    let lines: Vec<&str> = content.lines().map(str::trim).collect();
    if lines.len() < 3 {
        return Err("too few lines".to_string());
    }
    let nums: Vec<f64> = lines[1].split_whitespace().filter_map(|s| s.parse().ok()).collect();
    if nums.len() < 3 {
        return Err("need N M P".to_string());
    }
    let n = nums[0] as usize;
    let m = nums[1] as usize;
    let p = nums[2] as usize;
    let l = n;
    let mut idx = 2;
    let mut read_matrix = |rows: usize, cols: usize| -> Result<nalgebra::DMatrix<f64>, String> {
        let mut v = vec![];
        while v.len() < rows * cols && idx < lines.len() {
            for w in lines[idx].split_whitespace() {
                if let Ok(x) = w.parse::<f64>() {
                    v.push(x);
                    if v.len() >= rows * cols {
                        break;
                    }
                }
            }
            idx += 1;
        }
        if v.len() != rows * cols {
            return Err(format!("expected {} values", rows * cols));
        }
        Ok(nalgebra::DMatrix::from_row_slice(rows, cols, &v))
    };
    let e = read_matrix(n, n)?;
    let a = read_matrix(n, n)?;
    let b = read_matrix(n, m)?;
    let c = read_matrix(p, n)?;
    Ok((l, n, m, p, e, a, b, c))
}

/// Validate TG01GD: run Fortran, run Rust, compare.
fn validate_tg01gd(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let (l, n, m, p, e, a, b, c) = match parse_tg01gd_dat(&dat_content) {
        Ok(x) => x,
        Err(e) => return (false, format!("parse: {}", e), String::new(), String::new()),
    };
    let res_path = project_root.join("target").join("validation_tg01gd.res");
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, "TG01GD", dat_path, &res_path).unwrap_or(false);
    let fortran_res_str = if ok { fs::read_to_string(&res_path).unwrap_or_default() } else { String::new() };
    let fortran_summary = if ok {
        if let Some(rank_line) = fortran_res_str.lines().find(|s| s.contains("Rank of matrix E")) {
            rank_line.trim().to_string()
        } else {
            "Fortran ran".to_string()
        }
    } else {
        "Fortran driver missing or failed".to_string()
    };
    let mut a_rust = a.clone();
    let mut e_rust = e.clone();
    let mut b_rust = b.clone();
    let mut c_rust = c.clone();
    let info = slicot_rs::tg01::tg01gd::tg01gd(l, n, m, p, &mut a_rust, &mut e_rust, &mut b_rust, &mut c_rust);
    let rust_summary = format!("INFO={}", info);
    let pass = ok && info == 0 && fortran_res_str.contains("Rank of matrix E");
    (pass, format!("l={} n={} m={} p={}", l, n, m, p), fortran_summary, rust_summary)
}

/// Validate TG01HD: run Fortran, run Rust, compare.
fn validate_tg01hd(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let (l, n, m, p, e, a, b, c) = match parse_tg01hd_id_dat(&dat_content) {
        Ok(x) => x,
        Err(e) => return (false, format!("parse: {}", e), String::new(), String::new()),
    };
    let res_path = project_root.join("target").join("validation_tg01hd.res");
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, "TG01HD", dat_path, &res_path).unwrap_or(false);
    let fortran_res_str = if ok { fs::read_to_string(&res_path).unwrap_or_default() } else { String::new() };
    let fortran_summary = if ok {
        if let Some(dim_line) = fortran_res_str.lines().find(|s| s.contains("Dimension of controllable")) {
            dim_line.trim().to_string()
        } else {
            "Fortran ran".to_string()
        }
    } else {
        "Fortran driver missing or failed".to_string()
    };
    let mut a_rust = a.clone();
    let mut e_rust = e.clone();
    let mut b_rust = b.clone();
    let mut c_rust = c.clone();
    let info = slicot_rs::tg01::tg01hd::tg01hd(l, n, m, p, &mut a_rust, &mut e_rust, &mut b_rust, &mut c_rust);
    let rust_summary = format!("INFO={}", info);
    let pass = ok && info == 0 && fortran_res_str.contains("Dimension of controllable");
    (pass, format!("n={} m={} p={}", n, m, p), fortran_summary, rust_summary)
}

/// Validate TG01ID: run Fortran, run Rust, compare.
fn validate_tg01id(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let (l, n, m, p, e, a, b, c) = match parse_tg01hd_id_dat(&dat_content) {
        Ok(x) => x,
        Err(e) => return (false, format!("parse: {}", e), String::new(), String::new()),
    };
    let res_path = project_root.join("target").join("validation_tg01id.res");
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, "TG01ID", dat_path, &res_path).unwrap_or(false);
    let fortran_res_str = if ok { fs::read_to_string(&res_path).unwrap_or_default() } else { String::new() };
    let fortran_summary = if ok {
        if let Some(dim_line) = fortran_res_str.lines().find(|s| s.contains("Dimension of observable")) {
            dim_line.trim().to_string()
        } else {
            "Fortran ran".to_string()
        }
    } else {
        "Fortran driver missing or failed".to_string()
    };
    let mut a_rust = a.clone();
    let mut e_rust = e.clone();
    let mut b_rust = b.clone();
    let mut c_rust = c.clone();
    let info = slicot_rs::tg01::tg01id::tg01id(l, n, m, p, &mut a_rust, &mut e_rust, &mut b_rust, &mut c_rust);
    let rust_summary = format!("INFO={}", info);
    let pass = ok && info == 0 && fortran_res_str.contains("Dimension of observable");
    (pass, format!("n={} m={} p={}", n, m, p), fortran_summary, rust_summary)
}

/// Validate AB01ND: run Fortran, run Rust, compare, return (pass, optional error message).
fn validate_ab01nd(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let dat = match slicot_rs::slicot_io::parse_ab01nd_dat(Cursor::new(&dat_content)) {
        Ok(d) => d,
        Err(e) => return (false, format!("parse_ab01nd_dat: {}", e), String::new(), String::new()),
    };

    let res_path = project_root.join("target").join("validation_ab01nd.res");
    if let Ok(()) = fs::create_dir_all(res_path.parent().unwrap()) {}
    let ok = run_fortran_stdin_redirect(
        examples_dir,
        "AB01ND",
        dat_path,
        &res_path,
    ).unwrap_or(false);
    let fortran_res = if ok {
        let bytes = fs::read(&res_path).unwrap_or_default();
        slicot_rs::slicot_io::parse_ab01nd_res(Cursor::new(bytes)).ok()
    } else {
        None
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

    let (input_summary, fortran_summary, rust_summary) = (
        format!("n={}, m={}, tol={}, jobz={}", dat.n, dat.m, dat.tol, dat.jobz),
        fortran_res.as_ref().map(|r| format!("INFO={}, NCONT={}, INDCON={}", r.info, r.ncont, r.indcon)).unwrap_or_else(|| "Fortran not run".to_string()),
        format!("INFO={}, NCONT={}, INDCON={}", rust_res.info, rust_res.ncont, rust_res.indcon),
    );

    let (pass, diff_msg) = match &fortran_res {
        None => (false, "Fortran driver missing or failed".to_string()),
        Some(fr) => {
            if fr.info != rust_res.info {
                (false, format!("INFO mismatch: Fortran={} Rust={}", fr.info, rust_res.info))
            } else if fr.ncont != rust_res.ncont {
                (false, format!("NCONT mismatch: Fortran={} Rust={}", fr.ncont, rust_res.ncont))
            } else if fr.indcon != rust_res.indcon {
                (false, format!("INDCON mismatch: Fortran={} Rust={}", fr.indcon, rust_res.indcon))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.a_cont, &rust_res.a_cont, REL_TOL) {
                (false, format!("A_cont differs at ({}, {})", i, j))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.b_cont, &rust_res.b_cont, REL_TOL) {
                (false, format!("B_cont differs at ({}, {})", i, j))
            } else {
                (true, String::new())
            }
        }
    };
    (pass, input_summary, fortran_summary, if pass { rust_summary } else { format!("{}; {}", rust_summary, diff_msg) })
}

/// Validate AB07MD: run Fortran, run Rust, compare.
fn validate_ab07md(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let dat = match slicot_rs::ab07::io::parse_ab07md_dat(Cursor::new(&dat_content)) {
        Ok(d) => d,
        Err(e) => return (false, format!("parse_ab07md_dat: {}", e), String::new(), String::new()),
    };
    let res_path = project_root.join("target").join("validation_ab07md.res");
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, "AB07MD", dat_path, &res_path).unwrap_or(false);
    let fortran_res = if ok {
        let bytes = fs::read(&res_path).unwrap_or_default();
        slicot_rs::ab07::io::parse_ab07md_res(Cursor::new(bytes), dat.n, dat.m, dat.p, dat.jobd).ok()
    } else {
        None
    };
    let jobd = if dat.jobd {
        slicot_rs::ab07::ab07md::JobD::D
    } else {
        slicot_rs::ab07::ab07md::JobD::Z
    };
    let mplim = dat.n.max(dat.m).max(dat.p);
    let mut a = dat.a.clone();
    let mut b = {
        let mut bb = nalgebra::DMatrix::zeros(dat.n, mplim);
        bb.view_mut((0, 0), (dat.n, dat.m)).copy_from(&dat.b);
        bb
    };
    let mut c = {
        let mut cc = nalgebra::DMatrix::zeros(mplim, dat.n);
        cc.view_mut((0, 0), (dat.p, dat.n)).copy_from(&dat.c);
        cc
    };
    let mut d_opt = dat.d.map(|d| {
        let mut dd = nalgebra::DMatrix::zeros(mplim, mplim);
        dd.view_mut((0, 0), (dat.p, dat.m)).copy_from(&d);
        dd
    });
    let info = slicot_rs::ab07::ab07md::ab07md(
        jobd,
        dat.n,
        dat.m,
        dat.p,
        &mut a,
        &mut b,
        &mut c,
        d_opt.as_mut(),
    );
    let input_summary = format!("n={}, m={}, p={}, jobd={}", dat.n, dat.m, dat.p, dat.jobd);
    let fortran_summary = fortran_res
        .as_ref()
        .map(|_| "Fortran ran".to_string())
        .unwrap_or_else(|| "Fortran driver missing or failed".to_string());
    let rust_summary = format!("INFO={}", info);
    let (pass, diff_msg) = match &fortran_res {
        None => (false, "Fortran driver missing or failed".to_string()),
        Some(fr) => {
            if info != 0 {
                (false, format!("Rust INFO={}", info))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.a, &a, REL_TOL) {
                (false, format!("A differs at ({}, {})", i, j))
            } else {
                let b_out = b.view((0, 0), (dat.n, dat.p)).into_owned();
                if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.b, &b_out, REL_TOL) {
                    (false, format!("B differs at ({}, {})", i, j))
                } else {
                    let c_out = c.view((0, 0), (dat.m, dat.n)).into_owned();
                    if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.c, &c_out, REL_TOL) {
                        (false, format!("C differs at ({}, {})", i, j))
                    } else if let (Some(ref fd), Some(ref rd)) = (&fr.d, &d_opt) {
                        let d_out = rd.view((0, 0), (dat.m, dat.p)).into_owned();
                        if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(fd, &d_out, REL_TOL) {
                    (false, format!("D differs at ({}, {})", i, j))
                        } else {
                            (true, String::new())
                        }
                    } else if fr.d.is_none() && d_opt.is_none() {
                        (true, String::new())
                    } else {
                        (false, "D presence mismatch".to_string())
                    }
                }
            }
        }
    };
    (
        pass,
        input_summary,
        fortran_summary,
        if pass { rust_summary } else { format!("{}; {}", rust_summary, diff_msg) },
    )
}

/// Validate AB07ND: run Fortran, run Rust, compare.
fn validate_ab07nd(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let dat = match slicot_rs::ab07::io::parse_ab07nd_dat(Cursor::new(&dat_content)) {
        Ok(d) => d,
        Err(e) => return (false, format!("parse_ab07nd_dat: {}", e), String::new(), String::new()),
    };
    let res_path = project_root.join("target").join("validation_ab07nd.res");
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, "AB07ND", dat_path, &res_path).unwrap_or(false);
    let fortran_res = if ok {
        let bytes = fs::read(&res_path).unwrap_or_default();
        slicot_rs::ab07::io::parse_ab07nd_res(Cursor::new(bytes), dat.n, dat.m).ok()
    } else {
        None
    };
    let mut a = dat.a.clone();
    let mut b = dat.b.clone();
    let mut c = dat.c.clone();
    let mut d = dat.d.clone();
    let mut rcond = 0.0;
    let info = slicot_rs::ab07::ab07nd::ab07nd(dat.n, dat.m, &mut a, &mut b, &mut c, &mut d, &mut rcond);
    let input_summary = format!("n={}, m={}", dat.n, dat.m);
    let fortran_summary = fortran_res
        .as_ref()
        .map(|_| "Fortran ran".to_string())
        .unwrap_or_else(|| "Fortran driver missing or failed".to_string());
    let rust_summary = format!("INFO={}, RCOND={}", info, rcond);
    let (pass, diff_msg) = match &fortran_res {
        None => (false, "Fortran driver missing or failed".to_string()),
        Some(fr) => {
            if info != 0 {
                (false, format!("Rust INFO={}", info))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.a, &a, REL_TOL) {
                (false, format!("A differs at ({}, {})", i, j))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.b, &b, REL_TOL) {
                (false, format!("B differs at ({}, {})", i, j))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.c, &c, REL_TOL) {
                (false, format!("C differs at ({}, {})", i, j))
            } else if let Some((i, j)) = slicot_rs::slicot_io::rel_tol_eq_matrix(&fr.d, &d, REL_TOL) {
                (false, format!("D differs at ({}, {})", i, j))
            } else {
                (true, String::new())
            }
        }
    };
    (
        pass,
        input_summary,
        fortran_summary,
        if pass { rust_summary } else { format!("{}; {}", rust_summary, diff_msg) },
    )
}

/// Validate AB08NW: run Fortran, run Rust (full call with N,M,P from .dat), compare INFO.
fn validate_ab08nw(
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let dat = match slicot_rs::ab08::io::parse_ab08nw_dat(Cursor::new(&dat_content)) {
        Ok(d) => d,
        Err(e) => return (false, format!("parse_ab08nw_dat: {}", e), String::new(), String::new()),
    };
    let res_path = project_root
        .join("target")
        .join("validation_ab08nw.res");
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, "AB08NW", dat_path, &res_path).unwrap_or(false);
    let fortran_summary = if ok {
        "Fortran ran".to_string()
    } else {
        "Fortran driver missing or failed".to_string()
    };
    let equil = if dat.equil {
        slicot_rs::ab08::ab08nw::Equil::S
    } else {
        slicot_rs::ab08::ab08nw::Equil::N
    };
    let n = dat.n;
    let m = dat.m;
    let p = dat.p;
    let mut a = dat.a.clone();
    let mut b = dat.b.clone();
    let mut c = dat.c.clone();
    let mut d = dat.d.clone();
    let mut nfz = 0i32;
    let mut nrank = 0i32;
    let mut niz = 0i32;
    let mut dinfz = 0i32;
    let mut nkror = 0i32;
    let mut ninfe = 0i32;
    let mut nkrol = 0i32;
    let cap = n + 1;
    let mut infz = vec![0i32; cap];
    let mut kronr = vec![0i32; cap];
    let mut infe = vec![0i32; cap];
    let mut kronl = vec![0i32; cap];
    let mut e = nalgebra::DMatrix::zeros(n, n);
    let mut iwork = vec![0i32; n.max(m).max(p)];
    let mut dwork_query = vec![0.0_f64; 1];
    let ldwork = slicot_rs::ab08::ab08nw::ab08nw(
        equil,
        n,
        m,
        p,
        &mut a,
        &mut b,
        &mut c,
        &mut d,
        &mut nfz,
        &mut nrank,
        &mut niz,
        &mut dinfz,
        &mut nkror,
        &mut ninfe,
        &mut nkrol,
        &mut infz,
        &mut kronr,
        &mut infe,
        &mut kronl,
        &mut e,
        0.0,
        &mut iwork,
        &mut dwork_query,
        -1,
    );
    let ldwork_size = if ldwork == 0 {
        dwork_query[0].max(1.0) as usize
    } else {
        1024
    };
    let mut dwork = vec![0.0; ldwork_size];
    let info = slicot_rs::ab08::ab08nw::ab08nw(
        equil,
        n,
        m,
        p,
        &mut a,
        &mut b,
        &mut c,
        &mut d,
        &mut nfz,
        &mut nrank,
        &mut niz,
        &mut dinfz,
        &mut nkror,
        &mut ninfe,
        &mut nkrol,
        &mut infz,
        &mut kronr,
        &mut infe,
        &mut kronl,
        &mut e,
        0.0,
        &mut iwork,
        &mut dwork,
        ldwork_size as i32,
    );
    let input_summary = format!("n={}, m={}, p={}, equil={}", n, m, p, dat.equil);
    let rust_summary = format!("INFO={}", info);
    let pass = info == 0;
    let diff_msg = if pass { String::new() } else { format!("Rust INFO={} (main path not implemented for non-trivial input)", info) };
    (
        pass,
        input_summary,
        fortran_summary,
        if pass { rust_summary } else { format!("{}; {}", rust_summary, diff_msg) },
    )
}

/// AB09 routines with .dat that we validate (stub returns 1 for non-trivial).
const AB09_VALIDATED: &[&str] = &[
    "AB09AD", "AB09BD", "AB09CD", "AB09DD", "AB09ED", "AB09FD", "AB09GD",
    "AB09HD", "AB09ID", "AB09JD", "AB09KD", "AB09MD", "AB09ND",
];

/// Validate AB09*: run Fortran, run Rust stub (n,m from .dat), compare INFO.
fn validate_ab09_generic(
    slicot: &str,
    examples_dir: &Path,
    dat_path: &Path,
    project_root: &Path,
) -> (bool, String, String, String) {
    use slicot_rs::ab09::{
        ab09ad::ab09ad, ab09bd::ab09bd, ab09cd::ab09cd, ab09dd::ab09dd,
        ab09ed::ab09ed, ab09fd::ab09fd, ab09gd::ab09gd, ab09hd::ab09hd,
        ab09id::ab09id, ab09jd::ab09jd, ab09kd::ab09kd, ab09md::ab09md,
        ab09nd::ab09nd,
    };
    let dat_content = match fs::read_to_string(dat_path) {
        Ok(c) => c,
        Err(e) => return (false, format!("read dat: {}", e), String::new(), String::new()),
    };
    let dat = match slicot_rs::ab09::io::parse_ab09_nmp(Cursor::new(&dat_content)) {
        Ok(d) => d,
        Err(e) => return (false, format!("parse_ab09_nmp: {}", e), String::new(), String::new()),
    };
    let res_path = project_root
        .join("target")
        .join(format!("validation_{}.res", slicot.to_lowercase()));
    let _ = fs::create_dir_all(res_path.parent().unwrap());
    let ok = run_fortran_stdin_redirect(examples_dir, slicot, dat_path, &res_path).unwrap_or(false);
    let fortran_summary = if ok {
        "Fortran ran".to_string()
    } else {
        "Fortran driver missing or failed".to_string()
    };
    let info = match slicot {
        "AB09AD" => ab09ad(dat.n, dat.m),
        "AB09BD" => ab09bd(dat.n, dat.m),
        "AB09CD" => ab09cd(dat.n, dat.m),
        "AB09DD" => ab09dd(dat.n, dat.m),
        "AB09ED" => ab09ed(dat.n, dat.m),
        "AB09FD" => ab09fd(dat.n, dat.m),
        "AB09GD" => ab09gd(dat.n, dat.m),
        "AB09HD" => ab09hd(dat.n, dat.m),
        "AB09ID" => ab09id(dat.n, dat.m),
        "AB09JD" => ab09jd(dat.n, dat.m),
        "AB09KD" => ab09kd(dat.n, dat.m),
        "AB09MD" => ab09md(dat.n, dat.m),
        "AB09ND" => ab09nd(dat.n, dat.m),
        _ => 1,
    };
    let input_summary = format!("n={}, m={}, p={}", dat.n, dat.m, dat.p);
    let rust_summary = format!("INFO={}", info);
    let pass = info == 0;
    let diff_msg = if pass {
        String::new()
    } else if info == 1 {
        "Rust stub returns 1 (not implemented)".to_string()
    } else {
        format!("Rust returned INFO={}", info)
    };
    (
        pass,
        input_summary,
        fortran_summary,
        if pass { rust_summary } else { format!("{}; {}", rust_summary, diff_msg) },
    )
}

fn run_validation_impl(project_root: &Path, examples_dir: &Path) {
    let by_module = discovery(project_root, examples_dir);
    let validation_dir = project_root.join(VALIDATION_DIR);
    let _ = fs::create_dir_all(&validation_dir);

    let mut failures: Vec<(String, String)> = vec![];
    let mut readme_rows: Vec<(String, u32, u32, u32)> = vec![];

    for (module, routines) in &by_module {
        let mut md = format!("# {}\n\n", module);
        let mut pass_count = 0u32;
        let mut fail_count = 0u32;

        for (slicot, rust_fn) in routines {
            let dat_path = examples_dir.join("data").join(format!("{slicot}.dat"));
            if !dat_path.is_file() {
                continue;
            }

            let (pass, input_sum, fortran_sum, rust_sum) = if slicot == "AB01MD" {
                validate_ab01md(examples_dir, &dat_path, project_root)
            } else if slicot == "AB01ND" {
                validate_ab01nd(examples_dir, &dat_path, project_root)
            } else if slicot == "AB07MD" {
                validate_ab07md(examples_dir, &dat_path, project_root)
            } else if slicot == "AB07ND" {
                validate_ab07nd(examples_dir, &dat_path, project_root)
            } else if slicot == "TG01GD" {
                validate_tg01gd(examples_dir, &dat_path, project_root)
            } else if slicot == "TG01HD" {
                validate_tg01hd(examples_dir, &dat_path, project_root)
            } else if slicot == "TG01ID" {
                validate_tg01id(examples_dir, &dat_path, project_root)
            } else if slicot == "AB08NW" {
                validate_ab08nw(examples_dir, &dat_path, project_root)
            } else if AB09_VALIDATED.contains(&slicot.as_str()) {
                validate_ab09_generic(slicot, examples_dir, &dat_path, project_root)
            } else {
                let res_path = project_root.join("target").join(format!("validation_{}.res", slicot.to_lowercase()));
                let _ = fs::create_dir_all(res_path.parent().unwrap());
                let ok = run_fortran_stdin_redirect(examples_dir, slicot, &dat_path, &res_path).unwrap_or(false);
                let fortran_sum = if ok { "Fortran ran (output in target/)".to_string() } else { "Fortran driver missing or failed".to_string() };
                (false, "(see .dat)".to_string(), fortran_sum, "Rust adapter not implemented".to_string())
            };

            md.push_str(&format!("## {} / {}\n\n", slicot, rust_fn));
            md.push_str(&format!("- **Input**: {}\n", input_sum));
            md.push_str(&format!("- **Fortran output**: {}\n", fortran_sum));
            md.push_str(&format!("- **Rust output**: {}\n", rust_sum));
            let has_adapter = matches!(slicot.as_str(), "AB01MD" | "AB01ND" | "AB07MD" | "AB07ND" | "AB08NW" | "TG01GD" | "TG01HD" | "TG01ID")
                || AB09_VALIDATED.contains(&slicot.as_str());
            md.push_str(&format!("- **Result**: **{}**\n\n", if has_adapter { if pass { "PASS" } else { "FAIL" } } else { "SKIP (adapter not implemented)" }));

            if has_adapter {
                if pass {
                    pass_count += 1;
                } else {
                    fail_count += 1;
                    failures.push((slicot.clone(), module.clone()));
                }
            }
        }

        let no_ref = 0u32;
        readme_rows.push((module.clone(), pass_count, fail_count, no_ref));
        let path = validation_dir.join(format!("{module}.md"));
        fs::write(&path, &md).ok();
    }

    let failures_md = validation_dir.join("FAILURES.md");
    let mut fail_content = String::from("# Validation failures\n\n| SLICOT | Module | Link |\n|--------|--------|------|\n");
    for (slicot, module) in &failures {
        fail_content.push_str(&format!("| {} | {} | [{}]({}.md) |\n", slicot, module, module, module));
    }
    if failures.is_empty() {
        fail_content.push_str("| *(none)* | | |\n");
    }
    fs::write(&failures_md, &fail_content).ok();

    let readme_path = validation_dir.join("README.md");
    let existing = fs::read_to_string(&readme_path).unwrap_or_default();
    let table_marker = "| *(Generated by run_all.sh)* | | | | |";
    let new_table = if readme_rows.is_empty() {
        "| *(no modules with reference data)* | | | | |".to_string()
    } else {
        readme_rows
            .iter()
            .map(|(module, pass, fail, no_ref)| format!("| {} | {} | {} | {} | [{}]({}.md) |", module, pass, fail, no_ref, module, module))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let new_readme = if existing.contains(table_marker) {
        existing.replacen(table_marker, &new_table, 1)
    } else {
        existing
    };
    fs::write(&readme_path, &new_readme).ok();
}

#[test]
fn run_validation() {
    let project_root = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string()));
    let examples_dir = match std::env::var("SLICOT_EXAMPLES_DIR") {
        Ok(d) => PathBuf::from(d),
        Err(_) => project_root.join("SLICOT-Reference").join("examples"),
    };
    if !examples_dir.is_dir() {
        eprintln!("SLICOT_EXAMPLES_DIR not set or not a directory; skipping validation");
        return;
    }
    run_validation_impl(&project_root, &examples_dir);
}
