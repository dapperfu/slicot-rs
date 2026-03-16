//! SLICOT-style .dat / .res I/O for fuzzer and tests.
//!
//! Pilot: AB01ND, AB01MD. Layout matches SLICOT-Reference examples/data and results.

use nalgebra::{DMatrix, DVector};
use std::io::{BufRead, BufReader, Read, Write};

/// AB01ND input: N, M, TOL, JOBZ, A (N×N), B (N×M).
#[derive(Clone, Debug)]
pub struct Ab01ndDat {
    pub n: usize,
    pub m: usize,
    pub tol: f64,
    pub jobz: char,
    pub a: DMatrix<f64>,
    pub b: DMatrix<f64>,
}

/// AB01ND output: INFO, NCONT, INDCON, A block, B block, Z (optional).
#[derive(Clone, Debug, Default)]
pub struct Ab01ndRes {
    pub info: i32,
    pub ncont: usize,
    pub indcon: usize,
    pub nblk: Vec<i32>,
    pub a_cont: DMatrix<f64>,
    pub b_cont: DMatrix<f64>,
    pub z: Option<DMatrix<f64>>,
}

/// Parse AB01ND .dat from a reader (optional title line, then N M TOL JOBZ, then A, then B).
pub fn parse_ab01nd_dat<R: Read>(r: R) -> Result<Ab01ndDat, String> {
    let mut lines = BufReader::new(r).lines();
    let mut next = |msg: &str| {
        lines
            .next()
            .ok_or_else(|| msg.to_string())?
            .map_err(|e| e.to_string())
    };
    // Optional title
    let first = next("missing first line")?;
    let first = first.trim();
    let (n, m, tol, jobz) = if first.starts_with("AB01ND") || first.starts_with(' ') {
        let line = next("missing N M TOL JOBZ line")?;
        parse_n_m_tol_jobz(&line)?
    } else {
        parse_n_m_tol_jobz(first)?
    };
    let mut a = vec![];
    for _ in 0..n {
        let line = next("missing A row")?;
        for w in line.split_whitespace() {
            a.push(w.parse::<f64>().map_err(|_| "bad A value")?);
        }
    }
    // B is stored in file as M lines of N values each (Fortran: (B(I,J), I=1,N), J=1,M)
    let mut b = vec![];
    for _ in 0..m {
        let line = next("missing B row")?;
        for w in line.split_whitespace() {
            b.push(w.parse::<f64>().map_err(|_| "bad B value")?);
        }
    }
    let a = DMatrix::from_row_slice(n, n, &a);
    // b is column-major from file: columns 1..M, each with N rows
    let b = DMatrix::from_column_slice(n, m, &b);
    Ok(Ab01ndDat { n, m, tol, jobz, a, b })
}

fn parse_n_m_tol_jobz(line: &str) -> Result<(usize, usize, f64, char), String> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 4 {
        return Err("need N M TOL JOBZ".to_string());
    }
    let n = parts[0].parse::<usize>().map_err(|_| "bad N")?;
    let m = parts[1].parse::<usize>().map_err(|_| "bad M")?;
    let tol = parts[2].parse::<f64>().map_err(|_| "bad TOL")?;
    let jobz = parts[3].chars().next().ok_or("missing JOBZ")?;
    Ok((n, m, tol, jobz))
}

/// Write AB01ND .dat to a writer.
pub fn write_ab01nd_dat<W: Write>(w: &mut W, d: &Ab01ndDat) -> std::io::Result<()> {
    writeln!(w, " AB01ND EXAMPLE PROGRAM DATA")?;
    writeln!(w, "   {}     {}     {}     {}", d.n, d.m, d.tol, d.jobz)?;
    for i in 0..d.n {
        for j in 0..d.n {
            write!(w, " {}", d.a[(i, j)])?;
        }
        writeln!(w)?;
    }
    for j in 0..d.m {
        for i in 0..d.n {
            write!(w, " {}", d.b[(i, j)])?;
        }
        writeln!(w)?;
    }
    Ok(())
}

/// Parse AB01ND .res from a reader.
pub fn parse_ab01nd_res<R: Read>(r: R) -> Result<Ab01ndRes, String> {
    let lines: Vec<String> = BufReader::new(r).lines().collect::<Result<_, _>>().map_err(|e| e.to_string())?;
    let mut res = Ab01ndRes::default();
    let mut i = 0;
    while i < lines.len() {
        let s = lines[i].trim();
        if s.contains("INFO on exit") {
            let num: i32 = s.split('=').nth(1).and_then(|x| x.trim().parse().ok()).unwrap_or(0);
            res.info = num;
            i += 1;
            continue;
        }
        if s.contains("order of the controllable") {
            res.ncont = s.split('=').nth(1).and_then(|x| x.trim().parse().ok()).unwrap_or(0);
            i += 1;
            i += 1; // skip blank / "The transformed state dynamics"
            if i >= lines.len() {
                break;
            }
            let mut a_rows = vec![];
            while i < lines.len() && !lines[i].trim().is_empty() && !lines[i].contains("dimensions of its diagonal") {
                let row: Vec<f64> = lines[i].split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if !row.is_empty() {
                    a_rows.push(row);
                }
                i += 1;
            }
            if !a_rows.is_empty() {
                let cols = a_rows[0].len();
                let flat: Vec<f64> = a_rows.into_iter().flat_map(|r| r.into_iter()).collect();
                res.a_cont = DMatrix::from_row_slice(res.ncont, cols, &flat);
            }
            continue;
        }
        if s.contains("dimensions of its diagonal") {
            i += 1;
            if i < lines.len() {
                res.nblk = lines[i].split_whitespace().filter_map(|x| x.parse().ok()).collect();
                res.indcon = res.nblk.len();
            }
            i += 1;
            i += 1; // skip blank / "The transformed input/state"
            if i >= lines.len() {
                break;
            }
            let mut b_rows = vec![];
            while i < lines.len() && !lines[i].trim().is_empty() && !lines[i].contains("controllability index") {
                let row: Vec<f64> = lines[i].split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if !row.is_empty() {
                    b_rows.push(row);
                }
                i += 1;
            }
            if !b_rows.is_empty() {
                let cols = b_rows[0].len();
                let flat: Vec<f64> = b_rows.into_iter().flat_map(|r| r.into_iter()).collect();
                res.b_cont = DMatrix::from_row_slice(res.ncont, cols, &flat);
            }
            continue;
        }
        if s.contains("controllability index of the transformed") {
            res.indcon = s.split('=').nth(1).and_then(|x| x.trim().parse().ok()).unwrap_or(res.indcon);
            i += 1;
            i += 1;
            if i < lines.len() && lines[i].contains("similarity transformation matrix Z") {
                i += 1;
                let n = res.a_cont.nrows();
                let mut z_rows = vec![];
                for _ in 0..n {
                    i += 1;
                    if i < lines.len() {
                        let row: Vec<f64> = lines[i].split_whitespace().filter_map(|x| x.parse().ok()).collect();
                        if !row.is_empty() {
                            z_rows.push(row);
                        }
                    }
                }
                if !z_rows.is_empty() {
                    let nz = z_rows.len();
                    let cols = z_rows[0].len();
                    let flat: Vec<f64> = z_rows.into_iter().flat_map(|r| r.into_iter()).collect();
                    res.z = Some(DMatrix::from_row_slice(nz, cols, &flat));
                }
            }
            break;
        }
        i += 1;
    }
    Ok(res)
}

/// Write AB01ND .res to a writer (same style as Fortran output).
pub fn write_ab01nd_res<W: Write>(w: &mut W, res: &Ab01ndRes) -> std::io::Result<()> {
    writeln!(w, " AB01ND EXAMPLE PROGRAM RESULTS")?;
    writeln!(w)?;
    if res.info != 0 {
        writeln!(w, " INFO on exit from AB01ND = {}", res.info)?;
        return Ok(());
    }
    writeln!(w, " The order of the controllable state-space representation = {}", res.ncont)?;
    writeln!(w)?;
    writeln!(w, " The transformed state dynamics matrix of a controllable realization is ")?;
    for i in 0..res.a_cont.nrows() {
        for j in 0..res.a_cont.ncols() {
            write!(w, " {:8.4}", res.a_cont[(i, j)])?;
        }
        writeln!(w)?;
    }
    writeln!(w)?;
    writeln!(w, " and the dimensions of its diagonal blocks are ")?;
    for k in &res.nblk {
        write!(w, " {:2}", k)?;
    }
    writeln!(w)?;
    writeln!(w)?;
    writeln!(w, " The transformed input/state matrix B of a controllable realization is ")?;
    for i in 0..res.b_cont.nrows() {
        for j in 0..res.b_cont.ncols() {
            write!(w, " {:8.4}", res.b_cont[(i, j)])?;
        }
        writeln!(w)?;
    }
    writeln!(w)?;
    writeln!(w, " The controllability index of the transformed system representation = {}", res.indcon)?;
    if let Some(ref z) = res.z {
        writeln!(w)?;
        writeln!(w, " The similarity transformation matrix Z is ")?;
        for i in 0..z.nrows() {
            for j in 0..z.ncols() {
                write!(w, " {:8.4}", z[(i, j)])?;
            }
            writeln!(w)?;
        }
    }
    Ok(())
}

/// AB01MD input: N, TOL, JOBZ, A (N×N), B (N×1).
#[derive(Clone, Debug)]
pub struct Ab01mdDat {
    pub n: usize,
    pub tol: f64,
    pub jobz: char,
    pub a: DMatrix<f64>,
    pub b: DVector<f64>,
}

/// AB01MD output: NCONT, A (ncont×ncont), B (ncont×1), Z (N×N optional).
#[derive(Clone, Debug, Default)]
pub struct Ab01mdRes {
    pub ncont: usize,
    pub a_cont: DMatrix<f64>,
    pub b_cont: DMatrix<f64>,
    pub z: Option<DMatrix<f64>>,
}

/// Parse AB01MD .dat: title, then N TOL JOBZ, then A (N×N), then B (N×1).
pub fn parse_ab01md_dat<R: Read>(r: R) -> Result<Ab01mdDat, String> {
    let mut lines = BufReader::new(r).lines();
    let mut next = |msg: &str| {
        lines
            .next()
            .ok_or_else(|| msg.to_string())?
            .map_err(|e| e.to_string())
    };
    let first = next("missing first line")?;
    let first = first.trim();
    let (n, tol, jobz) = if first.starts_with("AB01MD") || first.starts_with(' ') {
        let line = next("missing N TOL JOBZ line")?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            return Err("need N TOL JOBZ".to_string());
        }
        let n = parts[0].parse::<usize>().map_err(|_| "bad N")?;
        let tol = parts[1].parse::<f64>().map_err(|_| "bad TOL")?;
        let jobz = parts[2].chars().next().ok_or("missing JOBZ")?;
        (n, tol, jobz)
    } else {
        let parts: Vec<&str> = first.split_whitespace().collect();
        if parts.len() < 3 {
            return Err("need N TOL JOBZ".to_string());
        }
        let n = parts[0].parse::<usize>().map_err(|_| "bad N")?;
        let tol = parts[1].parse::<f64>().map_err(|_| "bad TOL")?;
        let jobz = parts[2].chars().next().ok_or("missing JOBZ")?;
        (n, tol, jobz)
    };
    let mut a = vec![];
    for _ in 0..n {
        let line = next("missing A row")?;
        for w in line.split_whitespace() {
            a.push(w.parse::<f64>().map_err(|_| "bad A value")?);
        }
    }
    let mut b = vec![];
    while b.len() < n {
        let line = next("missing B row")?;
        for w in line.split_whitespace() {
            b.push(w.parse::<f64>().map_err(|_| "bad B value")?);
            if b.len() == n {
                break;
            }
        }
    }
    let a = DMatrix::from_row_slice(n, n, &a);
    let b = DVector::from_row_slice(&b);
    Ok(Ab01mdDat { n, tol, jobz, a, b })
}

/// Parse AB01MD .res: order of controllable, A matrix, B vector, Z matrix.
pub fn parse_ab01md_res<R: Read>(r: R) -> Result<Ab01mdRes, String> {
    let lines: Vec<String> = BufReader::new(r).lines().collect::<Result<_, _>>().map_err(|e| e.to_string())?;
    let mut res = Ab01mdRes::default();
    let mut i = 0;
    while i < lines.len() {
        let s = lines[i].trim();
        if s.contains("order of the controllable") {
            res.ncont = s.split('=').nth(1).and_then(|x| x.trim().parse().ok()).unwrap_or(0);
            i += 1;
            i += 1; // skip blank / "The state dynamics matrix"
            if i >= lines.len() {
                break;
            }
            let mut a_rows = vec![];
            while i < lines.len() && !lines[i].trim().is_empty() && !lines[i].contains("input/state vector") {
                let row: Vec<f64> = lines[i].split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if !row.is_empty() {
                    a_rows.push(row);
                }
                i += 1;
            }
            if !a_rows.is_empty() {
                let cols = a_rows[0].len();
                let flat: Vec<f64> = a_rows.into_iter().flat_map(|r| r.into_iter()).collect();
                res.a_cont = DMatrix::from_row_slice(res.ncont, cols, &flat);
            }
            continue;
        }
        if s.contains("input/state vector") || s.contains("input/state matrix") {
            i += 1;
            let mut b_rows = vec![];
            while i < lines.len() && !lines[i].trim().is_empty() && !lines[i].contains("similarity transformation") {
                let row: Vec<f64> = lines[i].split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if !row.is_empty() {
                    b_rows.push(row);
                }
                i += 1;
            }
            if !b_rows.is_empty() {
                let cols = b_rows[0].len();
                let flat: Vec<f64> = b_rows.into_iter().flat_map(|r| r.into_iter()).collect();
                res.b_cont = DMatrix::from_row_slice(res.ncont, cols, &flat);
            }
            continue;
        }
        if s.contains("similarity transformation matrix Z") {
            i += 1;
            let n = res.a_cont.nrows();
            let mut z_rows = vec![];
            for _ in 0..n {
                if i < lines.len() {
                    let row: Vec<f64> = lines[i].split_whitespace().filter_map(|x| x.parse().ok()).collect();
                    if !row.is_empty() {
                        z_rows.push(row);
                    }
                    i += 1;
                }
            }
            if !z_rows.is_empty() {
                let nz = z_rows.len();
                let cols = z_rows[0].len();
                let flat: Vec<f64> = z_rows.into_iter().flat_map(|r| r.into_iter()).collect();
                res.z = Some(DMatrix::from_row_slice(nz, cols, &flat));
            }
            break;
        }
        i += 1;
    }
    Ok(res)
}

/// Relative tolerance comparison: |a - b| <= rel_tol * max(|a|, |b|, 1.0).
#[inline]
pub fn rel_tol_eq(a: f64, b: f64, rel_tol: f64) -> bool {
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= rel_tol * scale
}

/// Compare two matrices element-wise with relative tolerance; returns first (row, col) where they differ, or None.
pub fn rel_tol_eq_matrix(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    rel_tol: f64,
) -> Option<(usize, usize)> {
    if a.nrows() != b.nrows() || a.ncols() != b.ncols() {
        return Some((0, 0));
    }
    for i in 0..a.nrows() {
        for j in 0..a.ncols() {
            if !rel_tol_eq(a[(i, j)], b[(i, j)], rel_tol) {
                return Some((i, j));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_ab01nd_dat_example() {
        // 8 lines: title, N M TOL JOBZ, 3×A rows, 3×B rows (trailing newline for 8th line)
        let raw = " AB01ND EXAMPLE PROGRAM DATA\n   3     2     0.0     I\n  -1.0   0.0   0.0\n  -2.0  -2.0  -2.0\n  -1.0   0.0  -3.0\n   1.0   0.0   0.0\n   0.0   2.0   1.0\n";
        let d = parse_ab01nd_dat(raw.as_bytes()).unwrap();
        assert_eq!(d.n, 3);
        assert_eq!(d.m, 2);
        assert_eq!(d.jobz, 'I');
        assert_eq!(d.a[(0, 0)], -1.0);
        assert_eq!(d.b[(0, 0)], 1.0);
    }

    #[test]
    fn roundtrip_ab01nd_dat() {
        let d = Ab01ndDat {
            n: 2,
            m: 1,
            tol: 0.0,
            jobz: 'N',
            a: DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]),
            b: DMatrix::from_row_slice(2, 1, &[0.0, 1.0]),
        };
        let mut out = vec![];
        write_ab01nd_dat(&mut out, &d).unwrap();
        let d2 = parse_ab01nd_dat(out.as_slice()).unwrap();
        assert_eq!(d2.n, d.n);
        assert_eq!(d2.m, d.m);
        assert_eq!(d2.a[(0, 0)], d.a[(0, 0)]);
    }
}
