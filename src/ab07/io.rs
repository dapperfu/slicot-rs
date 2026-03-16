//! I/O parsers for AB07 .dat and Fortran driver output.

use nalgebra::DMatrix;
use std::io::Read;

/// AB07MD .dat: title, then N M P JOBD, then A(N×N), B(N×M), C(P×N), D(P×M).
pub fn parse_ab07md_dat<R: Read>(mut r: R) -> Result<Ab07mdDat, String> {
    let mut s = String::new();
    r.read_to_string(&mut s).map_err(|e| e.to_string())?;
    let lines: Vec<&str> = s.lines().map(str::trim).collect();
    if lines.len() < 2 {
        return Err("ab07md dat: too few lines".to_string());
    }
    let toks: Vec<&str> = lines[1].split_whitespace().collect();
    if toks.len() < 4 {
        return Err("ab07md dat: need N M P JOBD".to_string());
    }
    let n: usize = toks[0].parse().map_err(|_| "ab07md dat: N")?;
    let m: usize = toks[1].parse().map_err(|_| "ab07md dat: M")?;
    let p: usize = toks[2].parse().map_err(|_| "ab07md dat: P")?;
    let jobd = toks[3] == "D" || toks[3] == "d";

    let mut vals: Vec<f64> = vec![];
    for line in lines.iter().skip(2) {
        for w in line.split_whitespace() {
            if let Ok(v) = w.parse::<f64>() {
                vals.push(v);
            }
        }
    }
    let need = n * n + n * m + p * n + (if jobd { p * m } else { 0 });
    if vals.len() < need {
        return Err(format!("ab07md dat: need {} values, got {}", need, vals.len()));
    }
    let mut idx = 0;
    let mut a = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = vals[idx];
            idx += 1;
        }
    }
    let mut b = DMatrix::zeros(n, m);
    for i in 0..n {
        for j in 0..m {
            b[(i, j)] = vals[idx];
            idx += 1;
        }
    }
    let mut c = DMatrix::zeros(p, n);
    for i in 0..p {
        for j in 0..n {
            c[(i, j)] = vals[idx];
            idx += 1;
        }
    }
    let d = if jobd {
        let mut dmat = DMatrix::zeros(p, m);
        for i in 0..p {
            for j in 0..m {
                dmat[(i, j)] = vals[idx];
                idx += 1;
            }
        }
        Some(dmat)
    } else {
        None
    };
    Ok(Ab07mdDat { n, m, p, jobd, a, b, c, d })
}

pub struct Ab07mdDat {
    pub n: usize,
    pub m: usize,
    pub p: usize,
    pub jobd: bool,
    pub a: DMatrix<f64>,
    pub b: DMatrix<f64>,
    pub c: DMatrix<f64>,
    pub d: Option<DMatrix<f64>>,
}

/// Parse AB07MD Fortran output given dimensions: A (n×n), B (n×p), C (m×n), D (m×p) if jobd.
pub fn parse_ab07md_res<R: Read>(mut r: R, n: usize, m: usize, p: usize, jobd: bool) -> Result<Ab07mdRes, String> {
    let mut s = String::new();
    r.read_to_string(&mut s).map_err(|e| e.to_string())?;
    let nums: Vec<f64> = s
        .lines()
        .flat_map(|l| l.split_whitespace().filter_map(|w| w.parse().ok()))
        .collect();
    let need = n * n + n * p + m * n + if jobd { m * p } else { 0 };
    if nums.len() < need {
        return Err(format!("ab07md res: need {} values, got {}", need, nums.len()));
    }
    let mut idx = 0;
    let mut a = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = nums[idx];
            idx += 1;
        }
    }
    let mut b = DMatrix::zeros(n, p);
    for i in 0..n {
        for j in 0..p {
            b[(i, j)] = nums[idx];
            idx += 1;
        }
    }
    let mut c = DMatrix::zeros(m, n);
    for i in 0..m {
        for j in 0..n {
            c[(i, j)] = nums[idx];
            idx += 1;
        }
    }
    let d = if jobd {
        let mut dmat = DMatrix::zeros(m, p);
        for i in 0..m {
            for j in 0..p {
                dmat[(i, j)] = nums[idx];
                idx += 1;
            }
        }
        Some(dmat)
    } else {
        None
    };
    Ok(Ab07mdRes { a, b, c, d })
}

pub struct Ab07mdRes {
    pub a: DMatrix<f64>,
    pub b: DMatrix<f64>,
    pub c: DMatrix<f64>,
    pub d: Option<DMatrix<f64>>,
}

/// AB07ND .dat: title, then N M, then A(N×N), B(N×M), C(M×N), D(M×M).
pub fn parse_ab07nd_dat<R: Read>(mut r: R) -> Result<Ab07ndDat, String> {
    let mut s = String::new();
    r.read_to_string(&mut s).map_err(|e| e.to_string())?;
    let lines: Vec<&str> = s.lines().map(str::trim).collect();
    if lines.len() < 2 {
        return Err("ab07nd dat: too few lines".to_string());
    }
    let nums: Vec<usize> = lines[1]
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    if nums.len() < 2 {
        return Err("ab07nd dat: need N M".to_string());
    }
    let n = nums[0];
    let m = nums[1];

    let mut vals: Vec<f64> = vec![];
    for line in lines.iter().skip(2) {
        for w in line.split_whitespace() {
            if let Ok(v) = w.parse::<f64>() {
                vals.push(v);
            }
        }
    }
    let need = n * n + n * m + m * n + m * m;
    if vals.len() < need {
        return Err(format!("ab07nd dat: need {} values, got {}", need, vals.len()));
    }
    let mut idx = 0;
    let mut a = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = vals[idx];
            idx += 1;
        }
    }
    let mut b = DMatrix::zeros(n, m);
    for i in 0..n {
        for j in 0..m {
            b[(i, j)] = vals[idx];
            idx += 1;
        }
    }
    let mut c = DMatrix::zeros(m, n);
    for i in 0..m {
        for j in 0..n {
            c[(i, j)] = vals[idx];
            idx += 1;
        }
    }
    let mut d = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            d[(i, j)] = vals[idx];
            idx += 1;
        }
    }
    Ok(Ab07ndDat { n, m, a, b, c, d })
}

pub struct Ab07ndDat {
    pub n: usize,
    pub m: usize,
    pub a: DMatrix<f64>,
    pub b: DMatrix<f64>,
    pub c: DMatrix<f64>,
    pub d: DMatrix<f64>,
}

/// Parse AB07ND Fortran output: A(N×N), B(N×M), C(M×N), D(M×M).
pub fn parse_ab07nd_res<R: Read>(mut r: R, n: usize, m: usize) -> Result<Ab07ndRes, String> {
    let mut s = String::new();
    r.read_to_string(&mut s).map_err(|e| e.to_string())?;
    let nums: Vec<f64> = s
        .lines()
        .flat_map(|l| l.split_whitespace().filter_map(|w| w.parse().ok()))
        .collect();
    let need = n * n + n * m + m * n + m * m;
    if nums.len() < need {
        return Err(format!("ab07nd res: need {} values, got {}", need, nums.len()));
    }
    let mut idx = 0;
    let mut a = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            a[(i, j)] = nums[idx];
            idx += 1;
        }
    }
    let mut b = DMatrix::zeros(n, m);
    for i in 0..n {
        for j in 0..m {
            b[(i, j)] = nums[idx];
            idx += 1;
        }
    }
    let mut c = DMatrix::zeros(m, n);
    for i in 0..m {
        for j in 0..n {
            c[(i, j)] = nums[idx];
            idx += 1;
        }
    }
    let mut d = DMatrix::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            d[(i, j)] = nums[idx];
            idx += 1;
        }
    }
    Ok(Ab07ndRes { a, b, c, d })
}

pub struct Ab07ndRes {
    pub a: DMatrix<f64>,
    pub b: DMatrix<f64>,
    pub c: DMatrix<f64>,
    pub d: DMatrix<f64>,
}
