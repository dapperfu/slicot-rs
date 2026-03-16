//! Benchmark runner for common control routines: runs each (routine, size) for up to 30s
//! and prints CSV (routine,n,time_us) to stdout for plotting.

use std::hint::black_box;
use std::time::{Duration, Instant};

use nalgebra::{DMatrix, DVector};
use slicot_rs::de01::de01od::{de01od, De01OdConv};
use slicot_rs::dlac::dlacpy_slc::{dlacpy_slc, DlacpyUplo};
use slicot_rs::ma02::ma02ed::{ma02ed, Ma02EdUplo};
use slicot_rs::ma02::ma02es::{ma02es, Ma02EsUplo};
use slicot_rs::mb01::mb01md::{mb01md, Mb01MdUplo};
use slicot_rs::tb01::tb01md::{tb01md, JobU, Uplo};

/// Max wall-clock time per (routine, size).
const BUDGET_PER_SIZE: Duration = Duration::from_secs(30);

/// Size ladder: up through 1024 (stay under 30s per size for typical routines).
fn size_ladder() -> Vec<usize> {
    vec![32, 64, 128, 256, 512, 1024]
}

fn matrix_nn(n: usize) -> DMatrix<f64> {
    DMatrix::from_fn(n, n, |i, j| (i + j) as f64 * 0.1)
}

fn matrix_nm(n: usize, m: usize) -> DMatrix<f64> {
    DMatrix::from_fn(n, m, |i, j| (i * 2 + j) as f64 * 0.1)
}

fn vector_n(n: usize) -> DVector<f64> {
    DVector::from_fn(n, |i, _| (i + 1) as f64 * 0.1)
}

fn vec_n(n: usize) -> Vec<f64> {
    (0..n).map(|i| (i + 1) as f64 * 0.1).collect()
}

fn m_from_n(n: usize) -> usize {
    (n / 2).max(1)
}

/// Run one routine at one size for up to BUDGET_PER_SIZE; return mean time per call in microseconds.
fn run_ma02ed(n: usize) -> f64 {
    let a = matrix_nn(n);
    let start = Instant::now();
    let mut count = 0u64;
    while start.elapsed() < BUDGET_PER_SIZE {
        let mut a = a.clone();
        let _ = ma02ed(black_box(Ma02EdUplo::Upper), &mut a);
        count += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    if count == 0 {
        elapsed * 1e6
    } else {
        (elapsed * 1e6) / (count as f64)
    }
}

fn run_ma02es(n: usize) -> f64 {
    let a = matrix_nn(n);
    let start = Instant::now();
    let mut count = 0u64;
    while start.elapsed() < BUDGET_PER_SIZE {
        let mut a = a.clone();
        let _ = ma02es(black_box(Ma02EsUplo::Upper), &mut a);
        count += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    if count == 0 {
        elapsed * 1e6
    } else {
        (elapsed * 1e6) / (count as f64)
    }
}

fn run_tb01md(n: usize) -> f64 {
    let m = m_from_n(n);
    let a = matrix_nn(n);
    let b = matrix_nm(n, m);
    let start = Instant::now();
    let mut count = 0u64;
    while start.elapsed() < BUDGET_PER_SIZE {
        let mut a = a.clone();
        let mut b = b.clone();
        let mut u = None;
        let _ = tb01md(black_box(JobU::No), black_box(Uplo::Upper), &mut a, &mut b, &mut u);
        count += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    if count == 0 {
        elapsed * 1e6
    } else {
        (elapsed * 1e6) / (count as f64)
    }
}

fn run_dlacpy_slc(n: usize) -> f64 {
    let a = matrix_nn(n);
    let mut b = DMatrix::zeros(n, n);
    let start = Instant::now();
    let mut count = 0u64;
    while start.elapsed() < BUDGET_PER_SIZE {
        b.fill(0.0);
        let _ = dlacpy_slc(black_box(DlacpyUplo::All), &a, &mut b);
        count += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    if count == 0 {
        elapsed * 1e6
    } else {
        (elapsed * 1e6) / (count as f64)
    }
}

fn run_mb01md(n: usize) -> f64 {
    let a = matrix_nn(n);
    let x = vector_n(n);
    let mut y = DVector::zeros(n);
    let start = Instant::now();
    let mut count = 0u64;
    while start.elapsed() < BUDGET_PER_SIZE {
        y.fill(0.0);
        let _ = mb01md(
            black_box(Mb01MdUplo::Upper),
            black_box(1.0),
            &a,
            &x,
            black_box(0.0),
            &mut y,
        );
        count += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    if count == 0 {
        elapsed * 1e6
    } else {
        (elapsed * 1e6) / (count as f64)
    }
}

fn run_de01od(n: usize) -> f64 {
    let a = vec_n(n);
    let b = vec_n(n);
    let start = Instant::now();
    let mut count = 0u64;
    while start.elapsed() < BUDGET_PER_SIZE {
        let mut a = a.clone();
        let mut b = b.clone();
        let _ = de01od(black_box(De01OdConv::Convolution), n, &mut a, &mut b);
        count += 1;
    }
    let elapsed = start.elapsed().as_secs_f64();
    if count == 0 {
        elapsed * 1e6
    } else {
        (elapsed * 1e6) / (count as f64)
    }
}

fn main() {
    let sizes = size_ladder();
    println!("routine,n,time_us");

    for &n in &sizes {
        if n < 2 {
            continue;
        }
        let t = run_ma02ed(n);
        println!("MA02ED,{},{}", n, t);
    }
    for &n in &sizes {
        if n < 2 {
            continue;
        }
        let t = run_ma02es(n);
        println!("MA02ES,{},{}", n, t);
    }
    for &n in &sizes {
        let t = run_tb01md(n);
        println!("TB01MD,{},{}", n, t);
    }
    for &n in &sizes {
        let t = run_dlacpy_slc(n);
        println!("DLACPY_SLC,{},{}", n, t);
    }
    for &n in &sizes {
        let t = run_mb01md(n);
        println!("MB01MD,{},{}", n, t);
    }
    for &n in &sizes {
        if !n.is_power_of_two() {
            continue;
        }
        let t = run_de01od(n);
        println!("DE01OD,{},{}", n, t);
    }
}
