//! Criterion benchmarks for every SLICOT routine with scaled problem sizes.
//! Each (routine, size) is a separate benchmark; sizes use a shared ladder (no tiny 2×2).

mod common;

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use nalgebra::{DMatrix, DVector};
use slicot_rs::de01::de01od::{de01od, De01OdConv};
use slicot_rs::dlac::dlacpy_slc::{dlacpy_slc, DlacpyUplo};
use slicot_rs::ma02::ma02ed::{ma02ed, Ma02EdUplo};
use slicot_rs::mb01::mb01md::{mb01md, Mb01MdUplo};
use slicot_rs::tb01::tb01md::{tb01md, JobU, Uplo};

use common::*;

fn bench_ma02ed(c: &mut Criterion) {
    let mut group = c.benchmark_group("ma02/ma02ed");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = ma02ed(black_box(Ma02EdUplo::Upper), &mut a);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_mb01md(c: &mut Criterion) {
    let mut group = c.benchmark_group("mb01/mb01md");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let x = vector_n(n);
        let mut y = DVector::zeros(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                y.fill(0.0);
                let info = mb01md(
                    black_box(Mb01MdUplo::Upper),
                    black_box(1.0),
                    &a,
                    &x,
                    black_box(0.0),
                    &mut y,
                );
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_tb01md(c: &mut Criterion) {
    let mut group = c.benchmark_group("tb01/tb01md");
    for &n in SIZE_LADDER_N {
        let m = m_from_n(n);
        let a = matrix_nn(n);
        let b = matrix_nm(n, m);
        group.throughput(Throughput::Elements((n * n + n * m) as u64));
        group.bench_function(format!("n{n}_m{m}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut b = b.clone();
                let info = tb01md(
                    black_box(JobU::No),
                    black_box(Uplo::Upper),
                    &mut a,
                    &mut b,
                    None,
                );
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_dlacpy_slc(c: &mut Criterion) {
    let mut group = c.benchmark_group("dlac/dlacpy_slc");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let b = DMatrix::zeros(n, n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut b = b.clone();
                let info = dlacpy_slc(black_box(DlacpyUplo::All), &a, &mut b);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_de01od(c: &mut Criterion) {
    let mut group = c.benchmark_group("de01/de01od");
    for &n in SIZE_LADDER_POW2 {
        let a = vec_n(n);
        let b = vec_n(n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut b = b.clone();
                let info = de01od(black_box(De01OdConv::Convolution), n, &mut a, &mut b);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

// --- MA02: additional matrix routines ---
fn bench_ma02es(c: &mut Criterion) {
    use slicot_rs::ma02::ma02es::{ma02es, Ma02EsUplo};
    let mut group = c.benchmark_group("ma02/ma02es");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = ma02es(black_box(Ma02EsUplo::Upper), &mut a);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02ad(c: &mut Criterion) {
    use slicot_rs::ma02::ma02ad::{ma02ad, Ma02AdJob};
    let mut group = c.benchmark_group("ma02/ma02ad");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let mut b = DMatrix::zeros(n, n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let info = ma02ad(black_box(Ma02AdJob::All), &a, &mut b);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02az(c: &mut Criterion) {
    use slicot_rs::ma02::ma02az::{ma02az, Ma02AzJob};
    let mut group = c.benchmark_group("ma02/ma02az");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let mut b = DMatrix::zeros(n, n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let info = ma02az(black_box('N'), black_box(Ma02AzJob::All), &a, &mut b);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02bd(c: &mut Criterion) {
    use slicot_rs::ma02::ma02bd::{ma02bd, Ma02BdSide};
    let mut group = c.benchmark_group("ma02/ma02bd");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = ma02bd(black_box(Ma02BdSide::Left), &mut a);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02bz(c: &mut Criterion) {
    use slicot_rs::ma02::ma02bz::{ma02bz, Ma02BzSide};
    let mut group = c.benchmark_group("ma02/ma02bz");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = ma02bz(black_box(Ma02BzSide::Left), &mut a);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02cd(c: &mut Criterion) {
    use slicot_rs::ma02::ma02cd::ma02cd;
    let mut group = c.benchmark_group("ma02/ma02cd");
    for &n in SIZE_LADDER_N {
        if n <= 1 {
            continue;
        }
        let a = matrix_nn(n);
        let kl = n / 4;
        let ku = n / 4;
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = ma02cd(&mut a, kl, ku);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02dd(c: &mut Criterion) {
    use slicot_rs::ma02::ma02dd::{ma02dd, Ma02DdJob, Ma02DdUplo};
    let mut group = c.benchmark_group("ma02/ma02dd");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let ap = vec![0.0; n * (n + 1) / 2];
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut ap = ap.clone();
                let info = ma02dd(
                    black_box(Ma02DdJob::Pack),
                    black_box(Ma02DdUplo::Upper),
                    &mut a,
                    &mut ap,
                );
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02gd(c: &mut Criterion) {
    use slicot_rs::ma02::ma02gd::ma02gd;
    let mut group = c.benchmark_group("ma02/ma02gd");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let ipiv: Vec<i32> = (0..n).map(|i| i as i32).collect();
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = ma02gd(&mut a, 0, n.saturating_sub(1), &ipiv);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02pd(c: &mut Criterion) {
    use slicot_rs::ma02::ma02pd::ma02pd;
    let mut group = c.benchmark_group("ma02/ma02pd");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut nzr = 0usize;
                let mut nzc = 0usize;
                let info = ma02pd(&a, &mut nzr, &mut nzc);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_ma02fd(c: &mut Criterion) {
    use slicot_rs::ma02::ma02fd::ma02fd;
    let mut group = c.benchmark_group("ma02/ma02fd");
    group.bench_function("scalar", |bencher| {
        bencher.iter(|| {
            let mut x1 = 3.0f64;
            let x2 = 1.0f64;
            let mut c = 0.0f64;
            let mut s = 0.0f64;
            let _ = ma02fd(&mut x1, x2, &mut c, &mut s);
        });
    });
    group.finish();
}

fn bench_ma02rd(c: &mut Criterion) {
    use slicot_rs::ma02::ma02rd::ma02rd;
    let mut group = c.benchmark_group("ma02/ma02rd");
    for &n in SIZE_LADDER_N {
        let d = vec_n(n);
        let e = vec_n(n.saturating_sub(1));
        if e.is_empty() {
            continue;
        }
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut d = d.clone();
                let mut e = e.clone();
                let _ = ma02rd(black_box('N'), &mut d, &mut e);
            });
        });
    }
    group.finish();
}

// --- MB01: additional matrix routines ---
fn bench_mb01xd(c: &mut Criterion) {
    use slicot_rs::mb01::mb01xd::{mb01xd, Mb01XdUplo};
    let mut group = c.benchmark_group("mb01/mb01xd");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = mb01xd(black_box(Mb01XdUplo::Upper), &mut a);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_mb01xy(c: &mut Criterion) {
    use slicot_rs::mb01::mb01xy::{mb01xy, Mb01XyUplo};
    let mut group = c.benchmark_group("mb01/mb01xy");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = mb01xy(black_box(Mb01XyUplo::Upper), &mut a);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_mb01ss(c: &mut Criterion) {
    use slicot_rs::mb01::mb01ss::{mb01ss, Mb01SsJobs, Mb01SsUplo};
    let mut group = c.benchmark_group("mb01/mb01ss");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let d = vec_n(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = mb01ss(
                    black_box(Mb01SsJobs::Scale),
                    black_box(Mb01SsUplo::Upper),
                    &mut a,
                    &d,
                );
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_mb01sd(c: &mut Criterion) {
    use slicot_rs::mb01::mb01sd::{mb01sd, Mb01SdJobs};
    let mut group = c.benchmark_group("mb01/mb01sd");
    for &n in SIZE_LADDER_N {
        let a = matrix_nn(n);
        let r = vec_n(n);
        let c = vec_n(n);
        group.throughput(Throughput::Elements((n * n) as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = mb01sd(black_box(Mb01SdJobs::Both), &mut a, &r, &c);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

// --- DLAT ---
fn bench_dlatzm(c: &mut Criterion) {
    use slicot_rs::dlat::dlatzm::{dlatzm, DlatzmSide};
    let mut group = c.benchmark_group("dlat/dlatzm");
    for &n in SIZE_LADDER_N {
        let m = n;
        let v: Vec<f64> = (0..m.saturating_sub(1)).map(|i| (i + 1) as f64 * 0.1).collect();
        let c = matrix_nm(m, n);
        group.throughput(Throughput::Elements((m * n) as u64));
        group.bench_function(format!("m{m}_n{n}"), |bencher| {
            bencher.iter(|| {
                let mut c = c.clone();
                let info = dlatzm(
                    black_box(DlatzmSide::Left),
                    &v,
                    1,
                    black_box(0.5),
                    &mut c,
                );
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

// --- AB01 (stubs: don't assert info==0) ---
fn bench_ab01nd(c: &mut Criterion) {
    use slicot_rs::ab01::ab01nd::{ab01nd, JobZ};
    let mut group = c.benchmark_group("ab01/ab01nd");
    for &n in SIZE_LADDER_N {
        let m = m_from_n(n);
        let a = matrix_nn(n);
        let b = matrix_nm(n, m);
        let mut ncont = 0usize;
        let mut indcon = 0usize;
        let mut nblk = vec![0i32; n + 1];
        group.throughput(Throughput::Elements((n * n + n * m) as u64));
        group.bench_function(format!("n{n}_m{m}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut b = b.clone();
                let _info = ab01nd(
                    black_box(JobZ::No),
                    n,
                    m,
                    &mut a,
                    &mut b,
                    &mut ncont,
                    &mut indcon,
                    &mut nblk,
                    None,
                    1e-10,
                );
            });
        });
    }
    group.finish();
}

fn bench_ab01od(c: &mut Criterion) {
    use slicot_rs::ab01::ab01od::{ab01od, JobUV, Stages};
    let mut group = c.benchmark_group("ab01/ab01od");
    for &n in SIZE_LADDER_N {
        let m = m_from_n(n);
        let a = matrix_nn(n);
        let b = matrix_nm(n, m);
        let mut ncont = 0usize;
        let mut indcon = 0usize;
        let mut kstair = vec![0i32; n + 1];
        group.throughput(Throughput::Elements((n * n + n * m) as u64));
        group.bench_function(format!("n{n}_m{m}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut b = b.clone();
                let _info = ab01od(
                    black_box(Stages::All),
                    black_box(JobUV::No),
                    black_box(JobUV::No),
                    n,
                    m,
                    &mut a,
                    &mut b,
                    None,
                    None,
                    &mut ncont,
                    &mut indcon,
                    &mut kstair,
                    1e-10,
                );
            });
        });
    }
    group.finish();
}

// --- AB08MD (state-space with workspace) ---
fn bench_ab08md(c: &mut Criterion) {
    use slicot_rs::ab08::ab08md::{ab08md, Ab08MdEquil};
    let mut group = c.benchmark_group("ab08/ab08md");
    for &n in SIZE_LADDER_N {
        let m = m_from_n(n);
        let p = p_from_n(n);
        let (a, b, c, d) = state_space_matrices(n, m, p);
        let mut a_vec = Vec::with_capacity(n * n);
        for j in 0..n {
            for i in 0..n {
                a_vec.push(a[(i, j)]);
            }
        }
        let mut b_vec = Vec::with_capacity(n * m);
        for j in 0..m {
            for i in 0..n {
                b_vec.push(b[(i, j)]);
            }
        }
        let mut c_vec = Vec::with_capacity(p * n);
        for j in 0..n {
            for i in 0..p {
                c_vec.push(c[(i, j)]);
            }
        }
        let mut d_vec = Vec::with_capacity(p * m);
        for j in 0..m {
            for i in 0..p {
                d_vec.push(d[(i, j)]);
            }
        }
        let lda = n.max(1);
        let ldb = n.max(1);
        let ldc = p.max(1);
        let ldd = p.max(1);
        let np = n + p;
        let nm = n + m;
        let kw = (np * nm)
            + (p.min(m) + (3 * m).saturating_sub(1).max(n)).max(1).max(
                p.min(n) + (3 * p).saturating_sub(1).max(np).max(nm),
            );
        let mut rank = 0i32;
        let mut iwork = vec![0i32; p.min(m) + n + p + m];
        let mut dwork = vec![0.0; kw];
        group.throughput(Throughput::Elements((n * n + n * m + p * n + p * m) as u64));
        group.bench_function(format!("n{n}_m{m}_p{p}"), |bencher| {
            bencher.iter(|| {
                let info = ab08md(
                    black_box(Ab08MdEquil::No),
                    n,
                    m,
                    p,
                    &a_vec,
                    lda,
                    &b_vec,
                    ldb,
                    &c_vec,
                    ldc,
                    &d_vec,
                    ldd,
                    &mut rank,
                    1e-10,
                    &mut iwork,
                    &mut dwork,
                    kw as i32,
                );
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

// --- MA01AD (scalar) ---
fn bench_ma01ad(c: &mut Criterion) {
    use slicot_rs::ma01::ma01ad::ma01ad;
    let mut group = c.benchmark_group("ma01/ma01ad");
    group.bench_function("scalar", |bencher| {
        bencher.iter(|| {
            let mut yr = 0.0f64;
            let mut yi = 0.0f64;
            let _ = ma01ad(black_box(4.0), black_box(0.0), &mut yr, &mut yi);
        });
    });
    group.finish();
}

// --- DG01 (FFT/signal, power-of-2) ---
fn bench_dg01md(c: &mut Criterion) {
    use slicot_rs::dg01::dg01md::{dg01md, Dg01MdIndi};
    let mut group = c.benchmark_group("dg01/dg01md");
    for &n in SIZE_LADDER_POW2 {
        let xr = vec_n(n);
        let xi = vec_n(n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut xr = xr.clone();
                let mut xi = xi.clone();
                let info = dg01md(black_box(Dg01MdIndi::Direct), n, &mut xr, &mut xi);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

fn bench_dg01nd(c: &mut Criterion) {
    use slicot_rs::dg01::dg01nd::{dg01nd, Dg01NdIndi};
    let mut group = c.benchmark_group("dg01/dg01nd");
    for &n in SIZE_LADDER_POW2 {
        let xr = vec_n(n + 1);
        let xi = vec_n(n + 1);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut xr = xr.clone();
                let mut xi = xi.clone();
                let info = dg01nd(black_box(Dg01NdIndi::Direct), n, &mut xr, &mut xi);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

// --- Stub routines: (n, m) -> i32. One group per module. ---
macro_rules! bench_stub_nm {
    ($c:expr, $mod_path:expr, $routines:expr) => {
        for (name, f) in $routines {
            let mut group = $c.benchmark_group(format!("{}/{}", $mod_path, name));
            for &n in common::SIZE_LADDER_N {
                let m = common::m_from_n(n);
                group.bench_function(format!("n{n}_m{m}"), |bencher| {
                    bencher.iter(|| f(black_box(n), black_box(m)));
                });
            }
            group.finish();
        }
    };
}

fn bench_stubs_ab09(c: &mut Criterion) {
    use slicot_rs::ab09::{
        ab09ad::ab09ad, ab09ax::ab09ax, ab09bd::ab09bd, ab09bx::ab09bx, ab09cd::ab09cd,
        ab09cx::ab09cx, ab09dd::ab09dd, ab09ed::ab09ed, ab09fd::ab09fd, ab09gd::ab09gd,
        ab09hd::ab09hd, ab09hx::ab09hx, ab09hy::ab09hy, ab09id::ab09id, ab09ix::ab09ix,
        ab09iy::ab09iy, ab09jd::ab09jd, ab09jv::ab09jv, ab09jw::ab09jw, ab09jx::ab09jx,
        ab09kd::ab09kd, ab09kx::ab09kx, ab09md::ab09md, ab09nd::ab09nd,
    };
    let routines: &[(&str, fn(usize, usize) -> i32)] = &[
        ("ab09ad", ab09ad),
        ("ab09ax", ab09ax),
        ("ab09bd", ab09bd),
        ("ab09bx", ab09bx),
        ("ab09cd", ab09cd),
        ("ab09cx", ab09cx),
        ("ab09dd", ab09dd),
        ("ab09ed", ab09ed),
        ("ab09fd", ab09fd),
        ("ab09gd", ab09gd),
        ("ab09hd", ab09hd),
        ("ab09hx", ab09hx),
        ("ab09hy", ab09hy),
        ("ab09id", ab09id),
        ("ab09ix", ab09ix),
        ("ab09iy", ab09iy),
        ("ab09jd", ab09jd),
        ("ab09jv", ab09jv),
        ("ab09jw", ab09jw),
        ("ab09jx", ab09jx),
        ("ab09kd", ab09kd),
        ("ab09kx", ab09kx),
        ("ab09md", ab09md),
        ("ab09nd", ab09nd),
    ];
    bench_stub_nm!(c, "ab09", routines);
}

fn bench_stubs_ab13(c: &mut Criterion) {
    use slicot_rs::ab13::{
        ab13dd::ab13dd, ab13ed::ab13ed, ab13fd::ab13fd, ab13hd::ab13hd, ab13md::ab13md,
    };
    let routines: &[(&str, fn(usize, usize) -> i32)] = &[
        ("ab13dd", ab13dd),
        ("ab13ed", ab13ed),
        ("ab13fd", ab13fd),
        ("ab13hd", ab13hd),
        ("ab13md", ab13md),
    ];
    bench_stub_nm!(c, "ab13", routines);
}

fn bench_stubs_ab07(c: &mut Criterion) {
    use slicot_rs::ab07::{ab07md::ab07md, ab07nd::ab07nd};
    let routines: &[(&str, fn(usize, usize) -> i32)] = &[
        ("ab07md", ab07md),
        ("ab07nd", ab07nd),
    ];
    bench_stub_nm!(c, "ab07", routines);
}

fn bench_stubs_ab08(c: &mut Criterion) {
    use slicot_rs::ab08::{
        ab08nd::ab08nd, ab08nw::ab08nw, ab08ny::ab08ny, ab08nz::ab08nz,
    };
    let routines: &[(&str, fn(usize, usize) -> i32)] = &[
        ("ab08nd", ab08nd),
        ("ab08nw", ab08nw),
        ("ab08ny", ab08ny),
        ("ab08nz", ab08nz),
    ];
    bench_stub_nm!(c, "ab08", routines);
}

fn bench_stubs_ab08_nmp(c: &mut Criterion) {
    use slicot_rs::ab08::{ab08mz::ab08mz, ab08nx::ab08nx};
    let mut group = c.benchmark_group("ab08/ab08mz");
    for &n in SIZE_LADDER_N {
        let m = m_from_n(n);
        let p = p_from_n(n);
        group.bench_function(format!("n{n}_m{m}_p{p}"), |bencher| {
            bencher.iter(|| ab08mz(black_box(n), black_box(m), black_box(p)));
        });
    }
    group.finish();
    let mut group = c.benchmark_group("ab08/ab08nx");
    for &n in SIZE_LADDER_N {
        let m = m_from_n(n);
        let p = p_from_n(n);
        group.bench_function(format!("n{n}_m{m}_p{p}"), |bencher| {
            bencher.iter(|| ab08nx(black_box(n), black_box(m), black_box(p)));
        });
    }
    group.finish();
}

fn bench_stubs_ag07_ag08_ab8n(c: &mut Criterion) {
    use slicot_rs::ag07::ag07bd::ag07bd;
    use slicot_rs::ag08::{ag08bd::ag08bd, ag08by::ag08by, ag08bz::ag08bz};
    use slicot_rs::ab8n::ab8nxz::ab8nxz;
    use slicot_rs::ag8b::ag8byz::ag8byz;
    let routines_ag07: &[(&str, fn(usize, usize) -> i32)] = &[("ag07bd", ag07bd)];
    bench_stub_nm!(c, "ag07", routines_ag07);
    let routines_ag08: &[(&str, fn(usize, usize) -> i32)] = &[
        ("ag08bd", ag08bd),
        ("ag08by", ag08by),
        ("ag08bz", ag08bz),
    ];
    bench_stub_nm!(c, "ag08", routines_ag08);
    let routines_ab8n: &[(&str, fn(usize, usize) -> i32)] = &[("ab8nxz", ab8nxz)];
    bench_stub_nm!(c, "ab8n", routines_ab8n);
    let routines_ag8b: &[(&str, fn(usize, usize) -> i32)] = &[("ag8byz", ag8byz)];
    bench_stub_nm!(c, "ag8b", routines_ag8b);
}

fn bench_stubs_fb01_fd01(c: &mut Criterion) {
    use slicot_rs::fb01::{
        fb01qd::fb01qd, fb01rd::fb01rd, fb01sd::fb01sd, fb01td::fb01td, fb01vd::fb01vd,
    };
    use slicot_rs::fd01::fd01ad::fd01ad;
    let routines: &[(&str, fn(usize, usize) -> i32)] = &[
        ("fb01qd", fb01qd),
        ("fb01rd", fb01rd),
        ("fb01sd", fb01sd),
        ("fb01td", fb01td),
        ("fb01vd", fb01vd),
    ];
    bench_stub_nm!(c, "fb01", routines);
    let mut group = c.benchmark_group("fd01/fd01ad");
    for &n in SIZE_LADDER_N {
        let m = m_from_n(n);
        group.bench_function(format!("n{n}_m{m}"), |bencher| {
            bencher.iter(|| fd01ad(black_box(n), black_box(m)));
        });
    }
    group.finish();
}

fn bench_stubs_ib01_ib03(c: &mut Criterion) {
    use slicot_rs::ib01::{
        ib01ad::ib01ad, ib01bd::ib01bd, ib01cd::ib01cd, ib01md::ib01md, ib01my::ib01my,
        ib01nd::ib01nd, ib01od::ib01od, ib01oy::ib01oy, ib01pd::ib01pd, ib01px::ib01px,
        ib01py::ib01py, ib01qd::ib01qd, ib01rd::ib01rd,
    };
    use slicot_rs::ib03::{ib03ad::ib03ad, ib03bd::ib03bd};
    let routines_ib01: &[(&str, fn(usize, usize) -> i32)] = &[
        ("ib01ad", ib01ad),
        ("ib01bd", ib01bd),
        ("ib01cd", ib01cd),
        ("ib01md", ib01md),
        ("ib01my", ib01my),
        ("ib01nd", ib01nd),
        ("ib01od", ib01od),
        ("ib01oy", ib01oy),
        ("ib01pd", ib01pd),
        ("ib01px", ib01px),
        ("ib01py", ib01py),
        ("ib01qd", ib01qd),
        ("ib01rd", ib01rd),
    ];
    bench_stub_nm!(c, "ib01", routines_ib01);
    let routines_ib03: &[(&str, fn(usize, usize) -> i32)] = &[
        ("ib03ad", ib03ad),
        ("ib03bd", ib03bd),
    ];
    bench_stub_nm!(c, "ib03", routines_ib03);
}

fn bench_stubs_dgeg(c: &mut Criterion) {
    use slicot_rs::dgeg::{dgegs::dgegs, dgegv::dgegv};
    let routines: &[(&str, fn(usize, usize) -> i32)] = &[("dgegs", dgegs), ("dgegv", dgegv)];
    bench_stub_nm!(c, "dgeg", routines);
}

// --- DE01PD, DF01MD, DG01OD (stubs with n and buffers) ---
fn bench_de01pd(c: &mut Criterion) {
    use slicot_rs::de01::de01pd::de01pd;
    let mut group = c.benchmark_group("de01/de01pd");
    for &n in SIZE_LADDER_POW2 {
        let a = vec_n(n);
        let b = vec_n(n);
        let w = vec_n(n);
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut b = b.clone();
                let mut w = w.clone();
                let _ = de01pd(black_box(true), black_box(false), n, &mut a, &mut b, &mut w);
            });
        });
    }
    group.finish();
}

fn bench_df01md(c: &mut Criterion) {
    use slicot_rs::df01::df01md::df01md;
    let mut group = c.benchmark_group("df01/df01md");
    for &n in SIZE_LADDER_N {
        let a = vec_n(n);
        let dwork = vec![0.0; n * 2];
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut dwork = dwork.clone();
                let _ = df01md(black_box(b'S'), n, 1.0, &mut a, &mut dwork);
            });
        });
    }
    group.finish();
}

fn bench_dg01od(c: &mut Criterion) {
    use slicot_rs::dg01::dg01od::dg01od;
    let mut group = c.benchmark_group("dg01/dg01od");
    for &n in SIZE_LADDER_N {
        let a = vec_n(n);
        let w = vec_n(n);
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let mut w = w.clone();
                let _ = dg01od(black_box(0), black_box(0), n, &mut a, &mut w);
            });
        });
    }
    group.finish();
}

fn bench_dk01md(c: &mut Criterion) {
    use slicot_rs::dk01::dk01md::{dk01md, Dk01MdType};
    let mut group = c.benchmark_group("dk01/dk01md");
    for &n in SIZE_LADDER_N {
        if n <= 1 {
            continue;
        }
        let a = vec_n(n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("n{n}"), |bencher| {
            bencher.iter(|| {
                let mut a = a.clone();
                let info = dk01md(black_box(Dk01MdType::Hann), n, &mut a);
                assert_eq!(info, 0);
            });
        });
    }
    group.finish();
}

// --- BB01-BB04, BD01, BD02 (benchmark drivers: fixed params) ---
fn bench_bb01ad(c: &mut Criterion) {
    use slicot_rs::bb01::bb01ad::bb01ad;
    let mut group = c.benchmark_group("bb01/bb01ad");
    group.bench_function("fixed", |bencher| {
        bencher.iter(|| {
            let mut n = 0usize;
            let mut m = 0usize;
            let mut p = 0usize;
            let mut vec = [false; 10];
            let mut a = [0.0; 4];
            let mut b = [0.0; 2];
            let mut c = [0.0; 4];
            let mut g = [0.0; 4];
            let mut q = [0.0; 4];
            let mut x = [0.0; 4];
            let _ = bb01ad(
                [1, 1],
                &mut n,
                &mut m,
                &mut p,
                &mut vec,
                &mut a,
                2,
                &mut b,
                2,
                &mut c,
                2,
                &mut g,
                2,
                &mut q,
                2,
                &mut x,
                2,
            );
        });
    });
    group.finish();
}

fn bench_bb02ad(c: &mut Criterion) {
    use slicot_rs::bb02::bb02ad::bb02ad;
    let mut group = c.benchmark_group("bb02/bb02ad");
    group.bench_function("fixed", |bencher| {
        bencher.iter(|| {
            let mut n = 0usize;
            let mut m = 0usize;
            let mut p = 0usize;
            let mut vec = [false; 10];
            let mut a = [0.0; 4];
            let mut b = [0.0; 2];
            let mut c = [0.0; 4];
            let mut g = [0.0; 4];
            let mut q = [0.0; 4];
            let mut x = [0.0; 4];
            let _ = bb02ad(
                [1, 1],
                &mut n,
                &mut m,
                &mut p,
                &mut vec,
                &mut a,
                2,
                &mut b,
                2,
                &mut c,
                2,
                &mut g,
                2,
                &mut q,
                2,
                &mut x,
                2,
            );
        });
    });
    group.finish();
}

fn bench_bb03ad(c: &mut Criterion) {
    use slicot_rs::bb03::bb03ad::bb03ad;
    let mut group = c.benchmark_group("bb03/bb03ad");
    group.bench_function("fixed", |bencher| {
        bencher.iter(|| {
            let mut n = 0usize;
            let mut m = 0usize;
            let mut p = 0usize;
            let mut vec = [false; 10];
            let mut a = [0.0; 4];
            let mut b = [0.0; 2];
            let mut c = [0.0; 4];
            let mut x = [0.0; 4];
            let _ = bb03ad(
                [1, 1],
                &mut n,
                &mut m,
                &mut p,
                &mut vec,
                &mut a,
                2,
                &mut b,
                2,
                &mut c,
                2,
                &mut x,
                2,
            );
        });
    });
    group.finish();
}

fn bench_bb04ad(c: &mut Criterion) {
    use slicot_rs::bb04::bb04ad::bb04ad;
    let mut group = c.benchmark_group("bb04/bb04ad");
    group.bench_function("fixed", |bencher| {
        bencher.iter(|| {
            let mut n = 0usize;
            let mut m = 0usize;
            let mut p = 0usize;
            let mut vec = [false; 10];
            let mut a = [0.0; 4];
            let mut b = [0.0; 2];
            let mut c = [0.0; 4];
            let mut x = [0.0; 4];
            let _ = bb04ad(
                [1, 1],
                &mut n,
                &mut m,
                &mut p,
                &mut vec,
                &mut a,
                2,
                &mut b,
                2,
                &mut c,
                2,
                &mut x,
                2,
            );
        });
    });
    group.finish();
}

fn bench_bd01ad(c: &mut Criterion) {
    use slicot_rs::bd01::bd01ad::bd01ad;
    let mut group = c.benchmark_group("bd01/bd01ad");
    group.bench_function("fixed", |bencher| {
        bencher.iter(|| {
            let mut n = 0usize;
            let mut m = 0usize;
            let mut p = 0usize;
            let mut vec = [false; 10];
            let mut e = [0.0; 1];
            let mut a = [0.0; 1];
            let mut b = [0.0; 1];
            let mut c = [0.0; 1];
            let mut d = [0.0; 1];
            let _ = bd01ad(
                [1, 1],
                &mut n,
                &mut m,
                &mut p,
                &mut vec,
                &mut e,
                1,
                &mut a,
                1,
                &mut b,
                1,
                &mut c,
                1,
                &mut d,
                1,
            );
        });
    });
    group.finish();
}

fn bench_bd02ad(c: &mut Criterion) {
    use slicot_rs::bd02::bd02ad::bd02ad;
    let mut group = c.benchmark_group("bd02/bd02ad");
    group.bench_function("fixed", |bencher| {
        bencher.iter(|| {
            let mut n = 0usize;
            let mut m = 0usize;
            let mut p = 0usize;
            let mut vec = [false; 10];
            let mut e = [0.0; 1];
            let mut a = [0.0; 1];
            let mut b = [0.0; 1];
            let mut c = [0.0; 1];
            let mut d = [0.0; 1];
            let _ = bd02ad(
                [1, 1],
                &mut n,
                &mut m,
                &mut p,
                &mut vec,
                &mut e,
                1,
                &mut a,
                1,
                &mut b,
                1,
                &mut c,
                1,
                &mut d,
                1,
            );
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_ma02ed,
    bench_ma02es,
    bench_ma02ad,
    bench_ma02az,
    bench_ma02bd,
    bench_ma02bz,
    bench_ma02cd,
    bench_ma02dd,
    bench_ma02gd,
    bench_ma02pd,
    bench_ma02fd,
    bench_ma02rd,
    bench_mb01md,
    bench_mb01xd,
    bench_mb01xy,
    bench_mb01ss,
    bench_mb01sd,
    bench_tb01md,
    bench_dlacpy_slc,
    bench_de01od,
    bench_dlatzm,
    bench_ab01nd,
    bench_ab01od,
    bench_ab08md,
    bench_ma01ad,
    bench_dg01md,
    bench_dg01nd,
    bench_stubs_ab09,
    bench_stubs_ab13,
    bench_stubs_ab07,
    bench_stubs_ab08,
    bench_stubs_ab08_nmp,
    bench_stubs_ag07_ag08_ab8n,
    bench_stubs_fb01_fd01,
    bench_stubs_ib01_ib03,
    bench_stubs_dgeg,
    bench_de01pd,
    bench_df01md,
    bench_dg01od,
    bench_dk01md,
    bench_bb01ad,
    bench_bb02ad,
    bench_bb03ad,
    bench_bb04ad,
    bench_bd01ad,
    bench_bd02ad,
);
criterion_main!(benches);
