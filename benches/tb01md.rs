//! Criterion benchmark for TB01MD (controller Hessenberg reduction).
//! Compare performance across matrix sizes.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::DMatrix;
use slicot_rs::tb01::tb01md::{tb01md, JobU, Uplo};

fn bench_tb01md(c: &mut Criterion) {
    let mut group = c.benchmark_group("tb01md");
    for n in [4usize, 8, 16, 32] {
        let m = n / 2;
        let a = DMatrix::from_fn(n, n, |i, j| (i + j) as f64 * 0.1);
        let b = DMatrix::from_fn(n, m, |i, j| (i * 2 + j) as f64 * 0.1);
        group.bench_function(format!("n{}_m{}", n, m), |bencher| {
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

criterion_group!(benches, bench_tb01md);
criterion_main!(benches);
