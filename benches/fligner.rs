use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_fligner::{Center, fligner};
use std::hint::black_box;

fn groups(per: usize, k: usize) -> Vec<Vec<f64>> {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..k)
        .map(|g| {
            let scale = 1.0 + g as f64;
            (0..per).map(|_| scale * next()).collect()
        })
        .collect()
}

fn bench_fligner(c: &mut Criterion) {
    let g = groups(1_000_000, 3);
    c.bench_function("fligner_median_3x1M", |b| {
        b.iter(|| {
            black_box(
                fligner(black_box(&g), Center::Median, 0.05)
                    .unwrap()
                    .statistic,
            )
        });
    });
    c.bench_function("fligner_mean_3x1M", |b| {
        b.iter(|| {
            black_box(
                fligner(black_box(&g), Center::Mean, 0.05)
                    .unwrap()
                    .statistic,
            )
        });
    });
}

criterion_group!(benches, bench_fligner);
criterion_main!(benches);
