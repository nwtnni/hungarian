#[macro_use]
extern crate criterion;
extern crate hungarian;
extern crate pathfinding;

use criterion::Criterion;
use hungarian::minimize;
use pathfinding::kuhn_munkres::kuhn_munkres_min;
use pathfinding::matrix::Matrix;

fn bench_hungarian(c: &mut Criterion) {
    c.bench_function_over_inputs(
        "hungarian_NxN",
        |b, &&max| {
            let mut matrix = vec![0; max * max];
            let mut n = 0;
            for i in 0..max {
                for j in 0..max {
                    matrix[max * i + j] = n;
                    n += 1;
                }
            }
            b.iter(move || minimize(&matrix, max, max))
        },
        &[5, 10, 25, 50, 100],
    );
}

fn bench_hungarian_worst_case(c: &mut Criterion) {
    c.bench_function_over_inputs(
        "hungarian_worst_case_NxN",
        |b, &&max| {
            let mut matrix = vec![0; max * max];
            for i in 0..max {
                for j in 0..max {
                    matrix[max * i + j] = ((i + 1) * (j + 1)) as i32;
                }
            }
            b.iter(move || minimize(&matrix, max, max))
        },
        &[5, 10, 25, 50],
    );
}

fn bench_pathfinding_hungarian(c: &mut Criterion) {
    c.bench_function_over_inputs(
        "pathfinding_hungarian_NxN",
        |b, &&max| {
            let mut matrix = Matrix::new(max, max, 0);
            let mut n: i32 = 0;
            for i in 0..max {
                for j in 0..max {
                    matrix[&(i, j)] = n;
                    n += 1;
                }
            }
            b.iter(move || kuhn_munkres_min(&matrix))
        },
        &[5, 10, 25, 50, 100],
    );
}

fn bench_pathfinding_hungarian_worst_case(c: &mut Criterion) {
    c.bench_function_over_inputs(
        "pathfinding_hungarian_worst_case_NxN",
        |b, &&max| {
            let mut matrix = Matrix::new(max, max, 0);
            for i in 0..max {
                for j in 0..max {
                    matrix[&(i, j)] = ((i + 1) * (j + 1)) as i32;
                }
            }
            b.iter(move || kuhn_munkres_min(&matrix))
        },
        &[5, 10, 25, 50],
    );
}

criterion_group!(
    benches,
    bench_hungarian,
    bench_hungarian_worst_case,
    bench_pathfinding_hungarian,
    bench_pathfinding_hungarian_worst_case,
);

criterion_main!(benches);
