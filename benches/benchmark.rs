#[macro_use]
extern crate criterion;
extern crate hungarian;

use criterion::Criterion;
use hungarian::minimize;

fn bench_hungarian(c: &mut Criterion) {
    c.bench_function_over_inputs("hungarian_NxN", |b, &&max| {
        let mut matrix = vec![0; max * max];
        let mut n: i32 = 0;
        
        for i in 0..max {
            for j in 0..max {
                matrix[max*i + j] = n;
                n += 1; 
            }
        }
        b.iter(move || minimize(&matrix, max, max))
    }, &[5, 10, 25, 50, 100, 250, 500]);
}

fn bench_hungarian_worst_case(c: &mut Criterion) {
    c.bench_function_over_inputs("hungarian_worst_case_NxN", |b, &&max| {
        let mut matrix = vec![0; max * max];
        
        for i in 0..max {
            for j in 0..max {
                matrix[max*i + j] = ((i + 1)*(j + 1)) as i32;
            }
        }

        b.iter(move || minimize(&matrix, max, max))
    }, &[5, 10, 25, 50, 100, 250, 500]);
}

criterion_group!(benches, bench_hungarian, bench_hungarian_worst_case);
criterion_main!(benches);
