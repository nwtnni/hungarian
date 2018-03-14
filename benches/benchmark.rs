#[macro_use]
extern crate criterion;
extern crate hungarian;

use criterion::Criterion;
use hungarian::hungarian;

fn bench_hungarian(c: &mut Criterion) {
    c.bench_function_over_inputs("hungarian_NxN", |b, &&max| {
        let mut matrix = vec![0; max * max];
        let mut n: u64 = 0;
        
        for i in 0..max {
            for j in 0..max {
                matrix[max*i + j] = n;
                n += 1; 
            }
        }
        b.iter(move || hungarian(&matrix, max, max))
    }, &[5, 10, 25, 50]);
}

criterion_group!(benches, bench_hungarian);
criterion_main!(benches);

