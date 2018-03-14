# hungarian

[![Build Status](https://travis-ci.org/nwtnni/hungarian.svg?branch=master)](https://travis-ci.org/nwtnni/hungarian)

A simple Rust implementation of the Hungarian (or Kuhn–Munkres) algorithm.

Derived and modified from [this great explanation](http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html).

## Usage

There's only one function:

```rust
hungarian(matrix: &[u64], w: usize, h: usize) -> Vec<usize>`
```

#### Parameters

- `matrix` is a 2D matrix flattened into a 1D slice. Must be square; rectangular matrices can be padded with zeros.
- `w` is the width of the matrix (i.e. number of columns).
- `h` is the hight of the matrix (i.e. number of rows).

#### Returns

Returns a `Vec` where `Vec[i]` is the assignment for row `i`.

#### Example

```rust
use hungarian::hungarian;

fn main() {
  
  let matrix = vec![
    0, 1, 1, 1,
    1, 1, 0, 1,
    1, 0, 1, 1,
    1, 1, 1, 0,
  ];

  assert_eq!(hungarian(&matrix, 4, 4), vec![0, 2, 1, 3]);
}

```

## Notes

There's only one dependency ([fixedbitset](https://github.com/bluss/fixedbitset)) and one source file in this crate.

Instead of using splitting logic across files and helper functions, I tried to simplify and
condense the above explanation into a single, simple function while maintaining correctness.
After trawling the web for test cases, I'm reasonably confident that my implementation works,
even though the end result looks fairly different.

Please let me know if you find any bugs!

## Performance

Benchmarks were obtained using Criterion.rs, with the following two
types of cost matrices:

```
Worst Case               |  Generic Case
                         |  
-------------            |  -------------
| 1 | 2 | 3 | ...        |  | 1 | 2 | 3 |
-------------            |  -------------
| 2 | 4 | 6 | ...        |  | 4 | 5 | 6 |
-------------            |  -------------
| 3 | 6 | 9 | ...        |  | 7 | 8 | 9 |
-------------            |  -------------
  .   .   .              |  
  .   .   .              |   
  .   .   .              |     
                         |  
C(i, j) = (i + 1)(j + 1) |  C(i, j) = (i * width) + j
```

#### Criterion Results

| Cost Matrix | Matrix Size | Average Runtime |
-----------------------------------------------
| Generic     |  5 x  5     |   2.54 us       |
| Generic     | 10 x 10     |  16.31 us       |
| Generic     | 25 x 25     | 355.27 us       |
| Generic     | 50 x 50     |   4.12 ms       |
| Worst-Case  |  5 x  5     |   3.21 us       |
| Worst-Case  | 10 x 10     |  35.77 us       |
| Worst-Case  | 25 x 25     |   1.14 ms       |
| Worst-Case  | 50 x 50     |  18.59 ms       |

```
Benchmarking hungarian_NxN/5
Benchmarking hungarian_NxN/5: Warming up for 3.0000 s
Benchmarking hungarian_NxN/5: Collecting 100 samples in estimated 5.0094 s (2.0M iterations)
Benchmarking hungarian_NxN/5: Analyzing
hungarian_NxN/5         time:   [2.5008 us 2.5366 us 2.5796 us]

Benchmarking hungarian_NxN/10
Benchmarking hungarian_NxN/10: Warming up for 3.0000 s
Benchmarking hungarian_NxN/10: Collecting 100 samples in estimated 5.0369 s (308k iterations)
Benchmarking hungarian_NxN/10: Analyzing
hungarian_NxN/10        time:   [16.196 us 16.311 us 16.447 us]

Benchmarking hungarian_NxN/25
Benchmarking hungarian_NxN/25: Warming up for 3.0000 s
Benchmarking hungarian_NxN/25: Collecting 100 samples in estimated 5.4580 s (15k iterations)
Benchmarking hungarian_NxN/25: Analyzing
hungarian_NxN/25        time:   [355.27 us 357.79 us 360.33 us]

Benchmarking hungarian_NxN/50
Benchmarking hungarian_NxN/50: Warming up for 3.0000 s
Benchmarking hungarian_NxN/50: Collecting 100 samples in estimated 21.105 s (5050 iterations)
Benchmarking hungarian_NxN/50: Analyzing
hungarian_NxN/50        time:   [4.1241 ms 4.1487 ms 4.1749 ms]

Benchmarking hungarian_worst_case_NxN/5
Benchmarking hungarian_worst_case_NxN/5: Warming up for 3.0000 s
Benchmarking hungarian_worst_case_NxN/5: Collecting 100 samples in estimated 5.0095 s (1.6M iterations)
Benchmarking hungarian_worst_case_NxN/5: Analyzing
hungarian_worst_case_NxN/5
                        time:   [3.2128 us 3.2256 us 3.2401 us]

Benchmarking hungarian_worst_case_NxN/10
Benchmarking hungarian_worst_case_NxN/10: Warming up for 3.0000 s
Benchmarking hungarian_worst_case_NxN/10: Collecting 100 samples in estimated 5.1088 s (141k iterations)
Benchmarking hungarian_worst_case_NxN/10: Analyzing
hungarian_worst_case_NxN/10
                        time:   [35.767 us 35.975 us 36.226 us]

Benchmarking hungarian_worst_case_NxN/25
Benchmarking hungarian_worst_case_NxN/25: Warming up for 3.0000 s
Benchmarking hungarian_worst_case_NxN/25: Collecting 100 samples in estimated 5.7591 s (5050 iterations)
Benchmarking hungarian_worst_case_NxN/25: Analyzing
hungarian_worst_case_NxN/25
                        time:   [1.1416 ms 1.1525 ms 1.1652 ms]

Benchmarking hungarian_worst_case_NxN/50
Benchmarking hungarian_worst_case_NxN/50: Warming up for 3.0000 s
Benchmarking hungarian_worst_case_NxN/50: Collecting 100 samples in estimated 91.351 s (5050 iterations)
Benchmarking hungarian_worst_case_NxN/50: Analyzing
hungarian_worst_case_NxN/50
                        time:   [18.592 ms 18.762 ms 18.942 ms]
```
