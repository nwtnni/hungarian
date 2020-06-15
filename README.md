# hungarian

[![Build Status](https://travis-ci.org/nwtnni/hungarian.svg?branch=master)](https://travis-ci.org/nwtnni/hungarian)
[![License](https://img.shields.io/github/license/nwtnni/hungarian.svg)](https://raw.githubusercontent.com/nwtnni/hungarian/master/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/hungarian.svg)](https://crates.io/crates/hungarian)
[![Rustdoc](https://docs.rs/hungarian/badge.svg)](https://docs.rs/hungarian/)
![Crates.io](https://img.shields.io/crates/d/hungarian.svg)

**IMPORTANT**: The [`pathfinding`](https://github.com/samueltardieu/pathfinding) crate has a significantly
faster [implementation](https://docs.rs/pathfinding/2.0/pathfinding/kuhn_munkres/index.html) of this
algorithm (benchmarks below), uses traits to abstract over cost matrices, and is also better maintained.
I recommend using it instead.

A simple Rust implementation of the Hungarian (or Kuhn–Munkres) algorithm.
Should run in `O(n^3)` time and take `O(m*n)` space, given an `m * n` rectangular
matrix (represented as a 1D slice).

Derived and modified from [this great explanation](http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html).

## Usage

Add the following to your `Cargo.toml` file:

```
[dependencies]
hungarian = "1.1.1"
```

Add the following to the top of your binary or library:

```
extern crate hungarian;

use hungarian::minimize;
```

And you should be good to go!
[For more information, check out the documentation.](https://docs.rs/hungarian/)

## Recent Changes

- 1.1.1
  - Version bump so the [`pathfinding`](https://github.com/samueltardieu/pathfinding) redirect appears on `crates.io`.
- 1.1.0
  - Greatly optimized performance (by a factor of 2-4 on benchmarks on matrices from 5x5 to 100x100)
  - Now uses `num-trait` to take generic integer weights
  - Now backed by `ndarray` to scale better with larger inputs
- 1.0.0 
  - Greatly improved source code documentation
  - Renamed `hungarian` function to `minimize`
  - Now handle arbitrary rectangular matrices
  - Added more test cases to cover non-square matrices
  - Now returns `Vec<Option<Usize>>` to handle when not all columns are assigned to rows
- 0.1.0
  - Initial release 
  - Working base algorithm, but only works for square matrices.
  - Not well documented

## Notes

Instead of using splitting logic across files and helper functions, I tried to simplify and
condense the above explanation into a single, simple function while maintaining correctness.
After trawling the web for test cases, I'm reasonably confident that my implementation works,
even though the end result looks fairly different.

Please let me know if you find any bugs!

## Performance

Benchmarks were obtained using Criterion.rs, with the following two
types of cost matrices:

```
     Worst Case           |       Generic Case
                          |
   -------------          |       -------------
   | 1 | 2 | 3 | ...      |       | 1 | 2 | 3 |
   -------------          |       -------------
   | 2 | 4 | 6 | ...      |       | 4 | 5 | 6 |
   -------------          |       -------------
   | 3 | 6 | 9 | ...      |       | 7 | 8 | 9 |
   -------------          |       -------------
     .   .   .            |
     .   .   .            |
     .   .   .            |
                          |
C(i, j) = (i + 1)(j + 1)  |  C(i, j) = (i * width) + j
```

#### Criterion Results

| Cost Matrix | Matrix Size | `hungarian` Average Runtime | `pathfinding` Average Runtime |
|:------------|------------:|----------------------------:|------------------------------:|
| Worst-Case  |     5 x   5 |                     2.42 us |                       1.19 us |
| Worst-Case  |    10 x  10 |                    20.38 us |                       4.24 us |
| Worst-Case  |    25 x  25 |                   546.88 us |                      59.66 us |
| Worst-Case  |    50 x  50 |                     6.97 ms |                     422.05 us |
| Generic     |     5 x   5 |                     1.75 us |                     871.24 ns |
| Generic     |    10 x  10 |                     7.49 us |                       3.50 us |
| Generic     |    25 x  25 |                    86.33 us |                      33.91 us |
| Generic     |    50 x  50 |                   556.48 us |                     285.69 us |
| Generic     |   100 x 100 |                     3.97 ms |                       1.93 ms |

Measured on a quad-core 2.6GHz Intel(R) i7-6700HQ with
16GB RAM; using Ubuntu 16.04 Linux x86\_64 4.8.0-53-generic
