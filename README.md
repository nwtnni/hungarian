# hungarian

[![Build Status](https://travis-ci.org/nwtnni/hungarian.svg?branch=master)](https://travis-ci.org/nwtnni/hungarian)
[![License](https://img.shields.io/github/license/nwtnni/hungarian.svg)](https://raw.githubusercontent.com/nwtnni/hungarian/master/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/hungarian.svg)](https://crates.io/crates/hungarian)
[![Rustdoc](https://docs.rs/hungarian/badge.svg)](https://docs.rs/hungarian/)
![Crates.io](https://img.shields.io/crates/d/hungarian.svg)

A simple Rust implementation of the Hungarian (or Kuhn–Munkres) algorithm.
Should run in `O(n^3)` time and take `O(m*n)` space, given an `m * n` rectangular
matrix (represented as a 1D slice).

Derived and modified from [this great explanation](http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html).

## Usage

Check out the [documentation!](https://docs.rs/hungarian/)

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

| Cost Matrix | Matrix Size | Average Runtime |  Iterations |
|:------------|------------:|----------------:|------------:|
| Worst-Case  |     5 x   5 |         3.13 us | 1\_600\_000 |
| Worst-Case  |    10 x  10 |        36.85 us |    141\_000 |
| Worst-Case  |    25 x  25 |         1.20 ms |      5\_050 |
| Worst-Case  |    50 x  50 |        19.62 ms |      5\_050 |
| Generic     |     5 x   5 |         2.22 us | 2\_100\_000 |
| Generic     |    10 x  10 |        12.69 us |    379\_000 |
| Generic     |    25 x  25 |       182.48 us |     30\_000 |
| Generic     |    50 x  50 |         1.84 ms |      5\_050 |
| Generic     |   100 x 100 |        19.17 ms |      5\_050 |

Measured on a quad-core 2.6GHz Intel(R) i7-6700HQ with
16GB RAM; using Ubuntu 16.04 Linux x86\_64 4.8.0-53-generic
