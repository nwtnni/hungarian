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

TODO
