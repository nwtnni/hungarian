extern crate fixedbitset;

use fixedbitset::FixedBitSet;

/// Internal helper macro for converting from a 2D index to a 1D index.
macro_rules! index {
    ($w:expr, $i:expr, $j:expr) => (($w*$i) + $j)
}

/// # Summary
///
/// Implementation of the Hungarian / Munkres Assignment Algorithm.
///
/// Given a cost matrix, this algorithm finds a perfect matching such
/// that the total cost is minimized. [Follows the general outline explained here.][1]
///
/// # Requires
///
/// - `matrix` is square (i.e. `w` == `h`). `matrix` can be padded with 0's if
///     non-square, without changing the solution.
/// - Should also work if `w` > `h`, but there are currently no test cases for this scenario.
///
/// # Takes
///
/// - `matrix: &[u64]`: a 1D slice in row-major order, representing a 2D matrix
/// - `w`: width of `matrix` (i.e. number of columns)
/// - `h`: height of `matrix` (i.e. number of rows)
///
/// # Returns
///
/// - `v`: A Vec such that `v[i]` is the column assigned to row `i`.
///
/// # Examples
///
/// ```rust
/// extern crate hungarian;
///
/// use hungarian::minimize;
///
/// fn main() {
///
///     let width = 3;
///     let height = 3;
///
///     let matrix = vec![
///         1, 2, 1,
///         4, 5, 6,
///         7, 8, 9,
///     ];
///
///     let assignment = minimize(&matrix, width, height);
///
///     assert_eq!(&assignment, &vec![2, 1, 0]);
///
///     let cost: u64 = assignment.iter()
///         .enumerate()
///         .map(|(i, &j)| matrix[i*width + j])
///         .sum();
///
///     assert_eq!(cost, 13);
/// }
///
/// [1]: http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
/// ```
pub fn minimize(matrix: &[u64], w: usize, h: usize) -> Vec<usize> {

    //********************************************//
    //                                            //
    //                   Step 0                   //
    //                                            //
    //********************************************//

    // Number of assignments to find
    let target = usize::min(w, h);

    // Clone of the provided slice
    let mut matrix = Vec::from(matrix);

    // The set of starred zero entries (flattened into 1D)
    let mut stars = FixedBitSet::with_capacity(w * h);

    // The set of primed zero entries (flattened into 1D)
    let mut primes = FixedBitSet::with_capacity(w * h);

    // The set of covered row indices
    let mut row_cover = FixedBitSet::with_capacity(h);

    // The set of covered column indices
    let mut col_cover = FixedBitSet::with_capacity(w);

    //********************************************//
    //                                            //
    //                   Step 1                   //
    //                                            //
    //********************************************//

    // Reduce each row by its smallest element
    for row in matrix.chunks_mut(w) {
        let min = row.iter().min().unwrap().clone();
        row.iter_mut().for_each(|cost| *cost -= min);
    }

    //********************************************//
    //                                            //
    //                   Step 2                   //
    //                                            //
    //********************************************//

    // Find a zero (Z):
    // - If there is no starred zero in its row or column, then star it.
    // - Use row_cover, col_cover to keep track of stars.
    for i in 0..h {
        for j in 0..w {
            let k = index!(w, i, j);
            if matrix[k] == 0 && !row_cover[i] && !col_cover[j] {
                stars.insert(k);
                row_cover.insert(i);
                col_cover.insert(j);
            }
        }
    }

    // Reset cover
    row_cover.clear();
    col_cover.clear();
    let mut verify = true;

    loop {

        //********************************************//
        //                                            //
        //                   Step 3                   //
        //                                            //
        //********************************************//

        // Cover each column with a starred zero.
        // If there are `target` starred zeros, then we're done.
        if verify {
            stars.ones().for_each(|k| col_cover.insert(k % w));
            if col_cover.count_ones(..) == target {
                return stars.ones().map(|k| k % w).collect()
            }
        }

        //********************************************//
        //                                            //
        //                   Step 4                   //
        //                                            //
        //********************************************//

        // Find an uncovered zero and prime it
        let mut uncovered = None;

        'outer : for i in 0..h {
            for j in 0..w {
                if row_cover[i] || col_cover[j] { continue }
                let k = index!(w, i, j);
                if matrix[k] == 0 {
                    uncovered = Some((i, j));
                    primes.insert(k);
                    break 'outer;
                }
            }
        }

        // No uncovered zeros left
        if let None = uncovered {

            //********************************************//
            //                                            //
            //                   Step 6                   //
            //                                            //
            //********************************************//

            // Find minimum uncovered value
            let mut min = u64::max_value();
            for i in 0..h {
                for j in 0..w {
                    if row_cover[i] || col_cover[j] { continue }
                    let value = matrix[index!(w, i, j)];
                    min = if value < min { value } else { min };
                }
            }

            // Add minimum to covered rows
            // Subtract minimum from uncovered columns
            for i in 0..h {
                for j in 0..w {
                    if row_cover[i] { matrix[index!(w, i, j)] += min }
                    if !col_cover[j] { matrix[index!(w, i, j)] -= min }
                }
            }

            // Return to Step 4
            // - Skip rest of this loop
            // - Set `verify` to false to skip Step 3
            verify = false;
            continue
        }

        // Find starred zero in the same row
        let (i, j) = uncovered.unwrap();
        let starred = (0..w).filter(|&j| {
            stars[index!(w, i, j)]
        }).next();

        // Starred zero exists:
        // - Cover row of uncovered zero from Step 4
        // - Uncover column of starred zero
        // - Repeat Step 4
        if let Some(j) = starred {
            row_cover.insert(i);
            col_cover.set(j, false);
            verify = false;
            continue
        }

        //********************************************//
        //                                            //
        //                   Step 5                   //
        //                                            //
        //********************************************//

        // Construct an alternating path of stars and primes
        let mut path = vec![(i, j)];
        loop {
            let (_, j) = path[path.len() - 1];

            // Find starred zero in same column
            let next_star = (0..h).filter(|&i| {
                stars[index!(w, i, j)]
            }).next();

            if let None = next_star { break }
            let i = next_star.unwrap();
            path.push((i, j));

            // Find primed zero in same row
            let next_prime = (0..w).filter(|&j| {
                primes[index!(w, i, j)]
            }).next();

            let j = next_prime.expect("Guaranteed to exist");
            path.push((i, j));
        }

        // Unstar each starred zero
        // Star each primed zero
        for (i, j) in path {
            let k = index!(w, i, j);
            stars.set(k, primes[k]);
        }

        // Reset cover
        row_cover.clear();
        col_cover.clear();

        // Erase primes and return to Step 3
        primes.clear();
        verify = true;
    }
}

#[cfg(test)]
mod tests {
    use minimize;

    #[test]
    fn test_basic() {
        let matrix = vec![
            1, 2, 2,
            2, 1, 2,
            2, 2, 1,
        ];
        assert_eq!(minimize(&matrix, 3, 3), vec![0, 1, 2])
    }

    #[test]
    fn test_increasing() {
        let matrix = vec![
            1,   2,  3,  1,
            5,   6,  7,  8,
            9,  10, 11, 12,
            13, 14, 15, 16,
        ];
        assert_eq!(minimize(&matrix, 4, 4), vec![3, 2, 1, 0])
    }

    // From http://www.math.harvard.edu/archive/20_spring_05/handouts/assignment_overheads.pdf
    #[test]
    fn test_sales_example() {
        let matrix = vec![
            250, 400, 350,
            400, 600, 350,
            200, 400, 250,
        ];
        assert_eq!(minimize(&matrix, 3, 3), vec![1, 2, 0]);
    }

    // From https://brilliant.org/wiki/hungarian-matching/
    #[test]
    fn test_party_example() {
        let matrix = vec![
            108, 125, 150,
            150, 135, 175,
            122, 148, 250,
        ];
        assert_eq!(minimize(&matrix, 3, 3), vec![2, 1, 0]);
    }

    // From https://en.wikipedia.org/wiki/Hungarian_algorithm#Matrix_interpretation
    #[test]
    fn test_wiki() {
        let matrix = vec![
            0, 1, 2, 3,
            4, 5, 6, 0,
            0, 2, 4, 5,
            3, 0, 0, 9,
        ];
        assert_eq!(minimize(&matrix, 4, 4), vec![1, 3, 0, 2]);
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_5() {
        let matrix = vec![
            12,  9, 27, 10, 23,
             7, 13, 13, 30, 19,
            25, 18, 26, 11, 26,
             9, 28, 26, 23, 13,
            16, 16, 24,  6,  9,
        ];
        assert_eq!(minimize(&matrix, 5, 5)
            .iter()
            .enumerate()
            .map(|(i, &j)| matrix[index!(5, i, j)])
            .sum::<u64>(), 51);
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_10() {
        let matrix = vec![
            37, 34, 29, 26, 19,  8,  9, 23, 19, 29,
             9, 28, 20,  8, 18, 20, 14, 33, 23, 14,
            15, 26, 12, 28,  6, 17,  9, 13, 21,  7,
             2,  8, 38, 36, 39,  5, 36,  2, 38, 27,
            30,  3, 33, 16, 21, 39,  7, 23, 28, 36,
             7,  5, 19, 22, 36, 36, 24, 19, 30,  2,
            34, 20, 13, 36, 12, 33,  9, 10, 23,  5,
             7, 37, 22, 39, 33, 39, 10,  3, 13, 26,
            21, 25, 23, 39, 31, 37, 32, 33, 38,  1,
            17, 34, 40, 10, 29, 37, 40,  3, 25,  3,
        ];
        assert_eq!(minimize(&matrix, 10, 10)
            .iter()
            .enumerate()
            .map(|(i, &j)| matrix[index!(10, i, j)])
            .sum::<u64>(), 66);
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_20() {
        let matrix = vec![
            5,  4,  3,  9,  8,  9,  3,  5,  6,  9,  4, 10,  3,  5,  6,  6,  1,  8, 10,  2,
            10, 9,  9,  2,  8,  3,  9,  9, 10,  1,  7, 10,  8,  4,  2,  1,  4,  8,  4,  8,
            10, 4,  4,  3,  1,  3,  5, 10,  6,  8,  6,  8,  4, 10,  7,  2,  4,  5,  1,  8,
            2,  1,  4,  2,  3,  9,  3,  4,  7,  3,  4,  1,  3,  2,  9,  8,  6,  5,  7,  8,
            3,  4,  4,  1,  4, 10,  1,  2,  6,  4,  5, 10,  2,  2,  3,  9, 10,  9,  9, 10,
            1, 10,  1,  8,  1,  3,  1,  7,  1,  1,  2,  1,  2,  6,  3,  3,  4,  4,  8,  6,
            1,  8,  7, 10, 10,  3,  4,  6,  1,  6,  6,  4,  9,  6,  9,  6,  4,  5,  4,  7,
            8, 10,  3,  9,  4,  9,  3,  3,  4,  6,  4,  2,  6,  7,  7,  4,  4,  3,  4,  7,
            1,  3,  8,  2,  6,  9,  2,  7,  4,  8, 10,  8, 10,  5,  1,  3, 10, 10,  2,  9,
            2,  4,  1,  9,  2,  9,  7,  8,  2,  1,  4, 10,  5,  2,  7,  6,  5,  7,  2,  6,
            4,  5,  1,  4,  2,  3,  3,  4,  1,  8,  8,  2,  6,  9,  5,  9,  6,  3,  9,  3,
            3,  1,  1,  8,  6,  8,  8,  7,  9,  3,  2,  1,  8,  2,  4,  7,  3,  1,  2,  4,
            5,  9,  8,  6, 10,  4, 10,  3,  4, 10, 10, 10,  1,  7,  8,  8,  7,  7,  8,  8,
            1,  4,  6,  1,  6,  1,  2, 10,  5, 10,  2,  6,  2,  4,  5,  5,  3,  5,  1,  5,
            5,  6,  9, 10,  6,  6, 10,  6,  4,  1,  5,  3,  9,  5,  2, 10,  9,  9,  5,  1,
            10, 9,  4,  6,  9,  5,  3,  7, 10,  1,  6,  8,  1,  1, 10,  9,  5,  7,  7,  5,
            2,  6,  6,  6,  6,  2,  9,  4,  7,  5,  3,  2, 10,  3,  4,  5, 10,  9,  1,  7,
            5,  2,  4,  9,  8,  4,  8,  2,  4,  1,  3,  7,  6,  8,  1,  6,  8,  8, 10, 10,
            9,  6,  3,  1,  8,  5,  7,  8,  7,  2,  1,  8,  2,  8,  3,  7,  4,  8,  7,  7,
            8,  4,  4,  9,  7, 10,  6,  2,  1,  5,  8,  5,  1,  1,  1,  9,  1,  3,  5,  3,
        ];
        assert_eq!(minimize(&matrix, 20, 20)
            .iter()
            .enumerate()
            .map(|(i, &j)| matrix[index!(20, i, j)])
            .sum::<u64>(), 22);
    }


    // From https://stackoverflow.com/questions/17419595/hungarian-kuhn-munkres-algorithm-oddity
    #[test]
    fn test_stack_overflow_a() {
        let matrix = vec![
            0, 7, 0, 0, 0,
            0, 8, 0, 0, 6,
            5, 0, 7, 3, 4,
            5, 0, 5, 9, 3,
            0, 4, 0, 0, 9,
        ];
        assert_eq!(minimize(&matrix, 5, 5)
            .iter()
            .enumerate()
            .map(|(i, &j)| matrix[index!(5, i, j)])
            .sum::<u64>(), 3);
        assert_eq!(minimize(&matrix, 5, 5), vec![4, 2, 3, 1, 0]);
    }

    // From https://stackoverflow.com/questions/26893961/cannot-solve-hungarian-algorithm
    #[test]
    fn test_stack_overflow_b() {
        let matrix = vec![
             0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,
            53, 207, 256, 207, 231, 348, 348, 348, 231, 244, 244,  0,  0,  0,
           240,  33,  67,  33,  56, 133, 133, 133,  56,  33,  33,  0,  0,  0,
           460, 107, 200, 107, 122, 324, 324, 324, 122,  33,  33,  0,  0,  0,
           167, 340, 396, 340, 422, 567, 567, 567, 422, 442, 442,  0,  0,  0,
           167, 367, 307, 367, 433, 336, 336, 336, 433, 158, 158,  0,  0,  0,
           160,  20,  37,  20,  31,  70,  70,  70,  31,  22,  22,  0,  0,  0,
           200, 307, 393, 307, 222, 364, 364, 364, 222, 286, 286,  0,  0,  0,
           33 , 153, 152, 153, 228, 252, 252, 252, 228,  78,  78,  0,  0,  0,
           93 , 140, 185, 140,  58, 118, 118, 118,  58,  44,  44,  0,  0,  0,
           0  ,   7,  22,   7,  19,  58,  58,  58,  19,   0,   0,  0,  0,  0,
           67 , 153, 241, 153, 128, 297, 297, 297, 128,  39,  39,  0,  0,  0,
           73 , 253, 389, 253, 253, 539, 539, 539, 253,  36,  36,  0,  0,  0,
           173, 267, 270, 267, 322, 352, 352, 352, 322, 231, 231,  0,  0,  0,
        ];
        assert_eq!(
            minimize(&matrix, 14, 14)
                .iter()
                .enumerate()
                .map(|(i, &j)| matrix[index!(14, i, j)])
                .sum::<u64>(),
            828
        );
    }

    // From https://stackoverflow.com/questions/46803600/hungarian-algorithm-wikipedia-method-doesnt-work-for-this-example
    #[test]
    fn test_stack_overflow_c() {
        let matrix = vec![
            35,  0,  0,  0,
            0 , 30,  0,  5,
            55,  5,  0, 10,
            0 , 45, 30, 45,
        ];
        assert_eq!(
            minimize(&matrix, 4, 4)
                .iter()
                .enumerate()
                .map(|(i, &j)| matrix[index!(4, i, j)])
                .sum::<u64>(),
            5
        );
        assert_eq!(minimize(&matrix, 4, 4), vec![1, 3, 2, 0]);
    }

    // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    #[test]
    fn test_stack_overflow_d() {
        let matrix = vec![
            2, 1, 0, 0, 0, 3,
            2, 0, 4, 5, 2, 7,
            0, 7, 0, 0, 0, 5,
            3, 2, 3, 1, 2, 0,
            0, 0, 6, 3, 3, 5,
            3, 4, 5, 2, 0, 3,
        ];
        assert_eq!(minimize(&matrix, 6, 6), vec![3, 1, 2, 5, 0, 4]);
    }

    // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    #[test]
    fn test_stack_overflow_e() {
        let matrix = vec![
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 2,
            0, 0, 3, 4,
        ];
        assert_eq!(minimize(&matrix, 4, 4), vec![3, 2, 1, 0]);
    }

    // From https://www.wikihow.com/Use-the-Hungarian-Algorithm
    #[test]
    fn test_wikihow_example() {
        let matrix = vec![
            10, 19, 8, 15, 19,
            10, 18, 7, 17, 19,
            13, 16, 9, 14, 19,
            12, 19, 8, 18, 19,
            14, 17, 10, 19, 19
        ];
        assert_eq!(minimize(&matrix, 5, 5), vec![0, 2, 3, 4, 1]);
    }

    #[test]
    fn test_stress() {
        for max in 1..100 {
            let mut matrix = vec![0; max * max];
            let mut n: u64 = 0;

            for i in 0..max {
                for j in 0..max {
                    matrix[index!(max, i, j)] = n;
                    n += 1;
                }
            }
            assert_eq!(minimize(&matrix, max, max), (0..max).rev().collect::<Vec<_>>());
        }
    }

    #[test]
    fn test_worst_case() {
        for max in 1..50 {
            let mut matrix = vec![0; max * max];

            for i in 0..max {
                for j in 0..max {
                    matrix[index!(max, i, j)] = ((i + 1)*(j + 1)) as u64;
                }
            }
            assert_eq!(minimize(&matrix, max, max), (0..max).rev().collect::<Vec<_>>());
        }
    }

    #[test]
    fn test_large() {
        let max = 250;
        let mut matrix = vec![0; max * max];
        let mut n: u64 = 0;

        for i in 0..max {
            for j in 0..max {
                matrix[index!(max, i, j)] = n;
                n += 1;
            }
        }
        assert_eq!(minimize(&matrix, max, max), (0..max).rev().collect::<Vec<_>>());
    }
}
