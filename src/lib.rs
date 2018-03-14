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
/// Given a rectangular cost matrix, this algorithm finds a maximal matching such
/// that the total cost is minimized. [Follows the general outline explained here.][1]
///
/// # Requires 
///
/// - `matrix` is rectangular (i.e. no ragged matrices)
///
/// # Takes
///
/// - `matrix: &[u64]`: a 1D slice in row-major order, representing a 2D matrix
/// - `height`: height of `matrix` (i.e. number of rows)
/// - `width`: width of `matrix` (i.e. number of columns)
///
/// # Returns
///
/// - `v`: A Vec where `v[i]` is:
///     - `Some(j)` if row `i` should be assigned to column `j`
///     - `None` if row `i` is not in the optimal assignment. Only possible if `width < height`.
///
/// # Panics
///
/// This function uses array indexing directly, and will almost definitely
/// panic with out-of-bounds if passed incorrect `width` or `height` arguments.
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
///     // Square matrix
///
///     let width = 3;
///     let height = 3;
///     let matrix = vec![
///         1, 2, 1,
///         4, 5, 6,
///         7, 8, 9,
///     ];
///     
///     let assignment = minimize(&matrix, height, width);
///
///     assert_eq!(&assignment, &vec![Some(2), Some(1), Some(0)]);
///
///     let cost: u64 = assignment.iter()
///         .enumerate()
///         .filter_map(|(i, &a)| {
///             a.map(|j| matrix[i*width + j])
///         })
///         .sum();
///
///     assert_eq!(cost, 13);
///
///     // Rectangular matrix (height < width)
///
///     let height = 2;
///     let width = 3;
///     let matrix = vec![
///         1, 0, 5,
///         2, 3, 1, 
///     ];
///
///     let assignment = minimize(&matrix, height, width);
///
///     assert_eq!(&assignment, &vec![Some(1), Some(2)]);
///
///     let cost: u64 = assignment.iter()
///         .enumerate()
///         .filter_map(|(i, &a)| {
///             a.map(|j| matrix[i*width + j])
///         })
///         .sum();
///
///     assert_eq!(cost, 1);
///
///     // Rectangular matrix (width < height)
///
///     let height = 3;
///     let width = 2;
///     let matrix = vec![
///         5, 5,
///         1, 0,
///         2, 3, 
///     ];
///
///     let assignment = minimize(&matrix, height, width);
///
///     assert_eq!(&assignment, &vec![None, Some(1), Some(0)]);
///
///     let cost: u64 = assignment.iter()
///         .enumerate()
///         .filter_map(|(i, &a)| {
///             a.map(|j| matrix[i*width + j])
///         })
///         .sum();
///
///     assert_eq!(cost, 2);
///
/// }
/// ```
///
/// [1]: http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
///
pub fn minimize(matrix: &[u64], height: usize, width: usize) -> Vec<Option<usize>> {

    // No possible assignment
    if height == 0 || width == 0 { return Vec::new() }

    //********************************************//
    //                                            //
    //                   Step 0                   //
    //                                            //
    //********************************************//

    // Rotate matrix if width < height
    let rotated = width < height;

    let mut matrix = if rotated {
        let mut result = vec![0; width * height];
        for i in 0..height {
            for j in 0..width {
                result[index!(height, width - 1 - j, i)] = matrix[index!(width, i, j)];
            } 
        }
        result
    } else {
        Vec::from(matrix)
    };
    
    // Swap dimensions if rotated
    let (w, h) = if rotated { (height, width) } else { (width, height) };

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

        // Cover each column with a starred zero. If all columns are covered, we're done.
        if verify {
            stars.ones().for_each(|k| col_cover.insert(k % w));
            if col_cover.count_ones(..) == h {
                let assign = stars.ones().map(|k| k % w);

                // Rotate results back
                if rotated {
                    let mut result = vec![None; w];
                    assign.enumerate().for_each(|(i, j)| result[j] = Some(h - i - 1));
                    return result
                } else {
                    return assign.map(|j| Some(j)).collect()
                }
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
                    let k = index!(w, i, j);
                    if  row_cover[i] { matrix[k] += min }
                    if !col_cover[j] { matrix[k] -= min }
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
    fn test_basic_0x0() {
        let matrix = Vec::new(); 
        assert_eq!(
            minimize(&matrix, 0, 0),
            Vec::new()
        )
    }

    #[test]
    fn test_basic_1x1() {
        let matrix = vec![1];
        assert_eq!(
            minimize(&matrix, 1, 1),
            vec![Some(0)]
        );
    }

    #[test]
    fn test_basic_1x2() {
        let matrix = vec![
            1, 2
        ];
        assert_eq!(
            minimize(&matrix, 1, 2),
            vec![Some(0)]
        );
    }

    #[test]
    fn test_basic_2x1() {
        let matrix = vec![
            1,
            2
        ];
        assert_eq!(
            minimize(&matrix, 2, 1),
            vec![Some(0), None]
        );
    }

    #[test]
    fn test_basic_2x2() {
        let matrix = vec![
            1, 2,
            2, 1,
        ];
        assert_eq!(
            minimize(&matrix, 2, 2),
            vec![Some(0), Some(1)]
        );
    }

    // From http://www.math.harvard.edu/archive/20_spring_05/handouts/assignment_overheads.pdf
    #[test]
    fn test_sales_3x3() {
        let matrix = vec![
            250, 400, 350,
            400, 600, 350,
            200, 400, 250,
        ];
        assert_eq!(
            minimize(&matrix, 3, 3),
            vec![Some(1), Some(2), Some(0)]
        );
    }

    // From https://brilliant.org/wiki/hungarian-matching/
    #[test]
    fn test_party_3x3() {
        let matrix = vec![
            108, 125, 150,
            150, 135, 175,
            122, 148, 250,
        ];
        assert_eq!(
            minimize(&matrix, 3, 3),
            vec![Some(2), Some(1), Some(0)]
        );
    }

    // From https://en.wikipedia.org/wiki/Hungarian_algorithm#Matrix_interpretation
    #[test]
    fn test_wikipedia_4x4() {
        let matrix = vec![
            0, 1, 2, 3,
            4, 5, 6, 0,
            0, 2, 4, 5,
            3, 0, 0, 9,
        ];
        assert_eq!(
            minimize(&matrix, 4, 4),
            vec![Some(1), Some(3), Some(0), Some(2)]
        );
    }

    // From https://www.wikihow.com/Use-the-Hungarian-Algorithm
    #[test]
    fn test_wikihow_5x5() {
        let matrix = vec![
            10, 19, 8, 15, 19,
            10, 18, 7, 17, 19,
            13, 16, 9, 14, 19,
            12, 19, 8, 18, 19,
            14, 17, 10, 19, 19
        ];
        assert_eq!(
            minimize(&matrix, 5, 5),
            vec![Some(0), Some(2), Some(3), Some(4), Some(1)]
        );
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_5x5() {
        let matrix = vec![
            12,  9, 27, 10, 23,
             7, 13, 13, 30, 19,
            25, 18, 26, 11, 26,
             9, 28, 26, 23, 13,
            16, 16, 24,  6,  9,
        ];
        assert_eq!(
            51,
            minimize(&matrix, 5, 5)
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| v.map(|j| matrix[index!(5, i, j)]))
                .sum::<u64>()
        );
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_10x10() {
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
        assert_eq!(
            66,
            minimize(&matrix, 10, 10)
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| v.map(|j| matrix[index!(10, i, j)]))
                .sum::<u64>()
        );
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_20x20() {
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
        assert_eq!(
            22,
            minimize(&matrix, 20, 20)
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| v.map(|j| matrix[index!(20, i, j)]))
                .sum::<u64>()
        );
    }

    // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    #[test]
    fn test_stack_overflow_4x4_zeros() {
        let matrix = vec![
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 1, 2,
            0, 0, 3, 4,
        ];
        assert_eq!(
            minimize(&matrix, 4, 4),
            vec![Some(3), Some(2), Some(1), Some(0)]
        );
    }

    // From https://stackoverflow.com/questions/46803600/hungarian-algorithm-wikipedia-method-doesnt-work-for-this-example
    #[test]
    fn test_stack_overflow_4x4() {
        let matrix = vec![
            35,  0,  0,  0,
            0 , 30,  0,  5,
            55,  5,  0, 10,
            0 , 45, 30, 45,
        ];
        assert_eq!(
            5,
            minimize(&matrix, 4, 4)
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| v.map(|j| matrix[index!(4, i, j)]))
                .sum::<u64>()
        );
        assert_eq!(
            minimize(&matrix, 4, 4),
            vec![Some(1), Some(3), Some(2), Some(0)]
        );
    }

    // From https://stackoverflow.com/questions/17419595/hungarian-kuhn-munkres-algorithm-oddity
    #[test]
    fn test_stack_overflow_5x5() {
        let matrix = vec![
            0, 7, 0, 0, 0,
            0, 8, 0, 0, 6,
            5, 0, 7, 3, 4,
            5, 0, 5, 9, 3,
            0, 4, 0, 0, 9,
        ];
        assert_eq!(
            3,
            minimize(&matrix, 5, 5)
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| v.map(|j| matrix[index!(5, i, j)]))
                .sum::<u64>()
        );
        assert_eq!(
            minimize(&matrix, 5, 5), 
            vec![Some(4), Some(2), Some(3), Some(1), Some(0)]
        );
    }

    // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    #[test]
    fn test_stack_overflow_6x6() {
        let matrix = vec![
            2, 1, 0, 0, 0, 3,
            2, 0, 4, 5, 2, 7,
            0, 7, 0, 0, 0, 5,
            3, 2, 3, 1, 2, 0,
            0, 0, 6, 3, 3, 5,
            3, 4, 5, 2, 0, 3,
        ];
        assert_eq!(
            minimize(&matrix, 6, 6),
            vec![Some(3), Some(1), Some(2), Some(5), Some(0), Some(4)]
        );
    }

    // From https://stackoverflow.com/questions/26893961/cannot-solve-hungarian-algorithm
    #[test]
    fn test_stack_overflow_14x11() {
        let matrix = vec![
             0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
            53, 207, 256, 207, 231, 348, 348, 348, 231, 244, 244,
           240,  33,  67,  33,  56, 133, 133, 133,  56,  33,  33,
           460, 107, 200, 107, 122, 324, 324, 324, 122,  33,  33,
           167, 340, 396, 340, 422, 567, 567, 567, 422, 442, 442,
           167, 367, 307, 367, 433, 336, 336, 336, 433, 158, 158,
           160,  20,  37,  20,  31,  70,  70,  70,  31,  22,  22,
           200, 307, 393, 307, 222, 364, 364, 364, 222, 286, 286,
           33 , 153, 152, 153, 228, 252, 252, 252, 228,  78,  78,
           93 , 140, 185, 140,  58, 118, 118, 118,  58,  44,  44,
           0  ,   7,  22,   7,  19,  58,  58,  58,  19,   0,   0,
           67 , 153, 241, 153, 128, 297, 297, 297, 128,  39,  39,
           73 , 253, 389, 253, 253, 539, 539, 539, 253,  36,  36,
           173, 267, 270, 267, 322, 352, 352, 352, 322, 231, 231,
        ];
        assert_eq!(
            828,
            minimize(&matrix, 14, 11)
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| v.map(|j| matrix[index!(11, i, j)]))
                .sum::<u64>()
        );
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

            let expected = (0..max).map(|i| Some(i)).rev().collect::<Vec<_>>();
            assert_eq!(minimize(&matrix, max, max), expected);
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

            let expected = (0..max).map(|i| Some(i)).rev().collect::<Vec<_>>();
            assert_eq!(minimize(&matrix, max, max), expected);
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

        let expected = (0..max).map(|i| Some(i)).rev().collect::<Vec<_>>();
        assert_eq!(minimize(&matrix, max, max), expected);
    }
}
