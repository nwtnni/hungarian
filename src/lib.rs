extern crate fixedbitset;
extern crate num_traits;
extern crate ndarray;

use fixedbitset::FixedBitSet;
use num_traits::{PrimInt, NumAssign};
use ndarray::prelude::Array2;

/// Internal macro for indexing an Array2 without the bounds check
macro_rules! get {
    ($m:expr, $i:expr, $j:expr) => (unsafe { *$m.uget(($i, $j)) })
}

/// Internal macro for mutating an Array2 without the bounds check
macro_rules! set {
    ($m:expr, $i:expr, $j:expr, $v: expr) => (unsafe { *$m.uget_mut(($i, $j)) = $v; })
}

/// Internal macro for querying a FixedBitSet.
/// Syntactic sugar for `s[i]`, but without the runtime overhead of the Index trait.
macro_rules! on {
    ($s:expr, $i:expr) => ($s.contains($i))
}

/// Internal macro for querying a FixedBitSet.
/// Syntactic sugar for `!s[i]`, but without the runtime overhead of the Index trait.
macro_rules! off {
    ($s:expr, $i:expr) => (!$s.contains($i))
}

/// Implementation of the Hungarian / Munkres Assignment Algorithm.
///
/// Given a rectangular cost matrix, this algorithm finds a maximal matching such
/// that the total cost is minimized. [Follows the general outline explained here.][1]
/// This implementation only works on integer costs (since checking if a float is 0 is
/// not a great idea).
///
/// The function will automatically clamp all costs in the
/// provided matrix to be greater or equal to zero, and as a result, won't correctly
/// calculate the minimum assignment for a matrix with negative entries.
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
/// This function uses unsafe array indexing directly in order to minimize,
/// runtime costs, and will eventually panic with out-of-bounds if passed
/// incorrect `width` or `height` arguments.
///
/// If given correct arguments, there is no way we can index out of bounds,
/// even in unsafe blocks: we only ever iterate between 0 and width/height.
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
pub fn minimize<N: NumAssign + PrimInt>(matrix: &[N], height: usize, width: usize) -> Vec<Option<usize>> {

    // No possible assignment
    if height == 0 || width == 0 { return Vec::new() }

    //********************************************//
    //                                            //
    //                   Step 0                   //
    //                                            //
    //********************************************//

    // Rotate matrix if width < height
    let rotated = width < height;
    let (w, h) = if rotated { (height, width) } else { (width, height) };
    let mut m = Array2::zeros((h, w));

    // Clamp matrix to be positive and rotate if necessary
    for i in 0..height {
        for j in 0..width {
            let cost = matrix[width * i + j];
            if cost < N::zero() {
                continue
            } else if rotated {
                set!(m, width - 1 - j, i, cost)
            } else {
                set!(m, i, j, cost)
            }
        }
    }

    // The set of starred zero entries
    let mut stars = Array2::from_elem((h, w), false);

    // The set of primed zero entries
    let mut primes = Array2::from_elem((h, w), false);

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
    for mut row in m.genrows_mut() {
        let min = row.iter().min().unwrap().clone();
        row.map_inplace(|v| *v -= min);
    }

    //********************************************//
    //                                            //
    //                   Step 2                   //
    //                                            //
    //********************************************//

    // Find a zero (Z):
    // - If there is no starred zero in its row or column, then star it.
    // - Use col_cover to keep track of stars.
    for i in 0..h {
        for j in 0..w {
            if on!(col_cover, j) { continue }
            if get!(m, i, j).is_zero() {
                set!(stars, i, j, true);
                col_cover.insert(j);
                break
            }
        }
    }

    // Reset cover
    col_cover.clear();
    let mut verify = true;

    loop {

        if verify {

            //********************************************//
            //                                            //
            //                   Step 3                   //
            //                                            //
            //********************************************//

            // Cover each column with a starred zero.
            stars.gencolumns()
                .into_iter()
                .enumerate()
                .for_each(|(j, col)| {
                    if col.iter().any(|&s| s) { col_cover.insert(j) }
                });

            // If the number of starred zeros equals the number of rows, we're done.
            if col_cover.count_ones(..) == h {

                let assign = stars.genrows().into_iter().map(|r| {
                    r.iter().enumerate()
                        .find(|&(_, &v)| v)
                        .map(|(i, _)| i)
                        .unwrap()
                });

                // Rotate results back if necessary
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

        let mut uncovered = None;

        // Find an uncovered zero and prime it
        'outer : for i in 0..h {
            if on!(row_cover, i) { continue }
            for j in 0..w {
                if on!(col_cover, j) { continue }
                if get!(m, i, j).is_zero() {
                    uncovered = Some((i, j));
                    set!(primes, i, j, true);
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
            let mut min = N::max_value();
            for i in 0..h {
                if on!(row_cover, i) { continue }
                for j in 0..w {
                    if on!(col_cover, j) { continue }
                    let value = get!(m, i, j);
                    min = if value < min { value } else { min };
                }
            }

            // Add minimum to covered rows
            for i in (0..h).filter(|&i| on!(row_cover, i)) {
                m.row_mut(i).map_inplace(|c| *c += min)
            }

            // Subtract minimum from uncovered columns
            for j in (0..w).filter(|&j| off!(col_cover, j)) {
                m.column_mut(j).map_inplace(|c| *c -= min)
            }

            // Return to [Step 4]
            // - Skip rest of this loop
            // - Set `verify` to false to skip [Step 3]
            verify = false;
            continue
        }

        let (i, j) = uncovered.unwrap();

        // If there's a starred zero in the same row
        // - Cover row of uncovered zero from [Step 4]
        // - Uncover column of starred zero
        // - Repeat [Step 4]
        if let Some(j) = (0..w).find(|&j| get!(stars, i, j)) {
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
            let next_star = (0..h).find(|&i| get!(stars, i, j));

            if let None = next_star { break }
            let i = next_star.unwrap();
            path.push((i, j));

            // Find primed zero in same row
            // Guaranteed to exist
            let j = (0..w).find(|&j| get!(primes, i, j)).unwrap();
            path.push((i, j));
        }

        // Unstar each starred zero
        // Star each primed zero
        for (i, j) in path {
            set!(stars, i, j, *primes.uget((i, j)));
        }

        // Reset cover
        row_cover.clear();
        col_cover.clear();

        // Erase primes and return to [Step 3]
        primes.map_inplace(|p| *p = false);
        verify = true;
    }
}
