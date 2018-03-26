extern crate fixedbitset;
extern crate num_traits;

use fixedbitset::FixedBitSet;
use num_traits::{PrimInt, NumAssign, NumOps};

use std::cmp::Ordering;
use std::fmt::Debug;

use std::collections::BinaryHeap;

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

struct Edge<N: PrimInt> {
    pub i: usize,
    pub j: usize,
    pub c: N,
}

impl<N: PrimInt> PartialEq for Edge<N> {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j 
    }
}

impl<N: PrimInt> Eq for Edge<N> {}

impl<N: PrimInt> PartialOrd for Edge<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<N: PrimInt> Ord for Edge<N> {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.c > other.c {
            Ordering::Less 
        } else if self.c < other.c {
            Ordering::Greater 
        } else {
            Ordering::Equal
        }
    }
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
/// // extern crate hungarian;
///
/// // use hungarian::minimize;
///
/// // fn main() {
///
/// //     // Square matrix
///
/// //     let width = 3;
/// //     let height = 3;
/// //     let matrix = vec![
/// //         1, 2, 1,
/// //         4, 5, 6,
/// //         7, 8, 9,
/// //     ];
///
/// //     let assignment = minimize(&matrix, height, width);
///
/// //     assert_eq!(&assignment, &vec![Some(2), Some(1), Some(0)]);
///
/// //     let cost: i32 = assignment.iter()
/// //         .enumerate()
/// //         .filter_map(|(i, &a)| {
/// //             a.map(|j| matrix[i*width + j])
/// //         })
/// //         .sum();
///
/// //     assert_eq!(cost, 13);
///
/// //     // Rectangular matrix (height < width)
///
/// //     let height = 2;
/// //     let width = 3;
/// //     let matrix = vec![
/// //         1, 0, 5,
/// //         2, 3, 1,
/// //     ];
///
/// //     let assignment = minimize(&matrix, height, width);
///
/// //     assert_eq!(&assignment, &vec![Some(1), Some(2)]);
///
/// //     let cost: i32 = assignment.iter()
/// //         .enumerate()
/// //         .filter_map(|(i, &a)| {
/// //             a.map(|j| matrix[i*width + j])
/// //         })
/// //         .sum();
///
/// //     assert_eq!(cost, 1);
///
/// //     // Rectangular matrix (width < height)
///
/// //     let height = 3;
/// //     let width = 2;
/// //     let matrix = vec![
/// //         5, 5,
/// //         1, 0,
/// //         2, 3,
/// //     ];
///
/// //     let assignment = minimize(&matrix, height, width);
///
/// //     assert_eq!(&assignment, &vec![None, Some(1), Some(0)]);
///
/// //     let cost: i32 = assignment.iter()
/// //         .enumerate()
/// //         .filter_map(|(i, &a)| {
/// //             a.map(|j| matrix[i*width + j])
/// //         })
/// //         .sum();
///
/// //     assert_eq!(cost, 2);
///
/// // }
/// ```
///
/// [1]: http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
///
pub fn minimize<N: NumAssign + PrimInt + Debug>(matrix: &[N], h: usize, w: usize) -> Vec<Option<usize>> {

    let mut p = vec![N::zero(); h];
    let mut q = vec![N::zero(); w];

    let mut mark_i: Vec<Option<usize>> = vec![None; h];
    let mut mark_j: Vec<Option<usize>> = vec![None; w];

    let mut tree_i = FixedBitSet::with_capacity(h);
    let mut tree_j = FixedBitSet::with_capacity(w);

    let mut path: Vec<Option<usize>> = vec![None; w];
    let mut queue: BinaryHeap<Edge<N>> = BinaryHeap::new();

    let mut si = vec![0; h];
    let mut sj = vec![0; w];
    let mut ds = Vec::new();

    // While not a perfect matching
    for _ in 0..h {

        tree_i.clear();
        tree_j.clear();

        queue.clear();
        ds.clear();
        ds.push(N::zero());
        let mut step = 0;


        // Initialize T to be the set of free vertices in X
        // for i in 0..h {
        //     if mark_i[i].is_none() {

        //         tree_i.insert(i);
        //         si[i] = 0;

        //         for j in 0..w {
                    
        //             sj[j] = 0;

        //             let c = matrix[i*w + j] - p[i] - q[j];

        //             queue.push(Edge { i, j, c });
                
        //         }

        //     }
        // }  

        // Initialize T to be a free vertex in X
        let i = (0..h).find(|&i| mark_i[i].is_none()).unwrap();
        tree_i.insert(i);
        si[i] = 0;

        // Update queue with all edges from X intersect T to Y \ T
        for j in 0..w {
            sj[j]  = 0;
            let c = matrix[i*w + j] - p[i] - q[j];
            queue.push(Edge { i, j, c });
        }

        for _ in 0..w {

            step += 1;
            println!("Step {}", step);

            let mut edge = queue.pop().unwrap();

            while on!(tree_j, edge.j) { edge = queue.pop().unwrap(); }

            let (i, j) = (edge.i, edge.j);

            // let pds = ds[si[i]].clone();
            // ds.push(pds + ((edge.c - pds) / (N::one() + N::one())));
           
            // p0 + ds(u) - q0 - d(s-1) + d(s - 1)?
            ds.push(edge.c);

            // Backtrack
            path[j] = Some(i);
            
            // j is a free vertex
            if let None = mark_j[j] {
                
                let mut y = Some(j);

                // Toggle edges
                while y.is_some() {

                    let x = path[y.unwrap()].unwrap();
                    let previous = mark_i[x];
                    mark_j[y.unwrap()] = Some(x);
                    mark_i[x] = y;
                    y = previous;

                }

                let mut sum = N::zero();
                for i in tree_i.ones() {
                    sum += p[i];
                }

                for j in tree_j.ones() {
                    sum += q[j];
                }

                println!("{:?}", mark_i);
                println!("{:?}", ds);
                println!("{:?}", p);
                println!("{:?}", q);
                println!("{:?}", sum);

                for i in tree_i.ones() {
                    p[i] += ds[ds.len() - 1] - ds[si[i]];
                    sum += p[i];
                }

                for j in tree_j.ones() {
                    q[j] += ds[sj[j]] - ds[ds.len() - 1];
                    sum += q[j];
                }

                // Finish phase
                break;
            }

            // X -> Y -> X -> Y

            // Tree growing step
            
            // Identify an edge x', y in M and add y, x' to T
            else if let Some(i) = mark_j[j] {

                // How do we know what edges are in M?

                let step = ds.len() - 1;

                tree_i.insert(i);
                tree_j.insert(j);

                si[i] = step;

                for j in 0..w {

                    if off!(tree_j, j)  {
                    
                        sj[j] = step;

                        let c = matrix[i*w + j] + ds[step] - q[j] - p[i];

                        queue.push(Edge { i, j, c });

                    }
                }
            }
        }

    }

    return mark_i;    
}

#[cfg(test)]
mod tests {

    /// Internal macro for converting from a 2D index to a 1D index.
    macro_rules! index {
        ($w:expr, $i:expr, $j:expr) => (($w*$i) + $j)
    }

    use minimize;

    // #[test]
    // fn test_basic_0x0() {
    //     let matrix: Vec<i32> = Vec::new();
    //     assert_eq!(
    //         minimize(&matrix, 0, 0),
    //         Vec::new()
    //     )
    // }

    // #[test]
    // fn test_basic_1x1() {
    //     let matrix = vec![1];
    //     assert_eq!(
    //         minimize(&matrix, 1, 1),
    //         vec![Some(0)]
    //     );
    // }

    // // #[test]
    // // fn test_basic_1x2() {
    // //     let matrix = vec![
    // //         1, 2
    // //     ];
    // //     assert_eq!(
    // //         minimize(&matrix, 1, 2),
    // //         vec![Some(0)]
    // //     );
    // // }

    // // #[test]
    // // fn test_basic_2x1() {
    // //     let matrix = vec![
    // //         1,
    // //         2
    // //     ];
    // //     assert_eq!(
    // //         minimize(&matrix, 2, 1),
    // //         vec![Some(0), None]
    // //     );
    // // }

    // #[test]
    // fn test_basic_2x2() {
    //     let matrix = vec![
    //         1, 2,
    //         2, 1,
    //     ];
    //     assert_eq!(
    //         minimize(&matrix, 2, 2),
    //         vec![Some(0), Some(1)]
    //     );
    // }

    // // From http://www.math.harvard.edu/archive/20_spring_05/handouts/assignment_overheads.pdf
    // #[test]
    // fn test_sales_3x3() {
    //     let matrix = vec![
    //         250, 400, 350,
    //         400, 600, 350,
    //         200, 400, 250,
    //     ];
    //     assert_eq!(
    //         minimize(&matrix, 3, 3),
    //         vec![Some(1), Some(2), Some(0)]
    //     );
    // }

    // // From https://brilliant.org/wiki/hungarian-matching/
    // #[test]
    // fn test_party_3x3() {
    //     let matrix = vec![
    //         108, 125, 150,
    //         150, 135, 175,
    //         122, 148, 250,
    //     ];
    //     assert_eq!(
    //         minimize(&matrix, 3, 3),
    //         vec![Some(2), Some(1), Some(0)]
    //     );
    // }

    // // From https://en.wikipedia.org/wiki/Hungarian_algorithm#Matrix_interpretation
    // #[test]
    // fn test_wikipedia_4x4() {
    //     let matrix = vec![
    //         0, 1, 2, 3,
    //         4, 5, 6, 0,
    //         0, 2, 4, 5,
    //         3, 0, 0, 9,
    //     ];
    //     assert_eq!(
    //         minimize(&matrix, 4, 4),
    //         vec![Some(1), Some(3), Some(0), Some(2)]
    //     );
    // }

    // // From https://www.wikihow.com/Use-the-Hungarian-Algorithm
    // #[test]
    // fn test_wikihow_5x5() {
    //     let matrix = vec![
    //         10, 19, 8, 15, 19,
    //         10, 18, 7, 17, 19,
    //         13, 16, 9, 14, 19,
    //         12, 19, 8, 18, 19,
    //         14, 17, 10, 19, 19
    //     ];
    //     assert_eq!(
    //         minimize(&matrix, 5, 5),
    //         vec![Some(0), Some(2), Some(3), Some(4), Some(1)]
    //     );
    // }

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
                .sum::<i32>()
        );
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    // #[test]
    // fn test_python_10x10() {
    //     let matrix = vec![
    //         37, 34, 29, 26, 19,  8,  9, 23, 19, 29,
    //          9, 28, 20,  8, 18, 20, 14, 33, 23, 14,
    //         15, 26, 12, 28,  6, 17,  9, 13, 21,  7,
    //          2,  8, 38, 36, 39,  5, 36,  2, 38, 27,
    //         30,  3, 33, 16, 21, 39,  7, 23, 28, 36,
    //          7,  5, 19, 22, 36, 36, 24, 19, 30,  2,
    //         34, 20, 13, 36, 12, 33,  9, 10, 23,  5,
    //          7, 37, 22, 39, 33, 39, 10,  3, 13, 26,
    //         21, 25, 23, 39, 31, 37, 32, 33, 38,  1,
    //         17, 34, 40, 10, 29, 37, 40,  3, 25,  3,
    //     ];
    //     assert_eq!(
    //         66,
    //         minimize(&matrix, 10, 10)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(10, i, j)]))
    //             .sum::<i32>()
    //     );
    // }

    // // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    // #[test]
    // fn test_python_20x20() {
    //     let matrix = vec![
    //         5,  4,  3,  9,  8,  9,  3,  5,  6,  9,  4, 10,  3,  5,  6,  6,  1,  8, 10,  2,
    //         10, 9,  9,  2,  8,  3,  9,  9, 10,  1,  7, 10,  8,  4,  2,  1,  4,  8,  4,  8,
    //         10, 4,  4,  3,  1,  3,  5, 10,  6,  8,  6,  8,  4, 10,  7,  2,  4,  5,  1,  8,
    //         2,  1,  4,  2,  3,  9,  3,  4,  7,  3,  4,  1,  3,  2,  9,  8,  6,  5,  7,  8,
    //         3,  4,  4,  1,  4, 10,  1,  2,  6,  4,  5, 10,  2,  2,  3,  9, 10,  9,  9, 10,
    //         1, 10,  1,  8,  1,  3,  1,  7,  1,  1,  2,  1,  2,  6,  3,  3,  4,  4,  8,  6,
    //         1,  8,  7, 10, 10,  3,  4,  6,  1,  6,  6,  4,  9,  6,  9,  6,  4,  5,  4,  7,
    //         8, 10,  3,  9,  4,  9,  3,  3,  4,  6,  4,  2,  6,  7,  7,  4,  4,  3,  4,  7,
    //         1,  3,  8,  2,  6,  9,  2,  7,  4,  8, 10,  8, 10,  5,  1,  3, 10, 10,  2,  9,
    //         2,  4,  1,  9,  2,  9,  7,  8,  2,  1,  4, 10,  5,  2,  7,  6,  5,  7,  2,  6,
    //         4,  5,  1,  4,  2,  3,  3,  4,  1,  8,  8,  2,  6,  9,  5,  9,  6,  3,  9,  3,
    //         3,  1,  1,  8,  6,  8,  8,  7,  9,  3,  2,  1,  8,  2,  4,  7,  3,  1,  2,  4,
    //         5,  9,  8,  6, 10,  4, 10,  3,  4, 10, 10, 10,  1,  7,  8,  8,  7,  7,  8,  8,
    //         1,  4,  6,  1,  6,  1,  2, 10,  5, 10,  2,  6,  2,  4,  5,  5,  3,  5,  1,  5,
    //         5,  6,  9, 10,  6,  6, 10,  6,  4,  1,  5,  3,  9,  5,  2, 10,  9,  9,  5,  1,
    //         10, 9,  4,  6,  9,  5,  3,  7, 10,  1,  6,  8,  1,  1, 10,  9,  5,  7,  7,  5,
    //         2,  6,  6,  6,  6,  2,  9,  4,  7,  5,  3,  2, 10,  3,  4,  5, 10,  9,  1,  7,
    //         5,  2,  4,  9,  8,  4,  8,  2,  4,  1,  3,  7,  6,  8,  1,  6,  8,  8, 10, 10,
    //         9,  6,  3,  1,  8,  5,  7,  8,  7,  2,  1,  8,  2,  8,  3,  7,  4,  8,  7,  7,
    //         8,  4,  4,  9,  7, 10,  6,  2,  1,  5,  8,  5,  1,  1,  1,  9,  1,  3,  5,  3,
    //     ];
    //     assert_eq!(
    //         22,
    //         minimize(&matrix, 20, 20)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(20, i, j)]))
    //             .sum::<i32>()
    //     );
    // }

    // // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    // #[test]
    // fn test_stack_overflow_4x4_zeros() {
    //     let matrix = vec![
    //         0, 0, 0, 0,
    //         0, 0, 0, 0,
    //         0, 0, 1, 2,
    //         0, 0, 3, 4,
    //     ];
    //     assert_eq!(
    //         minimize(&matrix, 4, 4),
    //         vec![Some(3), Some(2), Some(1), Some(0)]
    //     );
    // }

    // // From https://stackoverflow.com/questions/46803600/hungarian-algorithm-wikipedia-method-doesnt-work-for-this-example
    // #[test]
    // fn test_stack_overflow_4x4() {
    //     let matrix = vec![
    //         35,  0,  0,  0,
    //         0 , 30,  0,  5,
    //         55,  5,  0, 10,
    //         0 , 45, 30, 45,
    //     ];
    //     assert_eq!(
    //         5,
    //         minimize(&matrix, 4, 4)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(4, i, j)]))
    //             .sum::<i32>()
    //     );
    //     assert_eq!(
    //         minimize(&matrix, 4, 4),
    //         vec![Some(1), Some(3), Some(2), Some(0)]
    //     );
    // }

    // // From https://stackoverflow.com/questions/17419595/hungarian-kuhn-munkres-algorithm-oddity
    // #[test]
    // fn test_stack_overflow_5x5() {
    //     let matrix = vec![
    //         0, 7, 0, 0, 0,
    //         0, 8, 0, 0, 6,
    //         5, 0, 7, 3, 4,
    //         5, 0, 5, 9, 3,
    //         0, 4, 0, 0, 9,
    //     ];
    //     assert_eq!(
    //         3,
    //         minimize(&matrix, 5, 5)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(5, i, j)]))
    //             .sum::<i32>()
    //     );
    //     assert_eq!(
    //         minimize(&matrix, 5, 5),
    //         vec![Some(4), Some(2), Some(3), Some(1), Some(0)]
    //     );
    // }

    // // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    // #[test]
    // fn test_stack_overflow_6x6() {
    //     let matrix = vec![
    //         2, 1, 0, 0, 0, 3,
    //         2, 0, 4, 5, 2, 7,
    //         0, 7, 0, 0, 0, 5,
    //         3, 2, 3, 1, 2, 0,
    //         0, 0, 6, 3, 3, 5,
    //         3, 4, 5, 2, 0, 3,
    //     ];
    //     assert_eq!(
    //         minimize(&matrix, 6, 6),
    //         vec![Some(3), Some(1), Some(2), Some(5), Some(0), Some(4)]
    //     );
    // }

    // // From https://stackoverflow.com/questions/26893961/cannot-solve-hungarian-algorithm
    // #[test]
    // fn test_stack_overflow_14x11() {
    //     let matrix = vec![
    //          0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    //         53, 207, 256, 207, 231, 348, 348, 348, 231, 244, 244,
    //        240,  33,  67,  33,  56, 133, 133, 133,  56,  33,  33,
    //        460, 107, 200, 107, 122, 324, 324, 324, 122,  33,  33,
    //        167, 340, 396, 340, 422, 567, 567, 567, 422, 442, 442,
    //        167, 367, 307, 367, 433, 336, 336, 336, 433, 158, 158,
    //        160,  20,  37,  20,  31,  70,  70,  70,  31,  22,  22,
    //        200, 307, 393, 307, 222, 364, 364, 364, 222, 286, 286,
    //        33 , 153, 152, 153, 228, 252, 252, 252, 228,  78,  78,
    //        93 , 140, 185, 140,  58, 118, 118, 118,  58,  44,  44,
    //        0  ,   7,  22,   7,  19,  58,  58,  58,  19,   0,   0,
    //        67 , 153, 241, 153, 128, 297, 297, 297, 128,  39,  39,
    //        73 , 253, 389, 253, 253, 539, 539, 539, 253,  36,  36,
    //        173, 267, 270, 267, 322, 352, 352, 352, 322, 231, 231,
    //     ];
    //     assert_eq!(
    //         828,
    //         minimize(&matrix, 14, 11)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(11, i, j)]))
    //             .sum::<i32>()
    //     );
    // }

    // // From https://github.com/bmc/munkres/blob/master/munkres.py
    // #[test]
    // fn test_rectangle_3x4() {
    //     let matrix = vec![
    //         400, 150, 400, 1,
    //         400, 450, 600, 2,
    //         300, 225, 300, 3,
    //     ];
    //     assert_eq!(
    //         452,
    //         minimize(&matrix, 3, 4)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(4, i, j)]))
    //             .sum::<i32>()
    //     );
    //     assert_eq!(
    //         minimize(&matrix, 3, 4),
    //         vec![Some(1), Some(3), Some(0)]
    //     );
    // }

    // // Modified from http://www.hungarianalgorithm.com/examplehungarianalgorithm.php
    // #[test]
    // fn test_rectangle_4x5() {
    //     let matrix = vec![
    //         82, 83, 69, 92, 100,
    //         77, 37, 49, 92, 195,
    //         11, 69,  5, 86,  93,
    //          8,  9, 98, 23, 106,
    //     ];
    //     assert_eq!(
    //         140,
    //         minimize(&matrix, 4, 5)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(5, i, j)]))
    //             .sum::<i32>()
    //     );
    //     assert_eq!(
    //         minimize(&matrix, 4, 5),
    //         vec![Some(2), Some(1), Some(0), Some(3)]
    //     );
    // }

    // // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    // #[test]
    // fn test_rectangle_5x4() {
    //     let matrix = vec![
    //         34, 26, 17, 12,
    //         43, 43, 36, 10,
    //         97, 47, 66, 34,
    //         52, 42, 19, 36,
    //         15, 93, 55, 80
    //     ];
    //     assert_eq!(
    //         70,
    //         minimize(&matrix, 5, 4)
    //             .iter()
    //             .enumerate()
    //             .filter_map(|(i, &v)| v.map(|j| matrix[index!(4, i, j)]))
    //             .sum::<i32>()
    //     );
    //     assert_eq!(
    //         minimize(&matrix, 5, 4),
    //         vec![Some(1), Some(3), None, Some(2), Some(0)]
    //     );
    // }

    // #[test]
    // fn test_rectangle_nx1() {
    //     let max = 100;
    //     for n in 1..max {
    //         let matrix = (0..n as i32).rev().collect::<Vec<_>>();
    //         let mut expected = vec![None; n];
    //         expected[n - 1] = Some(0);
    //         assert_eq!(minimize(&matrix, n, 1), expected);
    //     }
    // }

    // #[test]
    // fn test_rectangle_1xn() {
    //     let max = 100;
    //     for n in 1..max {
    //         let matrix = (0..n as i32).rev().collect::<Vec<_>>();
    //         let expected = vec![Some(n - 1)];
    //         assert_eq!(minimize(&matrix, 1, n), expected);
    //     }
    // }

    // #[test]
    // fn test_stress() {
    //     for max in 1..100 {
    //         let mut matrix = vec![0; max * max];
    //         let mut n: i32 = 0;

    //         for i in 0..max {
    //             for j in 0..max {
    //                 matrix[index!(max, i, j)] = n;
    //                 n += 1;
    //             }
    //         }

    //         let expected = (0..max).map(|i| Some(i)).rev().collect::<Vec<_>>();
    //         assert_eq!(minimize(&matrix, max, max), expected);
    //     }
    // }

    // #[test]
    // fn test_worst_case() {
    //     for max in 1..50 {
    //         let mut matrix = vec![0; max * max];

    //         for i in 0..max {
    //             for j in 0..max {
    //                 matrix[index!(max, i, j)] = ((i + 1)*(j + 1)) as i32;
    //             }
    //         }

    //         let expected = (0..max).map(|i| Some(i)).rev().collect::<Vec<_>>();
    //         assert_eq!(minimize(&matrix, max, max), expected);
    //     }
    // }

    // #[test]
    // fn test_large() {
    //     let max = 1000;
    //     let mut matrix = vec![0; max * max];
    //     let mut n: i32 = 0;

    //     for i in 0..max {
    //         for j in 0..max {
    //             matrix[index!(max, i, j)] = n;
    //             n += 1;
    //         }
    //     }

    //     let expected = (0..max).map(|i| Some(i)).rev().collect::<Vec<_>>();
    //     assert_eq!(minimize(&matrix, max, max), expected);
    // }
}
