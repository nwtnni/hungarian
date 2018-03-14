extern crate fnv;
extern crate fixedbitset;

use fnv::FnvHashMap;
use fixedbitset::FixedBitSet;

#[derive(Debug, Eq, PartialEq)]
enum Type { Star, Prime }

pub fn hungarian(mut matrix: Vec<Vec<u64>>) -> Vec<usize> {
    let rows = matrix.len();
    let cols = matrix[0].len();
    let target = if rows < cols { rows } else { cols };

    // Reduce rows
    for mut row in &mut matrix {
        let min = row.iter().min().unwrap().clone();
        row.iter_mut().for_each(|cost| *cost -= min);
    }

    let mut mask = FnvHashMap::default();
    let mut row_cover = vec![false; rows];
    let mut col_cover = vec![false; cols];

    // Star zeros
    for row in 0..rows {
        for col in 0..cols {
            if matrix[row][col] == 0
            && !(row_cover[row] || col_cover[col]) {
                mask.insert((row, col), Type::Star);
                row_cover[row] = true;
                col_cover[col] = true;
            }
        }
    }

    // Reset cover
    row_cover.iter_mut().for_each(|cov| *cov = false);
    col_cover.iter_mut().for_each(|cov| *cov = false);
    let mut verify = true;

    loop {

        // Count cover
        if verify {
            for row in 0..rows {
                for col in 0..cols {
                    if let Some(&Type::Star) = mask.get(&(row, col)) {
                        col_cover[col] = true;
                    }
                }
            }

            if col_cover.iter().filter(|&&cov| cov).count() == target {
                let mut result = mask.into_iter()
                    .filter(|&(_, ref t)| t == &Type::Star)
                    .map(|(key, _)| key)
                    .collect::<Vec<_>>();
                result.sort_by_key(|&(row, _)| row);
                return result.into_iter()
                    .map(|(_, col)| col)
                    .collect()
            }
        }

        // Find uncovered zero
        let mut uncovered = None;
        for row in 0..rows {
            if uncovered != None { break }
            for col in 0..cols {
                if matrix[row][col] == 0
                && mask.get(&(row, col)).is_none()
                && !row_cover[row]
                && !col_cover[col] {
                    uncovered = Some((row, col));
                    break
                }
            }
        }

        // Add and subtract minimum uncovered value
        if let None = uncovered {
            let mut min = u64::max_value();
            for row in 0..rows {
                for col in 0..cols {
                    if row_cover[row] || col_cover[col] { continue }
                    let value = matrix[row][col];
                    min = if value < min { value } else { min };
                }
            }

            for row in 0..rows {
                for col in 0..cols {
                    if row_cover[row] { matrix[row][col] += min }
                    if !col_cover[col] { matrix[row][col] -= min }
                }
            }

            verify = false;
            continue
        }

        let (row, col) = uncovered.unwrap();
        mask.insert((row, col), Type::Prime);

        let starred = (0..cols).filter(|&col| {
            mask.get(&(row, col)) == Some(&Type::Star)
            && matrix[row][col] == 0
        }).nth(0);

        if let Some(adj) = starred {
            row_cover[row] = true;
            col_cover[adj] = false;
            verify = false;
            continue
        }

        // Alternating path of Stars and Primes
        let mut path = vec![(row, col)];
        loop {
            let (_, prev_col) = path[path.len() - 1];
            let next_star = (0..rows).filter(|&row| {
                mask.get(&(row, prev_col)) == Some(&Type::Star)
            }).nth(0);

            if let None = next_star { break }
            let star_row = next_star.unwrap();
            path.push((star_row, prev_col));

            let prime_col = (0..cols).filter(|&col| {
                mask.get(&(star_row, col)) == Some(&Type::Prime)
            }).nth(0).unwrap();
            path.push((star_row, prime_col));
        }

        // Augment path
        for (row, col) in path {
            match mask.get(&(row, col)) {
                None => continue,
                Some(&Type::Star) => mask.remove(&(row, col)),
                Some(&Type::Prime) => mask.insert((row, col), Type::Star),
            };
        }

        // Reset cover
        row_cover.iter_mut().for_each(|cov| *cov = false);
        col_cover.iter_mut().for_each(|cov| *cov = false);

        // Erase primes
        mask.retain(|_, t| t != &mut Type::Prime);
        verify = true;
    }
}

#[cfg(test)]
mod tests {
    use hungarian;

    #[test]
    fn test_basic() {
        let matrix = vec![
            vec![1, 2, 2],
            vec![2, 1, 2],
            vec![2, 2, 1],
        ];
        assert_eq!(hungarian(matrix), vec![0, 1, 2])
    }

    #[test]
    fn test_increasing() {
        let matrix = vec![
            vec![1,   2,  3,  4],
            vec![5,   6,  7,  8],
            vec![9,  10, 11, 12],
            vec![13, 14, 15, 16],
        ];
        assert_eq!(hungarian(matrix), vec![3, 2, 1, 0])
    }

    // From http://www.math.harvard.edu/archive/20_spring_05/handouts/assignment_overheads.pdf
    #[test]
    fn test_sales_example() {
        let matrix = vec![
            vec![250, 400, 350],
            vec![400, 600, 350],
            vec![200, 400, 250],
        ];
        assert_eq!(hungarian(matrix), vec![1, 2, 0]);
    }

    // From https://brilliant.org/wiki/hungarian-matching/
    #[test]
    fn test_party_example() {
        let matrix = vec![
            vec![108, 125, 150],
            vec![150, 135, 175],
            vec![122, 148, 250],
        ];
        assert_eq!(hungarian(matrix), vec![2, 1, 0]);
    }

    // From https://en.wikipedia.org/wiki/Hungarian_algorithm#Matrix_interpretation
    #[test]
    fn test_wiki() {
        let matrix = vec![
            vec![0, 1, 2, 3],
            vec![4, 5, 6, 0],
            vec![0, 2, 4, 5],
            vec![3, 0, 0, 9],
        ];
        assert_eq!(hungarian(matrix), vec![1, 3, 0, 2]);
    }

    // From https://www.wikihow.com/Use-the-Hungarian-Algorithm
    #[test]
    fn test_bulldozer() {
        let matrix = vec![
            vec![90, 75, 75, 80],
            vec![35, 85, 55, 65],
            vec![125, 95, 90, 105],
            vec![45, 110, 95, 115],
        ];
        hungarian(matrix);
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_5() {
        let matrix = vec![
            vec![12, 9, 27, 10, 23],
            vec![7, 13, 13, 30, 19],
            vec![25, 18, 26, 11, 26],
            vec![9, 28, 26, 23, 13],
            vec![16, 16, 24, 6, 9],
        ];
        assert_eq!(hungarian(matrix.clone()).iter()
            .enumerate()
            .map(|(i, &j)| matrix[i][j])
            .sum::<u64>(), 51);
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_10() {
        let matrix = vec![ 
            vec![37, 34, 29, 26, 19, 8, 9, 23, 19, 29],
            vec![9, 28, 20, 8, 18, 20, 14, 33, 23, 14],
            vec![15, 26, 12, 28, 6, 17, 9, 13, 21, 7],
            vec![2, 8, 38, 36, 39, 5, 36, 2, 38, 27],
            vec![30, 3, 33, 16, 21, 39, 7, 23, 28, 36],
            vec![7, 5, 19, 22, 36, 36, 24, 19, 30, 2],
            vec![34, 20, 13, 36, 12, 33, 9, 10, 23, 5],
            vec![7, 37, 22, 39, 33, 39, 10, 3, 13, 26],
            vec![21, 25, 23, 39, 31, 37, 32, 33, 38, 1],
            vec![17, 34, 40, 10, 29, 37, 40, 3, 25, 3],
        ];
        assert_eq!(hungarian(matrix.clone()).iter()
            .enumerate()
            .map(|(i, &j)| matrix[i][j])
            .sum::<u64>(), 66);
    }

    // From https://github.com/bmc/munkres/blob/master/test/test_munkres.py
    #[test]
    fn test_python_20() {
        let matrix = vec![
            vec![5, 4, 3, 9, 8, 9, 3, 5, 6, 9, 4, 10, 3, 5, 6, 6, 1, 8, 10, 2],
            vec![10, 9, 9, 2, 8, 3, 9, 9, 10, 1, 7, 10, 8, 4, 2, 1, 4, 8, 4, 8],
            vec![10, 4, 4, 3, 1, 3, 5, 10, 6, 8, 6, 8, 4, 10, 7, 2, 4, 5, 1, 8],
            vec![2, 1, 4, 2, 3, 9, 3, 4, 7, 3, 4, 1, 3, 2, 9, 8, 6, 5, 7, 8],
            vec![3, 4, 4, 1, 4, 10, 1, 2, 6, 4, 5, 10, 2, 2, 3, 9, 10, 9, 9, 10],
            vec![1, 10, 1, 8, 1, 3, 1, 7, 1, 1, 2, 1, 2, 6, 3, 3, 4, 4, 8, 6],
            vec![1, 8, 7, 10, 10, 3, 4, 6, 1, 6, 6, 4, 9, 6, 9, 6, 4, 5, 4, 7],
            vec![8, 10, 3, 9, 4, 9, 3, 3, 4, 6, 4, 2, 6, 7, 7, 4, 4, 3, 4, 7],
            vec![1, 3, 8, 2, 6, 9, 2, 7, 4, 8, 10, 8, 10, 5, 1, 3, 10, 10, 2, 9],
            vec![2, 4, 1, 9, 2, 9, 7, 8, 2, 1, 4, 10, 5, 2, 7, 6, 5, 7, 2, 6],
            vec![4, 5, 1, 4, 2, 3, 3, 4, 1, 8, 8, 2, 6, 9, 5, 9, 6, 3, 9, 3],
            vec![3, 1, 1, 8, 6, 8, 8, 7, 9, 3, 2, 1, 8, 2, 4, 7, 3, 1, 2, 4],
            vec![5, 9, 8, 6, 10, 4, 10, 3, 4, 10, 10, 10, 1, 7, 8, 8, 7, 7, 8, 8],
            vec![1, 4, 6, 1, 6, 1, 2, 10, 5, 10, 2, 6, 2, 4, 5, 5, 3, 5, 1, 5],
            vec![5, 6, 9, 10, 6, 6, 10, 6, 4, 1, 5, 3, 9, 5, 2, 10, 9, 9, 5, 1],
            vec![10, 9, 4, 6, 9, 5, 3, 7, 10, 1, 6, 8, 1, 1, 10, 9, 5, 7, 7, 5],
            vec![2, 6, 6, 6, 6, 2, 9, 4, 7, 5, 3, 2, 10, 3, 4, 5, 10, 9, 1, 7],
            vec![5, 2, 4, 9, 8, 4, 8, 2, 4, 1, 3, 7, 6, 8, 1, 6, 8, 8, 10, 10],
            vec![9, 6, 3, 1, 8, 5, 7, 8, 7, 2, 1, 8, 2, 8, 3, 7, 4, 8, 7, 7],
            vec![8, 4, 4, 9, 7, 10, 6, 2, 1, 5, 8, 5, 1, 1, 1, 9, 1, 3, 5, 3],
        ];
        assert_eq!(hungarian(matrix.clone()).iter()
            .enumerate()
            .map(|(i, &j)| matrix[i][j])
            .sum::<u64>(), 22);
    }


    // From https://stackoverflow.com/questions/17419595/hungarian-kuhn-munkres-algorithm-oddity
    #[test]
    fn test_stack_overflow_a() {
        let matrix = vec![ 
            vec![0, 7, 0, 0, 0],
            vec![0, 8, 0, 0, 6], 
            vec![5, 0, 7, 3, 4], 
            vec![5, 0, 5, 9, 3], 
            vec![0, 4, 0, 0, 9],
        ];
        assert_eq!(hungarian(matrix.clone()).iter()
            .enumerate()
            .map(|(i, &j)| matrix[i][j])
            .sum::<u64>(), 3);
        assert_eq!(hungarian(matrix), vec![4, 2, 3, 1, 0]);
    }

    // From https://stackoverflow.com/questions/26893961/cannot-solve-hungarian-algorithm
    //TODO: figure out why this works for padded square input but not rectangle
    #[test]
    fn test_stack_overflow_b() {
        let matrix = vec![
           vec![  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0],
           vec![ 53,207,256,207,231,348,348,348,231,244,244, 0, 0, 0],
           vec![240, 33, 67, 33, 56,133,133,133, 56, 33, 33, 0, 0, 0],
           vec![460,107,200,107,122,324,324,324,122, 33, 33, 0, 0, 0],
           vec![167,340,396,340,422,567,567,567,422,442,442, 0, 0, 0],
           vec![167,367,307,367,433,336,336,336,433,158,158, 0, 0, 0],
           vec![160, 20, 37, 20, 31, 70, 70, 70, 31, 22, 22, 0, 0, 0],
           vec![200,307,393,307,222,364,364,364,222,286,286, 0, 0, 0],
           vec![33 ,153,152,153,228,252,252,252,228, 78, 78, 0, 0, 0],
           vec![93 ,140,185,140, 58,118,118,118, 58, 44, 44, 0, 0, 0],
           vec![0  ,  7, 22,  7, 19, 58, 58, 58, 19,  0,  0, 0, 0, 0],
           vec![67 ,153,241,153,128,297,297,297,128, 39, 39, 0, 0, 0],
           vec![73 ,253,389,253,253,539,539,539,253, 36, 36, 0, 0, 0],
           vec![173,267,270,267,322,352,352,352,322,231,231, 0, 0, 0],
        ];
        assert_eq!(hungarian(matrix.clone()).iter()
            .enumerate()
            .map(|(i, &j)| matrix[i][j])
            .sum::<u64>(), 828);
    }

    // From https://stackoverflow.com/questions/46803600/hungarian-algorithm-wikipedia-method-doesnt-work-for-this-example
    #[test]
    fn test_stack_overflow_c() {
        let matrix = vec![
            vec![35, 0, 0, 0],
            vec![0 ,30, 0, 5],
            vec![55, 5, 0,10],
            vec![0 ,45,30,45],
        ];
        assert_eq!(hungarian(matrix.clone()).iter()
            .enumerate()
            .map(|(i, &j)| matrix[i][j])
            .sum::<u64>(), 5);
        assert_eq!(hungarian(matrix), vec![1, 3, 2, 0]);
    }

    // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    #[test]
    fn test_stack_overflow_d() {
        let matrix = vec![
            vec![2,1,0,0,0,3],
            vec![2,0,4,5,2,7],
            vec![0,7,0,0,0,5],
            vec![3,2,3,1,2,0],
            vec![0,0,6,3,3,5],
            vec![3,4,5,2,0,3],
        ];
        assert_eq!(hungarian(matrix), vec![3, 1, 2, 5, 0, 4]);
    }

    // From https://stackoverflow.com/questions/37687045/hungarian-algorithm-dead-end
    #[test]
    fn test_stack_overflow_e() {
        let matrix = vec![
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 0],
            vec![0, 0, 1, 2],
            vec![0, 0, 3, 4],
        ];
        assert_eq!(hungarian(matrix), vec![3, 2, 1, 0]);
    }

    // From https://www.wikihow.com/Use-the-Hungarian-Algorithm
    #[test]
    fn test_wikihow_example() {
        let matrix = vec![
            vec![10, 19, 8, 15, 19],
            vec![10, 18, 7, 17, 19],
            vec![13, 16, 9, 14, 19],
            vec![12, 19, 8, 18, 19],
            vec![14, 17, 10, 19, 19]
        ];
        assert_eq!(hungarian(matrix), vec![0, 2, 3, 4, 1]);
    }

    #[test]
    fn test_worst_case() {
        for max in 1..100 {
            let mut matrix = vec![vec![0; max]; max];
            let mut n: u64 = 0;
            
            for i in 0..max {
                for j in 0..max {
                    matrix[i][j] = n;
                    n += 1; 
                }
            }
            assert_eq!(hungarian(matrix), (0..max).rev().collect::<Vec<_>>());
        }
    }

    #[test]
    fn test_stress() {
        for max in 1..100 {
            let mut matrix = vec![vec![0; max]; max];
            
            for i in 0..max {
                for j in 0..max {
                    matrix[i][j] = (i*j) as u64;
                }
            }
            assert_eq!(hungarian(matrix), (0..max).rev().collect::<Vec<_>>());
        }
    }

    #[test]
    fn test_large() {
        let max = 250;
        let mut matrix = vec![vec![0; max]; max];
        
        for i in 0..max {
            for j in 0..max {
                matrix[i][j] = (i*j) as u64;
            }
        }
        assert_eq!(hungarian(matrix), (0..max).rev().collect::<Vec<_>>());
    }
}
