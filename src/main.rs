extern crate hungarian;
extern crate cpuprofiler;

use hungarian::minimize;
use cpuprofiler::PROFILER;

pub fn main() {
    let max = 1000;
    let mut matrix = vec![0; max * max];
    
    for i in 0..max {
        for j in 0..max {
            matrix[max*i + j] = ((i + 1)*(j + 1)) as i32;
        }
    }

    // Unlock the mutex and start the profiler
    PROFILER.lock().unwrap().start("./my-prof.profile").expect("Couldn't start");
    println!("{:?}", minimize(&matrix, max, max));
    PROFILER.lock().unwrap().stop().expect("Couldn't stop");
}
