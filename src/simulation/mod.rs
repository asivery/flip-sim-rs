mod cell;

use cell::Cell;

#[derive(Debug)]
pub struct Simulation {
    width: usize,
    height: usize,

    grid: Vec<Cell>,

}

impl Simulation {
    pub fn new(width: usize, height: usize) -> Simulation {
        let total_cells = width * height;

        Simulation {
            width,
            height,

            grid: vec![Cell::default(); total_cells]
        }
    }

    pub fn step(&mut self) {

    }
}
