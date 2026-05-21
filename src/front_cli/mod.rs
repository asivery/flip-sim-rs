use crate::simulation::Simulation; 

use std::{thread, time::Duration};

#[derive(Debug)]
pub struct FrontCLI {
    sim: Simulation,

    frames_per_second: f32
}

impl FrontCLI {
    pub fn new(sim: Simulation) -> Self {
        Self {
            sim,
            frames_per_second: 2.0
        }
    }

    pub fn run(&mut self) {
        loop {
            println!{"Frame"};
            self.sim.step();
            thread::sleep(Duration::from_secs_f32(1.0 / self.frames_per_second));
        }
    }
}
