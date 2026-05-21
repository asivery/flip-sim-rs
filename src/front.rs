use crate::simulation::Simulation; 

pub trait FrontEnd {
    fn new(sim: Simulation) -> Self;
    fn init(&mut self);
    fn run(&mut self);
}
