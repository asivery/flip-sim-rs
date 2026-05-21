use flip_sim_rs::simulation::*;
use flip_sim_rs::front_cli::*;

fn main() {
    let sim = Simulation::new(8, 8); 
    let mut front = FrontCLI::new(sim);

    front.run()
}
