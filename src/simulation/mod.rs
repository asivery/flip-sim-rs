pub mod config;
pub mod cell;

use config::Config;
use cell::*;

use rand::prelude::*;

#[derive(Debug)]
pub struct Simulation {
    pub config: Config,
    pub grid: Vec<Cell>,
    pub f_num_x: usize,
    pub f_num_y: usize,
    pub f_num_cells: usize,

    pub particle_pos: Vec<f32>,
    pub particle_vel: Vec<f32>,
    pub particle_density: Vec<f32>

}

impl Simulation {
    pub fn new(config: Config) -> Simulation {
        
        let total_cells = config.width * config.height;
        let f_num_x = (config.width as f32 / config.spacing).floor() as usize + 1;
        let f_num_y = (config.height as f32 / config.spacing).floor() as usize + 1;
        let f_num_cells = f_num_x * f_num_y;
        let particle_pos = vec![0.0; 2 * config.max_particles as usize];
        let particle_vel = vec![0.0; 2 * config.max_particles as usize];
        let particle_density = vec![0.0; f_num_cells];

        let mut sim = Simulation {
            config,
            grid: vec![Cell::default(); total_cells],
            f_num_x,
            f_num_y,
            f_num_cells,
            particle_pos,
            particle_vel,
            particle_density
        };

        sim.grid[0] = Cell { cell_type: CellTypes::Solid };

        sim
    }

    pub fn get_cell(&self, x: usize, y: usize) -> &Cell {
        &self.grid[y * self.config.height + x]
    }

    pub fn step(&mut self) {
        let mut rng = rand::rng();

        self.grid.shuffle(&mut rng);
    }

    pub fn integrate_particles(&mut self, dt:f32,gravity: (f32,f32)) {
        let (gx,gy) = gravity;
        for i in 0..self.config.max_particles{

            let idx_x = 2 * i as usize;
            let idx_y = (2 * i + 1 )as usize;

            self.particle_vel[idx_x] += dt * gx;
            self.particle_vel[idx_y] += dt * gy;

            self.particle_pos[idx_x] += self.particle_vel[idx_x] * dt;
            self.particle_pos[idx_y] += self.particle_vel[idx_y] * dt;
            
        }

    }

    // TODO: Push Particles Apart
    // TODO: Handle Particle Collisions
    // TODO: Transfer Velocities
    // TODO: Update particle density
    // TODO: solve incompressibility
    // TODO: update cell colors
    // TODO: Set Sci Color whatever that means
    // TODO: simulate
}
