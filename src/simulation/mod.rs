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
    pub h: f32,
    pub f_inv_spacing: f32,

    pub u: Vec<f32>,
    pub v: Vec<f32>,
    pub du: Vec<f32>,
    pub dv: Vec<f32>,
    pub prev_u: Vec<f32>,
    pub prev_v: Vec<f32>,
    pub p: Vec<f32>,
    pub s: Vec<f32>,
    pub cell_type: Vec<CellTypes>,
    pub cell_color: Vec<f32>,

    pub particle_pos: Vec<f32>,
    pub particle_vel: Vec<f32>,
    pub particle_density: Vec<f32>,
    pub particle_rest_density: f32,

    pub p_inv_spacing: f32,
    pub p_num_x: usize,
    pub p_num_y: usize,
    pub p_num_cells: usize,

    pub num_cell_particles: Vec<i32>,
    pub first_cell_particle: Vec<usize>,
    pub cell_particle_ids: Vec<i32>,

    pub num_particles: i32



}

impl Simulation {
    pub fn new(config: Config) -> Simulation {

        // let total_cells = config.width * config.height;
        let f_num_x = (config.width as f32 / config.spacing).floor() as usize + 1;
        let f_num_y = (config.height as f32 / config.spacing).floor() as usize + 1;
        let f_num_cells = f_num_x * f_num_y;
        let h = f32::max(config.width as f32 / f_num_x as f32, config.height as f32 / f_num_y as f32);
        let f_inv_spacing = 1.0 / h;

        let u = vec![0.0;f_num_cells];
        let v = vec![0.0;f_num_cells];
        let du = vec![0.0;f_num_cells];
        let dv = vec![0.0;f_num_cells];
        let prev_u = vec![0.0;f_num_cells];
        let prev_v = vec![0.0;f_num_cells];
        let p = vec![0.0;f_num_cells];
        let s = vec![0.0;f_num_cells];

        let particle_pos = vec![0.0; 2 * config.max_particles];
        let particle_vel = vec![0.0; 2 * config.max_particles];
        let particle_density = vec![0.0; f_num_cells];

        let p_inv_spacing = 1.0 / (2.2 * config.particle_radius);
        let p_num_x = (config.width as f32 * p_inv_spacing).floor() as usize + 1;
        let p_num_y = (config.height as f32 * p_inv_spacing).floor() as usize + 1;
        let p_num_cells = p_num_x * p_num_y;
        let first_cell_particle = vec![0;config.max_particles + 1];
        let cell_particle_ids = vec![0;config.max_particles as usize];

        let mut sim = Simulation {
            config,
            grid: vec![Cell::default(); f_num_cells],
            f_num_x,
            f_num_y,
            f_num_cells,
            h,
            f_inv_spacing,

            u,
            v,
            du,
            dv,
            prev_u,
            prev_v,
            p,
            s,
            cell_color: vec![0.0; 3*f_num_cells],
            cell_type: vec![CellTypes::Liquid; f_num_cells],

            particle_pos,
            particle_vel,
            particle_density,
            particle_rest_density: 0.0,
            p_inv_spacing,
            p_num_x,
            p_num_y,
            p_num_cells,
            num_cell_particles: vec![0;p_num_cells],
            first_cell_particle,
            cell_particle_ids,
            num_particles: 0
        };

        sim.grid[0] = Cell { cell_type: CellTypes::Solid, color: (0,0,0) };

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

    pub fn push_particles_apart(&mut self, num_iters: usize){
        self.num_cell_particles.fill(0);

        for i in 0..self.num_particles{
            let x = self.particle_pos[2 *i as usize];
            let y = self.particle_pos[2*i as usize+1];

            let xi = ((x*self.p_inv_spacing).floor() as i32).clamp(0, (self.p_num_x - 1) as i32) as usize;
            let yi = ((y*self.p_inv_spacing).floor() as i32).clamp(0, (self.p_num_y - 1) as i32) as usize;
            let cell_nr = xi * self.p_num_y + yi;
            self.num_cell_particles[cell_nr] += 1;
        }

        let mut first = 0;
        for i in 0..self.p_num_cells{
            first += self.num_cell_particles[i];
            self.first_cell_particle[i] = first as usize;
        }
        self.first_cell_particle[self.p_num_cells] = first as usize;

        for i in 0..self.num_particles {
            let x = self.particle_pos[2 * i as usize];
            let y = self.particle_pos[2 * i as usize + 1];

            let xi = ((x * self.p_inv_spacing).floor() as i32)
                .clamp(0, (self.p_num_x - 1) as i32) as usize;
            let yi = ((y * self.p_inv_spacing).floor() as i32)
                .clamp(0, (self.p_num_y - 1) as i32) as usize;
            let cell_nr = xi * self.p_num_y + yi;
            self.first_cell_particle[cell_nr] -= 1;
            let idx = self.first_cell_particle[cell_nr];
            self.cell_particle_ids[idx] = i;
        }

        let min_dist = 2.0 * self.config.particle_radius;
        let min_dist2 = min_dist * min_dist;

        for _iter in 0..num_iters{
            for i in 0..self.num_particles{
                let px = self.particle_pos[2*i as usize];
                let py = self.particle_pos[2*i as usize +1];

                let pxi = (px*self.p_inv_spacing).floor() as i32;
                let pyi = (py*self.p_inv_spacing).floor() as i32;
                let x0 = (pxi -1 ).max(0);
                let y0 = (pyi -1 ).max(0);
                let x1 = (pxi+1).min(self.p_num_x as i32 -1) as i32;
                let y1 = (pyi+1).min(self.p_num_y as i32 -1) as i32;

                for xi in x0..x1{
                    for yi in y0..y1{
                        let cell_nr = (xi * self.p_num_y as i32 + yi) as usize;
                        let first = self.first_cell_particle[cell_nr];
                        let last = self.first_cell_particle[cell_nr+1];
                        for j in first..last{
                            let id = self.cell_particle_ids[j];
                            if id == 1{
                                continue;
                            }
                            let qx = self.particle_pos[2*id as usize];
                            let qy = self.particle_pos[2*id as usize + 1];

                            let mut dx = qx - px;
                            let mut dy = qy - py;
                            let d2 = dx * dx + dy * dy;
                            if d2 > min_dist2 || d2 == 0.0{
                                continue;
                            }
                            let d = d2.sqrt();
                            let s = 0.5 * (min_dist - d) /d;
                            dx = dx * s;
                            dy = dy * s;
                            self.particle_pos[2*i as usize] -= dx;
                            self.particle_pos[2*i as usize + 1] -= dy;
                            self.particle_pos[2*id as usize] += dx;
                            self.particle_pos[2*id as usize + 1] += dx;

                            //here should be some color shit, and i cba doing that. as of writing this i have 0 clue if cell color is calculated by particle color, so if so we're fucked (maybe).
                        }
                    }
                }
            }
        }
    }

    pub fn update_particle_density(&mut self){
        let n = self.f_num_y;
        let h = self.h;
        let h1 = self.f_inv_spacing;
        let h2 = 0.5 * h;

        self.particle_density.fill(0.0);

        for i in 0..self.num_particles{
            let mut x = self.particle_pos[2 * i as usize];
            let mut y = self.particle_pos[2*i as usize + 1];

            x = x.clamp(h, (self.f_num_x - 1) as f32 * h);
            y = y.clamp(h, (self.f_num_y - 1) as f32 * h);

            let x0 = ((x - h2) * h1).floor() as i32;
            let tx = ((x - h2) - x0 as f32 * h) * h1;
            let x1 = (x0 + 1).min(self.f_num_y as i32-2);

            let y0 = ((y - h2) * h1).floor() as i32;
            let ty = ((y - h2) - (y0 as f32) * h) * h1;
            let y1 = (y0 + 1).min((self.f_num_y - 2) as i32);

            let sx = 1.0 - tx;
            let sy = 1.0 - ty;

            if x0 < self.f_num_x as i32 && y0 < self.f_num_y as i32 {self.particle_density[x0 as usize * n + y0 as usize] += sx*sy}
            if x1 < self.f_num_x as i32 && y0 < self.f_num_y as i32 {self.particle_density[x1 as usize * n + y0 as usize] += tx*sy}
            if x1 < self.f_num_x as i32 && y1 < self.f_num_y as i32 {self.particle_density[x1 as usize * n + y1 as usize] += tx*ty}
            if x0 < self.f_num_x as i32 && y1 < self.f_num_y as i32 {self.particle_density[x0 as usize * n + y1 as usize] += sx*ty}
        }

        if self.particle_rest_density == 0.0{
            let mut sum = 0.0;
            let mut num_fluid_cells = 0;

            for i in 0..self.f_num_cells{
                if self.cell_type[i] == CellTypes::Liquid{
                    sum += self.particle_density[i];
                    num_fluid_cells += 1;
                }
            }
            if num_fluid_cells > 0{
                self.particle_rest_density = sum / num_fluid_cells as f32;
            }
        }

// 			for (var xi = 1; xi < this.fNumX; xi++) {
// 				for (var yi = 1; yi < this.fNumY; yi++) {
// 					var cellNr = xi * n + yi;
// 					if (this.cellType[cellNr] != FLUID_CELL)
// 						continue;
// 					var hx = this.h;
// 					var hy = this.h;

// 					if (this.cellType[(xi - 1) * n + yi] == SOLID_CELL || this.cellType[(xi + 1) * n + yi] == SOLID_CELL)
// 						hx -= this.particleRadius;
// 					if (this.cellType[xi * n + yi - 1] == SOLID_CELL || this.cellType[xi * n + yi + 1] == SOLID_CELL)
// 						hy -= this.particleRadius;

// 					var scale = this.h * this.h / (hx * hy)
// 					d[cellNr] *= scale;
// 				}
// 			}
    }

    // TODO: Handle Particle Collisions, this might not be needed, idk really :3
    // TODO: Transfer Velocities
    // TODO: solve incompressibility
    // TODO: update cell colors
    // TODO: Set Sci Color whatever that means
    // TODO: simulate
}
