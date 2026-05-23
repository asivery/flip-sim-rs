pub mod config;
pub mod cell;
pub mod particle;

use config::Config;
use cell::*;
use particle::*;

use rand::prelude::*;


#[derive(Debug)]
pub struct Simulation {
    pub config: Config,
    
    // double buffering, to avoid weird borrow checker workarounds
    pub read_grid: Vec<Cell>,
    pub write_grid: Vec<Cell>,

    pub f_num_x: usize,
    pub f_num_y: usize,
    pub f_num_cells: usize,
    pub h: f32,
    pub f_inv_spacing: f32,

    pub particles: Vec<Particle>,

    pub particle_rest_density: f32,
    pub p_inv_spacing: f32,
    pub p_num_x: usize,
    pub p_num_y: usize,
    pub p_num_cells: usize,

    pub num_cell_particles: Vec<i32>,
    pub first_cell_particle: Vec<usize>,
    pub cell_particle_ids: Vec<usize>,

    pub num_particles: usize



}

impl Simulation {
    pub fn new(config: &Config) -> Simulation {

        // let total_cells = config.width * config.height;
        let f_num_x = (config.width as f32 / config.spacing).floor() as usize + 1;
        let f_num_y = (config.height as f32 / config.spacing).floor() as usize + 1;
        let f_num_cells = f_num_x * f_num_y;
        let h = f32::max(config.width as f32 / f_num_x as f32, config.height as f32 / f_num_y as f32);
        let f_inv_spacing = 1.0 / h;

        let p_inv_spacing = 1.0 / (2.2 * config.particle_radius);
        let p_num_x = (config.width as f32 * p_inv_spacing).floor() as usize + 1;
        let p_num_y = (config.height as f32 * p_inv_spacing).floor() as usize + 1;
        let p_num_cells = p_num_x * p_num_y;

        let mut read_grid = vec![Cell::default(); f_num_cells];
        read_grid[0].cell_type = CellTypes::Solid; // TODO: remove after testing 

        let write_grid = read_grid.clone();

        Simulation {
            config: config.clone(),
            read_grid,
            write_grid,
            particles: vec![Particle::default(); config.max_particles],

            f_num_x,
            f_num_y,
            f_num_cells,
            h,
            f_inv_spacing,

            particle_rest_density: 0.0,
            p_inv_spacing,
            p_num_x,
            p_num_y,
            p_num_cells,

            num_cell_particles: vec![0;p_num_cells],
            first_cell_particle: vec![0;config.max_particles + 1],
            cell_particle_ids: vec![0;config.max_particles],
            num_particles: 0
        }
    }

    pub fn get_cell(&self, x: usize, y: usize) -> &Cell {
        &self.read_grid[y * self.f_num_x + x]
    }

    pub fn step(&mut self) {
        let mut rng = rand::rng();

        self.read_grid.shuffle(&mut rng);
    }

    pub fn integrate_particles(&mut self, dt:f32,gravity: (f32,f32)) {
        let (gx,gy) = gravity;
        for part in self.particles.iter_mut() {
            part.vx += dt * gx;
            part.vy += dt * gy;
            part.x += part.vx * dt;
            part.y += part.vy * dt;
        }
    }

    pub fn push_particles_apart(&mut self, num_iters: usize){
        self.num_cell_particles.fill(0);

        for i in 0..self.num_particles{
            let p = &self.particles[i];
            let (x,y) = (p.x, p.y);

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
            let p = &self.particles[i];
            let (x,y) = (p.x, p.y);

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
                let p = &self.particles[i];
                let (px,py) = (p.x, p.y);

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
                            if id == i{
                                continue;
                            }
                            let q = &self.particles[id];
                            let (qx,qy) = (q.x, q.y);

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

                            
                            self.particles[i].x -= dx;
                            self.particles[i].y -= dy;
                            self.particles[id].x += dx;
                            self.particles[id].y += dy;

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

        self.write_grid.copy_from_slice(&self.read_grid);
        self.write_grid.iter_mut().for_each(|cell| {cell.particle_density = 0.0;});

        for i in 0..self.num_particles{
            let p = &self.particles[i];
            let x = p.x.clamp(h, (self.f_num_x as f32 - 1.0) * h);
            let y = p.y.clamp(h, (self.f_num_y as f32 - 1.0) * h);

            let x0 = (((x - h2) * h1).floor() as usize).min(self.f_num_x - 2);
            let tx = ((x - h2) - x0 as f32 * h) * h1;
            let x1 = (x0 + 1).min(self.f_num_x - 2);

            let y0 = (((y - h2) * h1).floor() as usize).min(self.f_num_y - 2);
            let ty = ((y - h2) - y0 as f32 * h) * h1;
            let y1 = (y0 + 1).min(self.f_num_y - 2);

            let sx = 1.0 - tx;
            let sy = 1.0 - ty;

            self.write_grid[x0 * n + y0].particle_density += sx * sy;
            self.write_grid[x1 * n + y0].particle_density += tx * sy;
            self.write_grid[x1 * n + y1].particle_density += tx * ty;
            self.write_grid[x0 * n + y1].particle_density += sx * ty;
        }

        if self.particle_rest_density == 0.0{
            let mut sum = 0.0;
            let mut num_fluid_cells = 0;

            for i in 0..self.f_num_cells{
                if self.write_grid[i].cell_type == CellTypes::Liquid{
                    sum += self.write_grid[i].particle_density;
                    num_fluid_cells += 1;
                }
            }
            if num_fluid_cells > 0{
                self.particle_rest_density = sum / num_fluid_cells as f32;
            }
        }

        std::mem::swap(&mut self.read_grid, &mut self.write_grid);
    }

    // TODO: Transfer Velocities
    // TODO: solve incompressibility
    // TODO: update cell colors
    // TODO: Set Sci Color whatever that means
    // TODO: simulate
}
