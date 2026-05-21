#[derive(Debug, Clone, Copy)]
pub struct Cell {
    test: f32,
}

impl Default for Cell {
    fn default() -> Self {
        Cell { test: 1.0 }
    }
}
