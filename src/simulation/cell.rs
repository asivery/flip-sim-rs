#[derive(Debug, Clone, Copy)]
enum CellType{Solid,Liquid,Gas}

pub struct Cell {
    pub type: CellType,
}

impl Default for Cell {
    fn default() -> Self {
        Cell { test: false }
    }
}
