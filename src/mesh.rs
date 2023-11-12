use std::fmt;
use crate::boundary::Boundary;


pub struct Node1D {
    x: f64,
    boundary: Boundary,
}

#[derive(Debug)]
pub struct Elem1D {
    nodes: [usize; 2],
}


#[derive(Debug)]
pub struct Mesh1D {
    pub elements: Vec<Elem1D>,
    pub nodes: Vec<Node1D>,
    pub extents: (f64, f64),
}


impl fmt::Display for Node1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.boundary {
            Boundary::Connected(conn) => write!(f, "{:.2}/Connected({})", self.x, conn),
            Boundary::Wall => write!(f, "{:.2}/Wall", self.x),
        }
    }
}
impl fmt::Debug for Node1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self)
    }
}


impl fmt::Display for Elem1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{},{}]", self.nodes[0], self.nodes[1])
    }
}


impl Mesh1D {
    pub fn new() -> Mesh1D {
        Mesh1D { 
            elements: Vec::new(),
            nodes: Vec::new(),
            extents: (0.0, 0.0),
        }
    }
    
    pub fn create_equidistant(bounds: (f64, f64), nr_cells: u32) -> Mesh1D {
        let mut mesh = Mesh1D::new();
        let step_size = (bounds.1 - bounds.0) / nr_cells as f64;
        for i in 0..nr_cells {
            // Create node boundaries
            let boundary_a = match i {
                0 => Boundary::Wall,
                _ => Boundary::Connected(mesh.nodes.len() - 1), // Last of previous element
            };
            let boundary_b = match i {
                _ if i == nr_cells - 1 => Boundary::Wall,
                _ => Boundary::Connected(mesh.nodes.len() + 2), // First of next element
            };

            // Create the nodes for the element
            mesh.nodes.push(
                Node1D{x: bounds.0 + ( i      as f64)*step_size, boundary: boundary_a}
            );
            mesh.nodes.push(
                Node1D{x: bounds.0 + ((i + 1) as f64)*step_size, boundary: boundary_b}
            );

            // Create and insert a new element
            mesh.elements.push(
                Elem1D {
                    nodes: [mesh.nodes.len() - 2, mesh.nodes.len() -1]
                }
            )
        }

        // Use the extents that were provided to the constructor
        mesh.extents = bounds;
        mesh
    }


    pub fn element_center(&self, element: usize) -> f64 {
        let [n1, n2] = self.elements[element].nodes;
        (self.nodes[n1].x + self.nodes[n2].x) * 0.5
    }

    pub fn element_extents(&self, element: usize) -> (f64, f64) {
        let [n1, n2] = self.elements[element].nodes;
        (self.nodes[n1].x, self.nodes[n2].x)
    }
}

impl fmt::Display for Mesh1D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Mesh1D {{\n  elements: ", )?;
        let mut it = self.elements.iter().peekable();
        while let Some(el) = it.next() {
            // Print a comma if an element follows.
            match it.peek() {
                Some(_) => write!(f, "{}, ", el)?,
                None    => write!(f, "{}"  , el)?,
            }
        }

        write!(f, "\n  nodes: ", )?;
        let mut it = self.nodes.iter().peekable();
        while let Some(node) = it.next() {
            // Print a comma if an element follows.
            match it.peek() {
                Some(_) => write!(f, "{}, ", node)?,
                None    => write!(f, "{}"  , node)?,
            }
        }
        write!(f, "\n}}")
    }

}