use crate::mesh;
use crate::state::Euler1DState;


#[derive(Debug)]
pub struct SolutionField {
    pub sol: Vec<Euler1DState>,
}


impl Euler1DState {
    pub fn zero() -> Euler1DState {
        Euler1DState { rho: 0.0, rho_v: 0.0, energy: 0.0 }
    }
}


impl SolutionField {
    pub fn new(mesh: &mesh::Mesh1D) -> SolutionField {
        SolutionField {
            sol: vec![Euler1DState::zero(); mesh.elements.len()],
        }
    }
}


pub fn initialize_shocktube(mesh: &mesh::Mesh1D, field: &mut SolutionField, left_state: Euler1DState, right_state: Euler1DState) {
    if mesh.elements.len() != field.sol.len() {
        eprintln!("Invalid field shape! Did you forget to initialize the field on your mesh?");
        return;
    }

    let threshold = (mesh.extents.1 - mesh.extents.0) * 0.5;

    for i in 0..mesh.elements.len() {
        field.sol[i] = match mesh.element_center(i) < threshold {
            true => left_state,
            false => right_state
        }
    }
}