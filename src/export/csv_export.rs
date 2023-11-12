use std::error::Error;

use crate::mesh::Mesh1D;
use crate::state::{Euler1DPrimState, GasLawInterface};
use crate::solution_field::{SolutionField};

use csv;


pub fn export_solution(filename: &str, mesh: &Mesh1D, field: &SolutionField, gas_law: &dyn GasLawInterface) -> Result<(), Box<dyn Error>> {
    let mut writer = csv::Writer::from_path(filename)?;
    if mesh.elements.len() != field.sol.len() {
        eprintln!("Invalid field shape! Did you forget to initialize the field on your mesh?");
        return Err(From::from("Invalid field shape"));
    }

    writer.write_record(&["elemIdx", "elemLeft", "elemRight", "elemCenter", "rho", "v", "p"])?;
    // Write to csv
    for i in 0..mesh.elements.len() {
        let (x0, x1) = mesh.element_extents(i);
        let x_center = mesh.element_center(i);
        let Euler1DPrimState { rho, v, p } = gas_law.cons_to_prim(field.sol[i]);
        writer.serialize((i, x0, x1, x_center, rho, v, p))?;
    }

    writer.flush()?;
    Ok(())
}