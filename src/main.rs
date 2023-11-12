#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

mod boundary;
mod mesh;
mod solution_field;
mod state;
mod export;
mod riemann_problem;

use riemann_problem::solve_riemann;
use solution_field::SolutionField;
use state::{Euler1DPrimState, STANDARD_GAS_MODEL, GasLawInterface};
use export::csv_export;

fn main() {
    let cfd_mesh = mesh::Mesh1D::create_equidistant((0.0, 1.0), 100);
    let mut cfd_solution = SolutionField::new(&cfd_mesh);

    // Initial shoch tube state
    // This is the standard SOD test problem
    let left_state = STANDARD_GAS_MODEL.prim_to_cons(Euler1DPrimState {
        rho: 5.99924,
        v: 19.5975,
        p: 460.894
    });
    let right_state = STANDARD_GAS_MODEL.prim_to_cons(Euler1DPrimState {
        rho: 5.99242,
        v: -6.19633,
        p: 46.0950
    });
    let t = 0.035;

    solution_field::initialize_shocktube(&cfd_mesh, &mut cfd_solution, left_state, right_state);
    let export_fname = "shocktube.csv";
    match csv_export::export_solution(export_fname, &cfd_mesh, &cfd_solution, &STANDARD_GAS_MODEL) {
        Err(msg) => eprintln!("Failed to export with error {msg}"),
        Ok(()) => println!("Exported solution to {export_fname}"),
    }


    let res = solve_riemann(&cfd_mesh, &mut cfd_solution, left_state, right_state, &STANDARD_GAS_MODEL, t);
    match res {
        Err(msg) => eprintln!("Failed to solve riemann problem: {msg}"),
        Ok(()) => (),
    }
    let export_fname = "solution.csv";
    match csv_export::export_solution(export_fname, &cfd_mesh, &cfd_solution, &STANDARD_GAS_MODEL) {
        Err(msg) => eprintln!("Failed to export with error {msg}"),
        Ok(()) => println!("Exported solution to {export_fname}"),
    }
}
