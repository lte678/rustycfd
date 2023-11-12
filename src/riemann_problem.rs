// Author: Leon Teichroeb
// Date:   11.11.2023
// Sources:
// - Eleuterio, F. T. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics (3rd ed.). Springer

use std::error::Error;

use crate::state::{Euler1DPrimState, Euler1DState, IdealGasLaw, GasLawInterface};
use crate::solution_field::SolutionField;
use crate::mesh::Mesh1D;


const MAX_ITERATIONS: i32 = 100;
const CONVERGENCE_CRITERIUM: f64 = 1e-6;


fn calc_flux(p_star: f64, state: Euler1DPrimState, a: f64, gamma: f64) -> (f64, f64) {
    // This function is implemented according to [1] p.119
    // K represents either the L or R state in the riemann problem
    // Parameters:
    // - p_star: pressure in star region (approximated in interative solver)
    // - state: state in the left or right region
    // - a: speed of sound in left or right region
    // - gamma: Ideal gas constant
    // Returns 
    // - Tuple(flux, derivative of flux)

    let Euler1DPrimState{rho, v: _, p} = state;

    if p_star > p {
        // This is a shock
        let A_K = 2.0/((gamma + 1.0)*rho);
        let B_K = ((gamma - 1.0)/(gamma + 1.0)) * p;

        let f = (p_star - p) * (A_K / (p_star + B_K)).sqrt();
        let f_prime = (A_K / (B_K + p_star)).sqrt() * (1.0 - (p_star - p)/(2.0 * (B_K + p_star)));

        return (f, f_prime);
    } else { 
        // This is a rarefaction
        let exp = (gamma - 1.0) / (2.0 * gamma);
        
        let f = (2.0*a)/(gamma - 1.0) * ((p_star / p).powf(exp) - 1.0);
        let f_prime = 1.0 / (rho * a) * (p_star / p).powf(-0.5*(gamma + 1.0)/gamma);

        return (f, f_prime);
    }
}


fn solve_p_star(left_state: Euler1DPrimState, right_state: Euler1DPrimState, eq_of_state: &IdealGasLaw) -> Result<(f64, f64), Box<dyn Error>> {
    // Calculates the pressure and velocity within the star region of the riemann problem iteratively.

    let a_L = eq_of_state.speed_of_sound(left_state);
    let a_R = eq_of_state.speed_of_sound(right_state);
    let gamma = eq_of_state.gamma;

    let du_crit = (2.0 * a_L) / (gamma - 1.0) + (2.0 * a_R) / (gamma - 1.0);
    let du = right_state.v - left_state.v;
    if du > du_crit {
        return Err(From::from("Pressure positivity condition is not satisfied."));
    }

    let mut p_star = 0.5 * (left_state.p + right_state.p);

    for i in 1..(MAX_ITERATIONS+1) {
        // Calculate and combine fluxes
        let (f_L, f_L_prime) = calc_flux(p_star, left_state, a_L, gamma);
        let (f_R, f_R_prime) = calc_flux(p_star, right_state, a_R, gamma);
        let f = f_L + f_R + du;
        let f_prime = f_L_prime + f_R_prime;

        // Perform Newton iteration
        let p_star_prev = p_star;
        p_star = p_star_prev - f / f_prime;

        // Check if the solution blew up :)
        if p_star.is_nan() || p_star.is_infinite() {
            return Err(From::from("Solution out of bounds."));
        }

        // Check convergence criteria
        let relative_change = ((p_star - p_star_prev) / (0.5 * (p_star + p_star_prev))).abs();

        //println!("Relative change in iteration {i}: {relative_change}");
        if relative_change < CONVERGENCE_CRITERIUM {
            // Also calculate the velocity
            let (f_L, _) = calc_flux(p_star, left_state, a_L, gamma);
            let (f_R, _) = calc_flux(p_star, right_state, a_R, gamma);
            let v_star = 0.5*(left_state.v + right_state.v) + 0.5*(f_R - f_L);

            println!("Solved p*={p_star:.4e} Pa, u*={v_star:.4e} m/s in {i} iterations");
            return Ok((p_star, v_star));
        }

        // Makes sure the derivative does not take us into the negative
        if p_star < 0.0 {
            p_star = 1e-6;
        }
    }

    Err(From::from("Solution did not converge."))
}


fn rarefaction_state(W: Euler1DPrimState, a: f64, gamma: f64, slope: f64) -> Euler1DPrimState {
    // For a left-travelling rarefaction, a should be positive, for a right travelling rarefaction, a should be negative.
    // This has little physical meaning, but satisfies the equations

    let k1 = (gamma-1.0)/((gamma + 1.0)*a);
    let k2 = 2.0/(gamma + 1.0);
    let k3 = (k2 + k1*(W.v - slope)).powf(2.0/(gamma - 1.0));
    let rho = W.rho * k3;
    let v = k2 * (a + (gamma - 1.0)/2.0 * W.v + slope);
    let p = W.p * k3.powf(gamma);

    Euler1DPrimState { rho, v, p }
}


pub fn solve_riemann(mesh: &Mesh1D, field: &mut SolutionField, left_state: Euler1DState, right_state: Euler1DState, eq_of_state: &IdealGasLaw, time: f64) -> Result<(), Box<dyn Error>> {
    if mesh.elements.len() != field.sol.len() {
        eprintln!("Invalid field shape! Did you forget to initialize the field on your mesh?");
        return Err(From::from("Invalid field shape."));
    }
    
    // For ideal gases, solves the shock tube problem for most left and right states
    let W_L = eq_of_state.cons_to_prim(left_state);
    let W_R = eq_of_state.cons_to_prim(right_state);
    let a_L = eq_of_state.speed_of_sound(W_L);
    let a_R = eq_of_state.speed_of_sound(W_R);
    let gamma = eq_of_state.gamma;
    // Iteratively solve for p_star
    let (p_star, v_star) = solve_p_star(W_L, W_R, eq_of_state)?;

    // Find the propagation speeds and rho_star for the left and right discontinuities
    let S_left_head: f64;
    let S_left_tail: f64;
    let star_L: Euler1DPrimState;
    if p_star > W_L.p {
        // Left shock
        S_left_head = W_L.v - a_L * ((gamma + 1.0)/(2.0*gamma) * p_star/W_L.p + (gamma - 1.0)/(2.0*gamma)).sqrt();
        S_left_tail = S_left_head;
        let k1 = p_star / W_L.p + (gamma - 1.0)/(gamma + 1.0);
        let k2 = ((gamma - 1.0)/(gamma + 1.0)) * p_star / W_L.p + 1.0;
        let rho_star_L = W_L.rho * (k1 / k2);
        star_L = Euler1DPrimState{rho: rho_star_L, v: v_star, p: p_star};
    } else {
        // Left rarefaction
        // Speed of sound after rarefaction
        let a_star_L = a_L * (p_star/W_L.p).powf((gamma - 1.0)/(2.0*gamma));
        S_left_head = W_L.v - a_L;
        S_left_tail = v_star - a_star_L;

        let rho_star_L = W_L.rho * (p_star / W_L.p).powf(1.0/gamma);
        star_L = Euler1DPrimState{rho: rho_star_L, v: v_star, p: p_star};
    }

    let S_right_head: f64;
    let S_right_tail: f64;
    let star_R: Euler1DPrimState;
    if p_star > W_R.p {
        // Right shock
        S_right_head = W_R.v + a_R * ((gamma + 1.0)/(2.0*gamma) * p_star/W_R.p + (gamma - 1.0)/(2.0*gamma)).sqrt();
        S_right_tail = S_right_head;
        let k1 = p_star / W_R.p + (gamma - 1.0)/(gamma + 1.0);
        let k2 = ((gamma - 1.0)/(gamma + 1.0)) * p_star / W_R.p + 1.0;
        let rho_star_R = W_R.rho * (k1 / k2);
        star_R = Euler1DPrimState{rho: rho_star_R, v: v_star, p: p_star};
    } else {
        // Right rarefaction
        // Speed of sound after rarefaction
        let a_star_R = a_R * (p_star/W_R.p).powf((gamma - 1.0)/(2.0*gamma));
        S_right_head = W_R.v + a_R;
        S_right_tail = v_star + a_star_R;

        let rho_star_R = W_R.rho * (p_star / W_R.p).powf(1.0/gamma);
        star_R = Euler1DPrimState{rho: rho_star_R, v: v_star, p: p_star};
    }

    println!("Left propagation speeds: {S_left_head} m/s, {S_left_tail} m/s");
    println!("Right propagation speeds: {S_right_head} m/s, {S_right_tail} m/s");

    // The center of the riemann problem
    let center = (mesh.extents.1 - mesh.extents.0) * 0.5;
    // Calculate the transition points
    let trans_left_head = S_left_head * time;
    let trans_left_tail = S_left_tail * time;
    let trans_discontinuity = v_star * time;
    let trans_right_head = S_right_head * time;
    let trans_right_tail = S_right_tail * time;

    // Build the field solution
    for i in 0..mesh.elements.len() {
        let x = mesh.element_center(i) - center;
        
        if x < trans_left_head {
            field.sol[i] = left_state;
        } else if x < trans_left_tail {
            field.sol[i] = eq_of_state.prim_to_cons(rarefaction_state(W_L, a_L, gamma, x / time));
        } else if x < trans_discontinuity {
            field.sol[i] = eq_of_state.prim_to_cons(star_L);
        } else if x < trans_right_tail {
            field.sol[i] = eq_of_state.prim_to_cons(star_R);
        } else if x < trans_right_head {
            field.sol[i] = eq_of_state.prim_to_cons(rarefaction_state(W_R, -a_R, gamma, x / time));
        } else {
            field.sol[i] = right_state;
        }
        
    }

    Ok(())
}