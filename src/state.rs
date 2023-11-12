// Typically using conservative quantities rho, rho*v, rho*e
#[derive(Clone, Copy, Debug)]
pub struct Euler1DState {
    pub rho:   f64,
    pub rho_v: f64,
    pub energy: f64,
}

// The same as the Euler1DState, but expressed in different variables.
// These two representations only become linked in the context of a gas law.
#[derive(Clone, Copy, Debug)]
pub struct Euler1DPrimState {
    pub rho: f64,  // Density
    pub v:   f64,    // Velocity, sometimes also denoted as u
    pub p:   f64,    // Pressure
}


pub trait GasLawInterface {
    fn gas_energy(&self, state: Euler1DPrimState) -> f64;
    fn gas_pressure(&self, state: Euler1DState) -> f64;
    fn speed_of_sound(&self, state: Euler1DPrimState) -> f64;
    fn cons_to_prim(&self, conservative: Euler1DState) -> Euler1DPrimState;
    fn prim_to_cons(&self, primitive: Euler1DPrimState) -> Euler1DState;
}


pub struct IdealGasLaw {
    pub gamma: f64,
}


impl GasLawInterface for IdealGasLaw {
    fn gas_energy(&self, state: Euler1DPrimState) -> f64 {
        let Euler1DPrimState {rho, v, p} = state;
        // Internal energy
        let e = p / ((self.gamma - 1.0) * rho);
        rho * (0.5*v*v + e)
    }

    fn gas_pressure(&self, state: Euler1DState) -> f64 {
        let Euler1DState { rho, rho_v, energy } = state;
        let v = rho_v / rho;
        let e = (energy / rho) - 0.5 * v*v;
        e * (self.gamma - 1.0) * rho
    }

    fn speed_of_sound(&self, state: Euler1DPrimState) -> f64 {
        ((self.gamma * state.p) / state.rho).sqrt()
    }

    fn cons_to_prim(&self, conservative: Euler1DState) -> Euler1DPrimState {
        let Euler1DState { rho, rho_v, energy: _ } = conservative;
        Euler1DPrimState {
            rho: rho,
            v:   rho_v / rho,
            p:   self.gas_pressure(conservative),
        }
    }

    fn prim_to_cons(&self, primitive: Euler1DPrimState) -> Euler1DState {
        let Euler1DPrimState {rho, v, p: _} = primitive;
        Euler1DState {
            rho: rho,
            rho_v: rho * v,
            energy: self.gas_energy(primitive), 
        }
    }
}

pub const STANDARD_GAS_MODEL: IdealGasLaw = IdealGasLaw{gamma: 1.4};
