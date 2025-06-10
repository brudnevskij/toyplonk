use ark_ff::Field;
use crate::fft::{fft, inverse_fft};
use crate::gate::Gate;
use crate::witness::Witness;

pub struct Circuit<F: Field>{
    pub gates: Vec<Gate<F>>,
    pub witness: Witness<F>,
    pub public_inputs: Vec<F>,
    pub domain: Vec<F>,
}


pub struct SelectorPolynomials<F: Field> {
    pub q_l: Vec<F>,
    pub q_r: Vec<F>,
    pub q_m: Vec<F>,
    pub q_o: Vec<F>,
    pub q_c: Vec<F>,
}

pub struct WitnessPolynomials<F: Field> {
    pub a: Vec<F>,
    pub b: Vec<F>,
    pub c: Vec<F>,
}

impl <F:Field>Circuit<F> {
    pub fn new(gates: Vec<Gate<F>>, witness: Witness<F>,public_inputs: Vec<F>,domain: Vec<F>) -> Circuit<F> {
        Circuit{
            gates,
            witness,
            public_inputs,
            domain,
        }
    }

    /// Extract and interpolate selector polynomials from gates
    pub fn get_selector_polynomials(&self) -> SelectorPolynomials<F> {
        let omega = self.domain[1]; // assumes domain = [1, ω, ω², ..., ωⁿ⁻¹]

        let mut q_l = Vec::new();
        let mut q_r = Vec::new();
        let mut q_m = Vec::new();
        let mut q_o = Vec::new();
        let mut q_c = Vec::new();

        for gate in &self.gates {
            q_l.push(gate.q_l);
            q_r.push(gate.q_r);
            q_m.push(gate.q_m);
            q_o.push(gate.q_o);
            q_c.push(gate.q_c);
        }

        SelectorPolynomials {
            q_l: inverse_fft(&q_l, omega),
            q_r: inverse_fft(&q_r, omega),
            q_m: inverse_fft(&q_m, omega),
            q_o: inverse_fft(&q_o, omega),
            q_c: inverse_fft(&q_c, omega),
        }
    }

    pub fn get_witness_polynomials(&self) -> WitnessPolynomials<F> {
        let omega = self.domain[1];
        WitnessPolynomials{
            a: inverse_fft(&self.witness.a, omega),
            b: inverse_fft(&self.witness.b, omega),
            c: inverse_fft(&self.witness.c, omega),
        }
    }
}