use crate::fft::{fft, inverse_fft, vec_to_poly};
use crate::gate::Gate;
use crate::witness::Witness;
use ark_ff::Field;
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;

pub struct Circuit<F: Field> {
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

impl<F: Field> Circuit<F> {
    pub fn new(
        gates: Vec<Gate<F>>,
        witness: Witness<F>,
        public_inputs: Vec<F>,
        domain: Vec<F>,
    ) -> Circuit<F> {
        Circuit {
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
        WitnessPolynomials {
            a: inverse_fft(&self.witness.a, omega),
            b: inverse_fft(&self.witness.b, omega),
            c: inverse_fft(&self.witness.c, omega),
        }
    }

    /// Checks if the constraint polynomial
    /// P(X) = QL(X)A(X) + QR(X)B(X) + QM(X)A(X)B(X) + QO(X)C(X) + QC(X)
    /// vanishes over evaluation domain H using pointwise operations
    pub fn is_gate_constraint_polynomial_zero_over_h(
        &self,
        selector: &SelectorPolynomials<F>,
        witness: &Witness<F>,
    ) -> bool {
        let omega = self.domain[1];
        let SelectorPolynomials {
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
        } = selector;
        let Witness { a, b, c } = witness;

        // getting polynomials in evaluation form
        let q_l = fft(q_l, omega);
        let q_r = fft(q_r, omega);
        let q_m = fft(q_m, omega);
        let q_o = fft(q_o, omega);
        let q_c = fft(q_c, omega);

        let a = fft(a, omega);
        let b = fft(b, omega);
        let c = fft(c, omega);

        let mut constraint_poly = Vec::with_capacity(a.len());
        for i in 0..a.len() {
            let term =
                q_l[i] * a[i] + q_r[i] * b[i] + q_m[i] * a[i] * b[i] + q_c[i] + q_o[i] * c[i];
            constraint_poly.push(term);
        }
        constraint_poly.iter().all(|x| x.is_zero())
    }

    fn get_gate_constraint_polynomial(
        &self,
        selector: &SelectorPolynomials<F>,
        witness: &Witness<F>,
    ) -> Vec<F> {
        let a_poly = vec_to_poly(witness.a.clone());
        let b_poly = vec_to_poly(witness.b.clone());
        let c_poly = vec_to_poly(witness.c.clone());

        let q_l_poly = vec_to_poly(selector.q_l.clone());
        let q_r_poly = vec_to_poly(selector.q_r.clone());
        let q_m_poly = vec_to_poly(selector.q_m.clone());
        let q_o_poly = vec_to_poly(selector.q_o.clone());
        let q_c_poly = vec_to_poly(selector.q_c.clone());

        // Now compute:
        let mut p_poly = q_l_poly.naive_mul(&a_poly);
        p_poly += &(q_r_poly.naive_mul(&b_poly));
        p_poly += &(q_m_poly.naive_mul(&a_poly.naive_mul(&b_poly)));
        p_poly += &(q_o_poly.naive_mul(&c_poly));
        p_poly += &q_c_poly;

        p_poly.coeffs.clone()
    }
}

#[cfg(test)]
mod tests {
    use crate::circuit::Circuit;
    use crate::gate::Gate;
    use crate::witness::Witness;
    use ark_bls12_381::Fr;
    use ark_ff::{FftField, Field};

    fn get_omega(n: usize) -> Fr {
        let generator = Fr::get_root_of_unity(n as u64).unwrap();
        generator
    }

    #[test]
    fn test_gate_constraint_is_zero_when_satisfied() {
        let n = 4;
        let omega = get_omega(n);
        let domain: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();

        // Construct simple addition gates: a + b = c
        let gates: Vec<Gate<Fr>> = vec![Gate::simple_addition_gate(); n];

        // Set witness such that a + b = c
        let witness = Witness {
            a: vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)],
            b: vec![Fr::from(9), Fr::from(8), Fr::from(7), Fr::from(6)],
            c: vec![Fr::from(10), Fr::from(10), Fr::from(10), Fr::from(10)],
        };

        let public_inputs = vec![];

        let circuit = Circuit::new(gates, witness.clone(), public_inputs, domain.clone());

        let selector = circuit.get_selector_polynomials();

        // Check that all evaluations are zero
        assert!(circuit.is_gate_constraint_polynomial_zero_over_h(&selector, &witness));
    }
}
