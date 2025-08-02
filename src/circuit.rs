use crate::fft::{fft, inverse_fft, pad_with_zeroes, vec_to_poly};
use crate::gate::Gate;
use crate::permutation::Permutation;
use crate::witness::Witness;
use ark_ff::Field;
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;
use itertools::Itertools;

/// A PLONK circuit containing gates, witness, public inputs, domain, and permutation
#[derive(Clone)]
pub struct Circuit<F: Field> {
    pub gates: Vec<Gate<F>>,
    pub witness: Witness<F>,
    pub public_inputs: Vec<F>,
    pub domain: Vec<F>,
    pub permutation: Permutation<F>,
}

/// Interpolated selector polynomials (q_L, q_R, q_M, q_O, q_C)
#[derive(Clone, Debug)]
pub struct SelectorPolynomials<F: Field> {
    pub q_l: DensePolynomial<F>,
    pub q_r: DensePolynomial<F>,
    pub q_m: DensePolynomial<F>,
    pub q_o: DensePolynomial<F>,
    pub q_c: DensePolynomial<F>,
}

/// Interpolated witness polynomials A(X), B(X), C(X)
#[derive(Clone, Debug)]
pub struct WitnessPolynomials<F: Field> {
    pub a: DensePolynomial<F>,
    pub b: DensePolynomial<F>,
    pub c: DensePolynomial<F>,
}

impl<F: Field> Circuit<F> {
    /// Constructs a new Circuit with gates, witness, and permutation wiring
    pub fn new(
        gates: Vec<Gate<F>>,
        witness: Witness<F>,
        public_inputs: Vec<F>,
        domain: Vec<F>,
        wiring: Vec<Vec<usize>>,
        k1: F,
        k2: F,
    ) -> Circuit<F> {
        let permutation = Permutation::new(witness.clone(), wiring, k1, k2);
        Circuit {
            gates,
            witness,
            public_inputs,
            domain,
            permutation,
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
            q_l: vec_to_poly(inverse_fft(&q_l, omega)),
            q_r: vec_to_poly(inverse_fft(&q_r, omega)),
            q_m: vec_to_poly(inverse_fft(&q_m, omega)),
            q_o: vec_to_poly(inverse_fft(&q_o, omega)),
            q_c: vec_to_poly(inverse_fft(&q_c, omega)),
        }
    }

    /// Interpolates witness wires A, B, C into DensePolynomials
    pub fn get_witness_polynomials(&self) -> WitnessPolynomials<F> {
        let omega = self.domain[1];
        WitnessPolynomials {
            a: vec_to_poly(inverse_fft(&self.witness.a, omega)),
            b: vec_to_poly(inverse_fft(&self.witness.b, omega)),
            c: vec_to_poly(inverse_fft(&self.witness.c, omega)),
        }
    }

    /// Checks if the constraint polynomial
    /// P(X) = QL(X)A(X) + QR(X)B(X) + QM(X)A(X)B(X) + QO(X)C(X) + QC(X)
    /// vanishes over evaluation domain H using point-wise operations
    pub fn is_gate_constraint_polynomial_zero_over_h(
        &self,
        selector: &SelectorPolynomials<F>,
        witness: &WitnessPolynomials<F>,
        public_input_polynomial: &DensePolynomial<F>,
    ) -> bool {
        let omega = self.domain[1];
        let SelectorPolynomials {
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
        } = selector;
        let WitnessPolynomials { a, b, c } = witness;

        // getting polynomials in evaluation form
        let q_l = fft(&pad_with_zeroes(&q_l, self.domain.len()), omega);
        let q_r = fft(&pad_with_zeroes(&q_r, self.domain.len()), omega);
        let q_m = fft(&pad_with_zeroes(&q_m, self.domain.len()), omega);
        let q_o = fft(&pad_with_zeroes(&q_o, self.domain.len()), omega);
        let q_c = fft(&pad_with_zeroes(&q_c, self.domain.len()), omega);

        let a = fft(a, omega);
        let b = fft(b, omega);
        let c = fft(c, omega);

        let pi = fft(
            &pad_with_zeroes(public_input_polynomial, self.domain.len()),
            omega,
        );

        let mut constraint_poly = Vec::with_capacity(a.len());
        for i in 0..a.len() {
            let term = q_l[i] * a[i]
                + q_r[i] * b[i]
                + q_m[i] * a[i] * b[i]
                + q_c[i]
                + q_o[i] * c[i]
                + pi[i];
            constraint_poly.push(term);
        }
        constraint_poly.iter().all(|x| x.is_zero())
    }

    /// Builds the full gate constraint polynomial symbolically
    pub fn generate_gate_constraint_polynomial(
        &self,
        selector: &SelectorPolynomials<F>,
        witness: &WitnessPolynomials<F>,
        public_input: DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let WitnessPolynomials { a, b, c } = witness;
        let SelectorPolynomials {
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
        } = selector;

        // Build each term safely
        let term_l = vec_to_poly(q_l.naive_mul(&a).coeffs);
        let term_r = vec_to_poly(q_r.naive_mul(&b).coeffs);
        let term_m = vec_to_poly(q_m.naive_mul(&vec_to_poly(a.naive_mul(&b).coeffs)).coeffs);
        let term_o = vec_to_poly(q_o.naive_mul(&c).coeffs);

        term_l + term_r + term_m + term_o + q_c.clone() + public_input
    }

    /// Shorthand for computing the constraint polynomial directly
    pub fn get_gate_constraint_polynomial(&self) -> DensePolynomial<F> {
        self.generate_gate_constraint_polynomial(
            &self.get_selector_polynomials(),
            &self.get_witness_polynomials(),
            self.compute_public_input_polynomial(),
        )
    }

    /// Constructs vanishing polynomial Z_H(X) = (X - ω^i)
    pub fn vanishing_poly(domain: &[F]) -> DensePolynomial<F> {
        let mut zh = DensePolynomial::from_coefficients_slice(&[F::one()]);
        for &root in domain {
            let x_minus_root = DensePolynomial::from_coefficients_slice(&[-root, F::one()]);
            zh = zh.naive_mul(&x_minus_root);
        }
        zh
    }

    /// Computes public input polynomial PI(X), where PI(ζ) = -Σ PIᵢ·Lᵢ(ζ)
    pub fn compute_public_input_polynomial(&self) -> DensePolynomial<F> {
        let mut evaluations = vec![F::zero(); self.domain.len()];
        for (i, &x) in self.public_inputs.iter().enumerate() {
            evaluations[i] = -x;
        }

        vec_to_poly(inverse_fft(&evaluations, self.domain[1]))
    }
}

#[cfg(test)]
mod tests {
    use crate::circuit::Circuit;
    use crate::gate::Gate;
    use crate::witness::Witness;
    use ark_bls12_381::Fr;
    use ark_ff::{FftField, Field, Zero};
    use ark_poly::univariate::DenseOrSparsePolynomial;
    use ark_poly::{EvaluationDomain, Polynomial};

    // helper functions
    fn get_omega(n: usize) -> Fr {
        let generator = Fr::get_root_of_unity(n as u64).unwrap();
        generator
    }
    fn fr(n: u64) -> Fr {
        Fr::from(n)
    }

    #[test]
    fn test_gate_constraint_is_zero_when_satisfied() {
        let n = 4;
        let omega = get_omega(n);
        let domain: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();

        // Construct simple addition gates: a + b = c
        let gates: Vec<Gate<Fr>> = vec![
            Gate::simple_addition_gate(),
            Gate::simple_addition_gate(),
            Gate::simple_addition_gate(),
            Gate::simple_mul_gate(),
        ];

        // Set witness such that a + b = c
        let witness_assignment = Witness {
            a: vec![Fr::from(1), Fr::from(3), Fr::from(3), Fr::from(5)],
            b: vec![Fr::from(9), Fr::from(8), Fr::from(6), Fr::from(5)],
            c: vec![Fr::from(10), Fr::from(11), Fr::from(9), Fr::from(25)],
        };

        let public_inputs = vec![];

        let circuit = Circuit::new(
            gates,
            witness_assignment.clone(),
            public_inputs,
            domain.clone(),
            Vec::new(),
            fr(2),
            fr(3),
        );

        // Selector and witness polynomials
        let selector = circuit.get_selector_polynomials();
        let witness = circuit.get_witness_polynomials();
        let pi_poly = circuit.compute_public_input_polynomial();

        // 1. Gate constraint polynomial should be zero over domain
        assert!(circuit.is_gate_constraint_polynomial_zero_over_h(&selector, &witness, &pi_poly));

        // 2. Extract P(x)
        let gate_poly = circuit.generate_gate_constraint_polynomial(&selector, &witness, pi_poly);
        assert!(gate_poly.coeffs.iter().filter(|&x| !x.is_zero()).count() > 0);

        // 3. P(x) should be divisible by Z_H(x)
        let zh = Circuit::vanishing_poly(&domain);
        let gate_dsp = DenseOrSparsePolynomial::from(gate_poly.clone());
        let zh_dsp = DenseOrSparsePolynomial::from(zh);
        let (_, remainder) = gate_dsp.divide_with_q_and_r(&zh_dsp).unwrap();
        assert!(remainder.is_zero());

        // 4. Evaluate P(x) over the domain to ensure all values are zero
        let evaluations_over_domain = domain
            .iter()
            .map(|x| gate_poly.evaluate(x))
            .collect::<Vec<Fr>>();
        assert!(evaluations_over_domain.iter().all(|x| x.is_zero()));
    }

    #[test]
    fn test_gate_constraint_simple_circuit() {
        let n = 2;
        let omega = get_omega(n);
        let domain: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();

        // Example: a + b = c  => q_l=1, q_r=1, q_m=0, q_o=-1, q_c=0
        // Constraint1: a + b - c = 0
        // Constraint2: ab - 13c = 0
        let gate1 = Gate::new(fr(1), fr(1), fr(0), -fr(1), fr(0)); // a + b - c = 0
        let gate2 = Gate::new(fr(0), fr(0), fr(1), -fr(13), fr(0)); // ab - 13c = 0
        let gates = vec![gate1, gate2];

        // Witness: a=3, b=4, c=7
        let witness = Witness {
            a: vec![fr(3), fr(13)],
            b: vec![fr(4), fr(2)],
            c: vec![fr(7), fr(2)],
        };

        let public_inputs = vec![];
        let circuit = Circuit::new(
            gates,
            witness,
            public_inputs,
            domain.clone(),
            vec![],
            fr(3),
            fr(3),
        );

        // Checking if CS is satisfied
        let selector = circuit.get_selector_polynomials();
        let witness = circuit.get_witness_polynomials();
        let pi_poly = circuit.compute_public_input_polynomial();

        assert!(circuit.is_gate_constraint_polynomial_zero_over_h(&selector, &witness, &pi_poly));

        // Checking if P(x) = 0 mod Zh(x)
        let gate_poly = circuit.generate_gate_constraint_polynomial(&selector, &witness, pi_poly);
        assert!(gate_poly.coeffs.iter().filter(|&x1| !x1.is_zero()).count() > 0);

        let zh = Circuit::vanishing_poly(&domain);
        let gate_dsp = DenseOrSparsePolynomial::from(gate_poly.clone());
        let zh_dsp = DenseOrSparsePolynomial::from(zh);
        let (t, remainder) = gate_dsp.divide_with_q_and_r(&zh_dsp).unwrap();
        assert!(remainder.is_zero());

        // Checking if P(x) is zero at w^i, i = [n]
        let evaluations_over_domain = domain
            .iter()
            .map(|x| gate_poly.evaluate(x))
            .collect::<Vec<Fr>>();
        assert!(evaluations_over_domain.iter().all(|x| x.is_zero()));
    }

    #[test]
    fn test_sum_and_product_with_two_public_inputs() {
        let n = 4;
        let omega = Fr::get_root_of_unity(n as u64).unwrap();
        let domain: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();

        let gate0 = Gate::new(fr(1), fr(0), fr(0), fr(0), fr(0)); // a = x₀
        let gate1 = Gate::new(fr(1), fr(0), fr(0), fr(0), fr(0)); // a = x₁
        let gate2 = Gate::new(fr(1), fr(1), fr(1), -fr(1), fr(0)); // a + b + ab = c
        let gate3 = Gate::new(fr(0), fr(0), fr(0), fr(0), fr(0)); // padding (optional)

        let gates = vec![gate0, gate1, gate2, gate3];

        // Witness
        let witness = Witness {
            a: vec![fr(3), fr(4), fr(3), fr(0)], // a[0]=x₀, a[1]=x₁, a[2]=x₀ for gate2
            b: vec![fr(0), fr(0), fr(4), fr(0)], // b[2]=x₁ for gate2
            c: vec![fr(0), fr(0), fr(19), fr(0)], // c[2] = 3 + 4 + 3*4 = 19
        };

        let public_inputs = vec![fr(3), fr(4)];

        let circuit = Circuit::new(
            gates,
            witness,
            public_inputs,
            domain.clone(),
            vec![vec![1, 6], vec![0, 2]],
            fr(3),
            fr(5),
        );

        let selector = circuit.get_selector_polynomials();
        let witness = circuit.get_witness_polynomials();
        let pi_poly = circuit.compute_public_input_polynomial();

        // Check constraints
        assert!(circuit.is_gate_constraint_polynomial_zero_over_h(&selector, &witness, &pi_poly));

        let gate_poly = circuit.generate_gate_constraint_polynomial(&selector, &witness, pi_poly);
        let zh = Circuit::vanishing_poly(&domain);
        let gate_dsp = DenseOrSparsePolynomial::from(gate_poly.clone());
        let zh_dsp = DenseOrSparsePolynomial::from(zh);
        let (_, remainder) = gate_dsp.divide_with_q_and_r(&zh_dsp).unwrap();
        assert!(remainder.is_zero());

        let evaluations = domain
            .iter()
            .map(|x| gate_poly.evaluate(x))
            .collect::<Vec<_>>();
        assert!(evaluations.iter().all(|x| x.is_zero()));
    }
}
