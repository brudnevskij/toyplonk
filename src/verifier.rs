use crate::circuit::Circuit;
use crate::fft::{inverse_fft, vec_to_poly};
use crate::prover::Proof;
use crate::transccript::hash_to_field;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, Zero};
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;

struct Challenges<E: Pairing> {
    alpha: E::ScalarField,
    beta: E::ScalarField,
    gamma: E::ScalarField,
    zeta: E::ScalarField,
    v: E::ScalarField,
    u: E::ScalarField,
}

impl<E: Pairing> Challenges<E> {
    pub fn new(proof: &Proof<E>) -> Self {
        let mut commitment_buffer = vec![proof.a, proof.b, proof.c];
        let beta = hash_to_field("beta", &commitment_buffer);
        let gamma = hash_to_field("gamma", &commitment_buffer);

        commitment_buffer.push(proof.z);
        let alpha = hash_to_field("alpha", &commitment_buffer);

        commitment_buffer.push(proof.t_lo);
        commitment_buffer.push(proof.t_mid);
        commitment_buffer.push(proof.t_hi);
        let zeta = hash_to_field("zeta", &commitment_buffer);
        let v = hash_to_field("v", &commitment_buffer);

        commitment_buffer.push(proof.w_zeta);
        commitment_buffer.push(proof.w_zeta_omega);
        let u = hash_to_field("u", &commitment_buffer);

        Challenges {
            alpha,
            beta,
            gamma,
            zeta,
            v,
            u,
        }
    }
}

pub struct VerifierPreprocessedInput<E: Pairing> {
    q_m: E::G1Affine,
    q_l: E::G1Affine,
    q_r: E::G1Affine,
    q_o: E::G1Affine,
    q_c: E::G1Affine,
    sigma_1: E::G1Affine,
    sigma_2: E::G1Affine,
    sigma_3: E::G1Affine,
    x: E::G2Affine,
}

pub fn verify_kzg_proof<E: Pairing>(
    proof: &Proof<E>,
    preprocessed_input: &VerifierPreprocessedInput<E>,
    public_input: &[E::ScalarField],
    domain: &[E::ScalarField],
) -> bool {
    // 1 - 3, validation is mostly done by type system

    // 4. challenges
    let Challenges {
        alpha,
        beta,
        gamma,
        zeta,
        v,
        u,
    } = Challenges::new(proof);

    // 5. zero polynomial eval
    let vanishing_polynomial = Circuit::vanishing_poly(domain);
    let vanishing_eval = vanishing_polynomial.evaluate(&zeta);

    // 6. Lagrange poly eval
    let l = compute_lagrange_polynomial_eval::<E>(domain.len(), domain[1], zeta);

    // 7. PI evaluation
    let pi = compute_public_input_polynomial::<E>(public_input, domain);
    let pi_eval = pi.evaluate(&zeta);

    true
}

fn compute_lagrange_polynomial_eval<E: Pairing>(
    n: usize,
    omega: E::ScalarField,
    zeta: E::ScalarField,
) -> E::ScalarField {
    // w(z^n - 1) / n(zeta-omega)
    let numerator = omega * (zeta.pow(&[n as u64]) - E::ScalarField::one());
    let denominator = E::ScalarField::from(n as u64) * (zeta - omega);
    numerator / denominator
}

pub fn compute_public_input_polynomial<E: Pairing>(
    public_inputs: &[E::ScalarField],
    domain: &[E::ScalarField],
) -> DensePolynomial<E::ScalarField> {
    let mut evaluations = vec![E::ScalarField::zero(); domain.len()];
    for (i, &x) in public_inputs.iter().enumerate() {
        evaluations[i] = -x;
    }

    vec_to_poly(inverse_fft(&evaluations, domain[1]))
}
